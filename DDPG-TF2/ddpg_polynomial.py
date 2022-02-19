"""
Implementation of DDPG-SOS program

The algorithm is tested on the Pendulum-v1 and Polynomial-v0 OpenAI gym task
and developed with SOSOPT + Tensorflow

Author: Hejun Huang
"""

import os
import tensorflow as tf
import numpy as np
import pandas as pd
import gym
import argparse
from scipy.io import savemat
import datetime
from absl import app
from replay_buffer import ReplayBuffer
from matplotlib import pyplot as plt
from sosp_cbf_learner import LEARNER
import matlab
import matlab.engine
import sys
from gym import spaces
from barrier_previous import BARRIER

sys.path.append(os.getcwd())
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "3"


# Define the input parameters
def parse_args():
    parser = argparse.ArgumentParser(description='provide arguments for DDPG agent')
    # agent parameters
    parser.add_argument('--actor-lr', help='actor network learning rate', default=0.0001)
    parser.add_argument('--critic-lr', help='critic network learning rate', default=0.0005)
    parser.add_argument('--mlp-lr', help='mlp model learning rate', default=0.0001)
    parser.add_argument('--gamma', help='discount factor for critic updates', default=0.99)
    parser.add_argument('--tau', help='soft target update parameter', default=0.001)
    parser.add_argument('--buffer-size', help='max size of the replay buffer', default=1000000)
    parser.add_argument('--minibatch-size', help='size of minibatch for minibatch-SGD', default=128)
    parser.add_argument('--env', help='choose the gym env- tested on Single Pendulum',
                        default='Pendulum-v1')  # Pendulum-v0, Pendulum-v1, Polynomial-v0
    parser.add_argument('--random-seed', help='random seed for repeatability', default=1754)
    parser.add_argument('--max-episodes', help='max num of episodes to do while training', default=50)  # 600
    parser.add_argument('--max-episode-len', help='max length of 1 episode', default=200)  # 200
    parser.add_argument('--render-env', help='render the gym env', action='store_false') # store_false
    parser.add_argument('--use-gym-monitor', help='record gym results', action='store_false') # store_false
    parser.add_argument('--monitor-dir', help='directory for storing gym results', default='./res/gym_ddpg')
    parser.add_argument('--summary-dir', help='directory for storing tensorboard info', default='./res/tf_ddpg')
    parser.add_argument('--mat-dir', help='directory for storing mat info', default='./res/mat')
    parser.add_argument('--csv-dir', help='directory for storing csv info', default='./res/csv')
    parser.add_argument('--mlp-mode', help='MLP model training model', default='original')  # original OR revised
    parser.add_argument('--mat-path', help='GP+SOSP store data', default='./res/mat/obs_env.mat')  # original OR revised
    parser.add_argument('--mat-path_bak', help='GP+SOSP store data',
                        default='./res/mat/obs_env_record.mat')  # original OR revised
    parser.set_defaults(render_env=True) # False
    parser.set_defaults(use_gym_monitor=True) # False
    args = vars(parser.parse_args())
    return args


def build_barrier(self, eng):
    # eng.cd('.\matlab_code', nargout=0)
    aa = self.firstIter
    c, d, e, f = eng.main_function(aa, nargout=4)
    self.a_bar = np.array(c)
    return c, d, e, f


def poly_ubar_match(b, s_input):
    x1 = s_input[0, 0]
    x2 = s_input[0, 1]
    fac = np.array([[x1 ** 2, x1 * x2, x1, x2 ** 2, x2, 1]])
    fun = np.dot(fac, np.array(b).T).reshape([1, ])
    return fun


def validate_safe_state(x1, x2, sys_bar, region, sys_d2, eng):
    # eng.cd('.\matlab_code', nargout=0)
    c = eng.check_state(x1, x2, sys_bar, region, sys_d2, nargout=1)
    return c


# ===========================
#   Actor and Critic DNNs
# ===========================

class ActorNetwork(object):
    """
    Input to the network is the state, output is the action
    under a deterministic policy.

    The output layer activation is a tanh to keep the action
    between -action_bound and action_bound
    """

    def __init__(self, state_dim, action_dim, action_bound, learning_rate, tau, batch_size):
        self.s_dim = state_dim
        self.a_dim = action_dim
        self.action_bound = action_bound  # tf.constant(action_bound)
        self.learning_rate = learning_rate
        self.tau = tau
        self.batch_size = batch_size

        # Actor Network
        self.actor_model = self.create_actor_network()
        # self.actor_model.summary()

        # Target Network
        self.target_actor_model = self.create_actor_network()

        # Initialize the optimizer of actor network
        self.actor_opt = tf.keras.optimizers.Adam(self.learning_rate)
        # Get the number of actor network and target actor network.
        self.num_trainable_vars = len(self.actor_model.trainable_variables) + len(
            self.target_actor_model.trainable_variables)

    def create_actor_network(self):
        inputs = tf.keras.Input(shape=(self.s_dim,))
        net = tf.keras.layers.Dense(units=256, use_bias=True)(inputs)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        net = tf.keras.layers.Dense(units=64, use_bias=True)(net)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        # Final layer weights are init to Uniform[-3e-3, 3e-3]
        w_init = tf.keras.initializers.RandomUniform(minval=-0.003, maxval=0.003)
        out = tf.keras.layers.Dense(units=self.a_dim, activation='tanh', use_bias=True, kernel_initializer=w_init)(net)
        # Scale output to -action_bound to action_bound
        scaled_out = tf.keras.layers.Multiply()([out, self.action_bound])
        actor_model = tf.keras.models.Model(inputs=inputs, outputs=scaled_out)
        return actor_model

    def train(self, inputs, a_gradient):
        # inputs: state
        # a_gradient: the graident of critic model loss
        with tf.GradientTape() as tape:
            unnormal_grads = tape.gradient(self.actor_model(inputs, training=True),
                                           self.actor_model.trainable_variables, -a_gradient)
            grads = list(map(lambda x: tf.divide(x, self.batch_size), unnormal_grads))
        self.actor_opt.apply_gradients(zip(grads, self.actor_model.trainable_variables))

    def predict(self, inputs):
        # actor_model
        # inputs: states
        return self.actor_model.predict(inputs)

    def predict_target(self, inputs):
        # actor_model
        # inputs: states
        return self.target_actor_model.predict(inputs)

    def update_target_network(self):
        actor_model_w = self.actor_model.get_weights()
        target_actor_model_w = self.target_actor_model.get_weights()
        update_w = []
        for src, dest in zip(actor_model_w, target_actor_model_w):
            u_w = (src * self.tau) + (dest * (1 - self.tau))
            update_w.append(u_w)
        self.target_actor_model.set_weights(update_w)

    def get_num_trainable_vars(self):
        return self.num_trainable_vars


class CriticNetwork(object):
    """
    Input to the network is the state and action, output is Q(s,a).
    The action must be obtained from the output of the Actor network.

    """

    def __init__(self, state_dim, action_dim, learning_rate, tau, gamma, num_actor_vars):
        self.s_dim = state_dim
        self.a_dim = action_dim
        self.learning_rate = learning_rate
        self.tau = tau
        self.gamma = gamma

        # Create the critic network
        self.critic_model = self.create_critic_network()
        # self.critic_model.summary()

        # Target Network
        self.target_critic_model = self.create_critic_network()

        # Define loss and optimization Op
        self.critic_opt = tf.keras.optimizers.Adam(self.learning_rate)

        # Get the gradient of the net w.r.t. the action.
        # For each action in the minibatch (i.e., for each x in xs),
        # this will sum up the gradients of each critic output in the minibatch
        # w.r.t. that action. Each output is independent of all
        # actions except for one.
        # self.action_grads = tf.gradients(self.out, self.action)

    def create_critic_network(self):
        inputs = tf.keras.Input(shape=(self.s_dim,))
        action = tf.keras.Input(shape=(self.a_dim,))

        net = tf.keras.layers.Dense(units=256, use_bias=True)(inputs)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        t1 = tf.keras.layers.Dense(units=64, use_bias=False)(net)
        t2 = tf.keras.layers.Dense(units=64, use_bias=True)(action)
        out = tf.keras.layers.Add()([t1, t2])
        out = tf.keras.layers.ReLU()(out)

        # linear layer connected to 1 output representing Q(s,a)
        # Weights are init to Uniform[-3e-3, 3e-3]
        w_init = tf.keras.initializers.RandomUniform(minval=-0.003, maxval=0.003)
        out = tf.keras.layers.Dense(units=1, use_bias=True, kernel_initializer=w_init)(out)
        critic_model = tf.keras.models.Model(inputs=[inputs, action], outputs=out)

        return critic_model

    def compute_loss(self, v_pred, predicted_q_value):
        mse = tf.keras.losses.MeanSquaredError()
        return mse(predicted_q_value, v_pred)

    def train(self, inputs, action, predicted_q_value):
        with tf.GradientTape() as tape:
            v_pred = self.critic_model([inputs, action], training=True)
            assert v_pred.shape == predicted_q_value.shape

            loss = self.compute_loss(v_pred, tf.stop_gradient(predicted_q_value))

        grads = tape.gradient(loss, self.critic_model.trainable_variables)
        self.critic_opt.apply_gradients(zip(grads, self.critic_model.trainable_variables))
        return v_pred

    def predict(self, inputs, action):
        return self.critic_model.predict([inputs, action])

    def predict_target(self, inputs, action):
        return self.target_critic_model.predict([inputs, action])

    def action_gradients(self, inputs, actions):
        actions = tf.convert_to_tensor(actions)
        with tf.GradientTape() as tape:
            tape.watch(actions)
            q_values_x = self.critic_model([inputs, actions])
            q_values = tf.squeeze(q_values_x)  # The same gradient can be obtained if remove it.

        return tape.gradient(q_values, actions)

    def update_target_network(self):
        critic_model_w = self.critic_model.get_weights()
        target_critic_model_w = self.target_critic_model.get_weights()
        update_w = []
        for src, dest in zip(critic_model_w, target_critic_model_w):
            u_w = (src * self.tau) + (dest * (1 - self.tau))
            update_w.append(u_w)
        self.target_critic_model.set_weights(update_w)


class OrnsteinUhlenbeckActionNoise:
    def __init__(self, mu, sigma=0.3, theta=0.1, dt=1e-2, x0=None):
        self.theta = theta
        self.mu = mu
        self.sigma = sigma
        self.dt = dt
        self.x0 = x0
        self.reset()

    def __call__(self):
        x = self.x_prev + self.theta * (self.mu - self.x_prev) * self.dt + \
            self.sigma * np.sqrt(self.dt) * np.random.normal(size=self.mu.shape)
        self.x_prev = x
        return x

    def reset(self):
        self.x_prev = self.x0 if self.x0 is not None else np.zeros_like(self.mu)

    def __repr__(self):
        return 'OrnsteinUhlenbeckActionNoise(mu={}, sigma={})'.format(self.mu, self.sigma)


# ===========================
#   Tensorflow Summary Ops
# ===========================

# def build_summaries():
#     episode_reward = tf.Variable(0.)
#     tf.summary.scalar("Reward", episode_reward)
#     episode_ave_max_q = tf.Variable(0.)
#     tf.summary.scalar("Qmax Value", episode_ave_max_q)

#     summary_vars = [episode_reward, episode_ave_max_q]
#     summary_ops = tf.summary.merge_all()

#     return summary_ops, summary_vars

# ===========================
#   Agent Training
# ===========================

def train(env, args, actor, critic, actor_noise, agent, eng):
    # Set up summary Ops
    sub_episodes = 6  # 20
    writer = tf.summary.create_file_writer(args['summary_dir'])
    reward_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    maxq_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    max_angle_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500

    # Initialize target network weights
    actor.update_target_network()
    critic.update_target_network()

    # Initialize replay memory
    replay_buffer = ReplayBuffer(int(args['buffer_size']), int(args['random_seed']))
    paths = list()

    s_record_sosp = []  # np.empty(shape=(0, 3))
    action_record_sosp = []  # np.array([])
    reward_boundary = []
    sample_len = 3

    for i in range(int(args['max_episodes'])):

        if (agent.firstIter == 1):
            obs, information, action_bar = [], [], []
            for j in range(int(args['max_episode_len']) * sample_len):
                s = env.reset()
                while (env.unwrapped.state[0] > 0.8 or env.unwrapped.state[0] < -0.8):
                    s = env.reset()
                    record = env.unwrapped.state
                a = actor.predict(np.reshape(s, (1, actor.s_dim))) + actor_noise()
                obs.append(s)
                information.append(record)
                action_bar.append(a)

            mat_path_bak = args['mat_path_bak']
            savemat(mat_path_bak, dict(data_base=information, action_base=action_bar))
            coe, sys_bar, region, sys_d2 = build_barrier(agent, eng)
        else:
            coe, sys_bar, region, sys_d2 = build_barrier(agent, eng)

        for el in range(sub_episodes):  # TODO 2
            # print('test369')
            obs, action, rewards, action_bar, action_BAR, action_RL_list, information = [], [], [], [], [], [], []
            # obs, action, rewards = [], [], []
            # print("Epoch: {}".format(i))
            s = env.reset()
            a_obs_raw = env.unwrapped.state

            kk = a_obs_raw.tolist()
            sta_1 = kk[0]
            sta_2 = kk[1]
            i = 0
            while i:
                i = validate_safe_state(sta_1, sta_2, sys_bar, region, sys_d2, eng)
                if i != 1:
                    s = env.reset()
                    a_obs_raw = env.unwrapped.state
                print(i)
                kk = a_obs_raw.tolist()
                sta_1 = kk[0]
                sta_2 = kk[1]

            ep_reward = 0
            ep_ave_max_q = 0
            r_lst = []

            for j in range(int(args['max_episode_len'])):
                s_cur_ = np.array([[a_obs_raw[0], a_obs_raw[1]]])
                # env.render()
                # Added exploration noise
                # a = actor.predict(np.reshape(s, (1, 3))) + (1. / (1. + i))
                a = actor.predict(np.reshape(s, (1, actor.s_dim))) + actor_noise()

                # Incorporate barrier function
                action_rl = a[0]  # (1, )  [[3.44]] => [3.44]

                # Utilize compensation barrier function
                if (agent.firstIter == 1):
                    u_BAR_ = [0]  # 在第一个ep时的取值，因此时未对agent.bar_comp的MLP训练
                else:
                    u_BAR_ = agent.bar_comp.get_action(s)[0]  # (1, )

                action_RL = action_rl + u_BAR_  # (1, )
                u_bar_ = poly_ubar_match(coe, s_cur_)  # TODO 3
                u_bar_ = u_bar_ / 10.
                # Test 1
                action_ = (action_RL + u_bar_) / 2.
                # Test 2
                # action_ = u_bar_ - action_rl + u_BAR_

                # Set threshold value
                # thres_hold = 10;
                # for i in range(len(u_bar_)):
                #     if u_bar_[i-1] >= thres_hold:
                #         u_bar_[i-1] = thres_hold
                #     elif u_bar_[i-1] <= -thres_hold:
                #         u_bar_[i - 1] = -thres_hold

                # Define Reward
                s2, r, terminal, info = env.step(a[0])

                # Define Reward
                # s2, r_cur, terminal, info = env.step(action_)
                # float_u_bar_ = u_bar_.astype(np.float64)
                # print(float_u_bar_[0])
                # r = r_cur - float_u_bar_[0]

                replay_buffer.add(np.reshape(s, (actor.s_dim,)), np.reshape(a, (actor.a_dim,)), r,
                                  terminal, np.reshape(s2, (actor.s_dim,)))

                if replay_buffer.size() > int(args['minibatch_size']):
                    s_batch, a_batch, r_batch, t_batch, s2_batch = \
                        replay_buffer.sample_batch(int(args['minibatch_size']))

                    # Calculate targets
                    target_q = critic.predict_target(
                        s2_batch, actor.predict_target(s2_batch))

                    y_i = []

                    for k in range(int(args['minibatch_size'])):
                        if t_batch[k]:
                            y_i.append(np.array([r_batch[k]]))
                        else:
                            y_i.append(r_batch[k] + critic.gamma * target_q[k])

                    predicted_q_value = critic.train(s_batch, a_batch,
                                                     np.reshape(y_i, (int(args['minibatch_size']), 1)))

                    ep_ave_max_q += np.amax(predicted_q_value)

                    # Update the actor policy using the sampled gradient
                    a_outs = actor.predict(s_batch)
                    s_grads = critic.action_gradients(s_batch, a_outs)

                    grads = np.array(s_grads).reshape((-1, actor.a_dim))
                    actor.train(s_batch, grads)

                    # Update target networks
                    actor.update_target_network()
                    critic.update_target_network()

                s = s2
                ep_reward += r
                r_lst.append(r)

                obs.append(s)
                information.append(a_obs_raw)
                rewards.append(r)
                action_bar.append(u_bar_)
                action_BAR.append(u_BAR_)
                action.append(action_)
                action_RL_list.append(action_rl)

                if terminal:
                    obs_x = np.concatenate(information).reshape((200, 2))
                    max_theta = np.max(obs_x[:, 0])
                    with writer.as_default():
                        print(action[len(action) - 1], '|', action_bar[len(action_bar) - 1], '|',
                              action_RL_list[len(action_RL_list) - 1], '|', action_BAR[len(action_BAR) - 1])
                        tf.summary.scalar("Reward", ep_reward, step=i * sub_episodes + el)
                        tf.summary.scalar("Qmax Value", ep_ave_max_q / float(j), step=i * sub_episodes + el)
                        tf.summary.scalar("Max Angle", max_theta, step=i * sub_episodes + el)
                    writer.flush()

                    print('| Reward: {:d} | Episode: {:d} | Qmax: {:.4f}'.format(int(ep_reward), \
                                                                                 i, (ep_ave_max_q / float(j))))
                    reward_result[i * sub_episodes + el] = ep_reward
                    maxq_result[i * sub_episodes + el] = ep_ave_max_q / float(j)
                    max_angle_result[i * sub_episodes + el] = max_theta

                    path = {"Observation": np.concatenate(obs).reshape((200, 3)),
                            "Action": np.concatenate(action),
                            "Action_bar": np.concatenate(action_bar),
                            "Action_BAR": np.concatenate(action_BAR),
                            "Reward": np.asarray(rewards)}
                    paths.append(path)

                    if len(action_record_sosp) < sample_len:
                        cur_s_reshape = np.concatenate(information).reshape((200, 2))
                        s_record_sosp.append(cur_s_reshape)
                        cur_action_reshape = np.concatenate(action).reshape((200, 1))
                        action_record_sosp.append(cur_action_reshape)
                        reward_boundary.append(ep_reward)
                        mat_path = args['mat_path']

                        if len(reward_boundary) == sample_len:
                            save_s_sosp = np.reshape(s_record_sosp, (200 * sample_len, 2))
                            save_action_sosp = np.reshape(action_record_sosp, (200 * sample_len, 1))
                            savemat(mat_path, dict(data=save_s_sosp, action=save_action_sosp))
                    else:
                        mat_path = args['mat_path']
                        min_boundary = np.min(reward_boundary)
                        if ep_reward > min_boundary:
                            min_idx = np.argmin(reward_boundary)
                            reward_boundary[min_idx] = ep_reward
                            cur_s_reshape = np.concatenate(information).reshape((200, 2))
                            s_record_sosp[min_idx] = cur_s_reshape
                            cur_action_reshape = np.concatenate(action).reshape((200, 1))
                            action_record_sosp[min_idx] = cur_action_reshape

                            save_s_sosp = np.reshape(s_record_sosp, (200 * sample_len, 2))
                            save_action_sosp = np.reshape(action_record_sosp, (200 * sample_len, 1))
                            savemat(mat_path, dict(data=save_s_sosp, action=save_action_sosp))
                    """
                    if sample_cnt >= (sub_episodes - sample_len):
                        cur_s_reshape = np.concatenate(obs).reshape((200, 3))
                        s_record_sosp = np.append(s_record_sosp, cur_s_reshape, axis=0) # (400, 3)
                        cur_action_reshape = np.concatenate(action)
                        action_record_sosp = np.append(action_record_sosp, cur_action_reshape, axis=0)

                        if sample_cnt == (sub_episodes - 1):
                            action_record_sosp = np.reshape(action_record_sosp, (200*sample_len, 1))

                            # mat_path = os.path.join(args['mat_dir'],
                            #                         "obs_env.mat".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M")))
                            savemat(mat_path, dict(data=s_record_sosp, action=action_record_sosp))

                    sample_cnt = sample_cnt + 1
                    """
                    cur_iter = i * sub_episodes + el + 1
                    os.makedirs(args['csv_dir'], exist_ok=True)
                    rec_pd = pd.DataFrame({'Reward': reward_result[:cur_iter],
                                           'Qmax': maxq_result[:cur_iter],
                                           'AngleMax': max_angle_result[:cur_iter]},
                                          columns=['Reward', 'Qmax', 'AngleMax'])
                    rec_pd.to_csv(os.path.join(args['csv_dir'], 'log.csv'), index=False)

                    # Trend graph
                    fig = plt.figure()
                    plt.plot(np.arange(len(maxq_result[:cur_iter])), maxq_result[:cur_iter])
                    plt.xlabel('Episode')
                    plt.ylabel('Max Q value')
                    pic_name = "Convergence of Max Q Value"
                    plt.title(pic_name)
                    img_path = os.path.join(args['csv_dir'], "MaxQ.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    fig = plt.figure()
                    plt.plot(np.arange(len(reward_result[:cur_iter])), reward_result[:cur_iter], label='ddpg-sosp',
                             color='r')
                    plt.xlabel('Episode')
                    plt.ylabel('Episode Reward')
                    pic_name = "Reward vs. Episode"
                    plt.title(pic_name)
                    plt.legend()
                    img_path = os.path.join(args['csv_dir'], "Reward.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    fig = plt.figure()
                    plt.plot(np.arange(len(max_angle_result[:cur_iter])), np.ones(len(max_angle_result[:cur_iter])),
                             label='Safe Boundary', color='k', linestyle='--')
                    plt.plot(np.arange(len(max_angle_result[:cur_iter])), max_angle_result[:cur_iter],
                             label='ddpg-sosp', color='r')
                    plt.xlabel('Episode')
                    plt.ylabel('Max Angle (rad)')
                    pic_name = "Safety Violation"
                    plt.title(pic_name)
                    plt.legend()
                    img_path = os.path.join(args['csv_dir'], "MaxAngle.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    break

        plt.close('all')
        agent.firstIter = 0  # 在第一个ep后变为0
        agent.bar_comp.get_training_rollouts(paths)
        barr_loss = agent.bar_comp.train()

    return [paths, reward_result]


def main(_argv):
    # gpu_devices = tf.config.experimental.list_physical_devices('GPU')  # Search available GPU
    # tf.debugging.set_log_device_placement(True)
    # if len(gpu_devices) > 0:
    #     # os.environ['CUDA_VISIBLE_DEVICES'] = "0"
    #     tf.config.set_visible_devices(devices=gpu_devices, device_type='GPU')  # Select GPU
    #     for single_gpu in gpu_devices:
    #         tf.config.experimental.set_memory_growth(single_gpu, True)  # Set the graphics memory usage
    # else:
    #     print("Cannot find the GPU device.")
    #     os.system("pause")

    # Setup matlab engine
    matlab_eng_id = matlab.engine.find_matlab()
    eng = matlab.engine.connect_matlab(matlab_eng_id[0])

    args = parse_args()

    env = gym.make(args['env'])
    np.random.seed(int(args['random_seed']))
    tf.random.set_seed(int(args['random_seed']))
    env.seed(int(args['random_seed']))

    # Set environment parameters for pendulum
    env.unwrapped.max_torque = 15.
    env.unwrapped.max_speed = 60.
    env.unwrapped.action_space = spaces.Box(low=-env.unwrapped.max_torque, high=env.unwrapped.max_torque, shape=(1,))
    high = np.array([1., 1., env.unwrapped.max_speed])
    env.unwrapped.observation_space = spaces.Box(low=-high, high=high)

    state_dim = env.observation_space.shape[0]
    action_dim = env.action_space.shape[0]
    action_bound = env.action_space.high

    # Ensure action bound is symmetric
    assert (env.action_space.high == -env.action_space.low)

    actor = ActorNetwork(state_dim, action_dim, action_bound,
                         float(args['actor_lr']), float(args['tau']),
                         int(args['minibatch_size']))

    critic = CriticNetwork(state_dim, action_dim,
                           float(args['critic_lr']), float(args['tau']),
                           float(args['gamma']),
                           actor.get_num_trainable_vars())

    actor_noise = OrnsteinUhlenbeckActionNoise(mu=np.zeros(action_dim))

    # Setup CBF Agent
    agent = LEARNER(env, eng)
    agent.max_episode_len = int(args['max_episode_len'])
    agent.obs_file_path = './res/mat/obs_save.mat'
    agent.file_path = 'obs.mat'
    agent.firstIter = 1

    agent.bar_comp = BARRIER(state_dim, action_dim, float(args['mlp_lr']),
                             args['mlp_mode'])  # 3: state_dim, 1: action_dim

    # Train
    [paths, reward_result] = train(env, args, actor, critic, actor_noise, agent, eng)

    # Save Result
    os.makedirs(args['mat_dir'], exist_ok=True)
    mat_path_test = os.path.join(args['mat_dir'],
                                 "data4_{}.mat".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M")))
    savemat(mat_path_test, dict(data=paths, reward=reward_result))
    return 0


if __name__ == '__main__':

    try:
        app.run(main)
    except SystemExit:
        pass
