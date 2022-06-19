"""
Implementation of DDPG-CBF on the Pendulum-v1 OpenAI gym task

Author: rcheng805
"""

import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
import numpy as np

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
import pandas as pd
import gym
from gym import wrappers
# import tflearn
import argparse
import pprint as pp
from scipy.io import savemat
import random
from absl import app
from matplotlib import pyplot as plt
from replay_buffer import ReplayBuffer

from learner import LEARNER
from barrier_comp import BARRIER
import cbf
import dynamics_gp
import datetime
from gym import spaces
import sys

sys.path.append('../matlab_demo/matlab_rl/matlab_code')
import logging

logging.getLogger("tensorflow").setLevel(logging.WARNING)
import matlab
import matlab.engine


def parse_args():
    cur_time = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    parser = argparse.ArgumentParser(description='provide arguments for DDPG agent')

    # agent parameters
    parser.add_argument('--actor-lr', help='actor network learning rate', default=0.0001)
    parser.add_argument('--critic-lr', help='critic network learning rate', default=0.001)
    parser.add_argument('--mlp-lr', help='mlp model learning rate', default=0.0001)
    parser.add_argument('--gamma', help='discount factor for critic updates', default=0.99)
    parser.add_argument('--tau', help='soft target update parameter', default=0.001)
    parser.add_argument('--buffer-size', help='max size of the replay buffer', default=1000000)
    parser.add_argument('--minibatch-size', help='size of minibatch for minibatch-SGD', default=64)

    # run parameters
    parser.add_argument('--env', help='choose the gym env- tested on {Pendulum-v0}',
                        default='Pendulum-v1')  # Pendulum-v0
    parser.add_argument('--random-seed', help='random seed for repeatability', default=1234)
    parser.add_argument('--max-episodes', help='max num of episodes to do while training', default=150)  # 150
    parser.add_argument('--max-episode-len', help='max length of 1 episode', default=200)  #
    parser.add_argument('--render-env', help='render the gym env', action='store_false')
    parser.add_argument('--use-gym-monitor', help='record gym results', action='store_false')
    parser.add_argument('--monitor-dir', help='directory for storing gym results', default='./res/QP/gym_ddpg/{}'.format(cur_time))
    parser.add_argument('--summary-dir', help='directory for storing tensorboard info', default='./res/QP/tensorboard/{}'.format(cur_time))
    parser.add_argument('--mat-dir', help='directory for storing mat info', default='./res/mat')
    parser.add_argument('--csv-dir', help='directory for storing csv info', default='./res/QP/csv/{}'.format(cur_time))
    parser.add_argument('--mlp-mode', help='MLP model training model', default='original')  # original OR revised
    parser.add_argument('--method', help='QP/SOSP/ONLY(ddpg)', default='QP')

    parser.set_defaults(render_env=False)
    parser.set_defaults(use_gym_monitor=False)

    args = vars(parser.parse_args())
    pp.pprint(args)

    return args


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
        self.action_bound = action_bound
        self.learning_rate = learning_rate
        self.tau = tau
        self.batch_size = batch_size

        # Actor Network
        # self.inputs, self.out, self.scaled_out = self.create_actor_network()
        self.actor_model = self.create_actor_network()
        self.actor_model.summary()
        # model.trainable_variables, model.get_weights()
        # self.network_params = tf.trainable_variables()

        # Target Network
        self.target_actor_model = self.create_actor_network()

        # Initialize the optimizer of actor network
        self.actor_opt = tf.keras.optimizers.Adam(self.learning_rate)
        # Get the number of actor network and target actor network.
        self.num_trainable_vars = len(self.actor_model.trainable_variables) + len(self.target_actor_model.trainable_variables)

    def create_actor_network(self):
        """
        Create actor network
        :param input_state_shape: state
        :return: act
        """
        inputs = tf.keras.Input(shape=(self.s_dim,))
        net = tf.keras.layers.Dense(units=400, use_bias=True)(inputs)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        net = tf.keras.layers.Dense(units=300, use_bias=True)(net)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        # Final layer weights are init to Uniform[-3e-3, 3e-3]
        w_init = tf.keras.initializers.RandomUniform(minval=-0.003, maxval=0.003)
        out = tf.keras.layers.Dense(units=self.a_dim, activation='tanh', use_bias=True, kernel_initializer=w_init)(net)
        # Scale output to -action_bound to action_bound (-15, 15)
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
        self.critic_model.summary()
        # model.trainable_variables, model.get_weights()
        # self.network_params = tf.trainable_variables()[num_actor_vars:]

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

        net = tf.keras.layers.Dense(units=400, use_bias=True)(inputs)
        net = tf.keras.layers.BatchNormalization()(net)
        net = tf.keras.layers.ReLU()(net)
        t1 = tf.keras.layers.Dense(units=300, use_bias=False)(net)
        t2 = tf.keras.layers.Dense(units=300, use_bias=True)(action)
        out = tf.keras.layers.Add()([t1, t2])
        out = tf.keras.layers.ReLU()(out)

        # linear layer connected to 1 output representing Q(s,a)
        # Weights are init to Uniform[-3e-3, 3e-3]
        w_init = tf.keras.initializers.RandomUniform(minval=-0.003, maxval=0.003)
        out = tf.keras.layers.Dense(units=1, use_bias=True, kernel_initializer=w_init)(out)
        critic_model = tf.keras.models.Model(inputs=[inputs, action], outputs=out)
        # (1, ) => (1batch, dim1)
        # print(np.shape(out))
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
        return v_pred  # Q-value

    def predict(self, inputs, action):
        # inputs: states
        # action
        return self.critic_model.predict([inputs, action])

    def predict_target(self, inputs, action):
        return self.target_critic_model.predict([inputs, action])

    def action_gradients(self, inputs, actions):
        actions = tf.convert_to_tensor(actions)
        with tf.GradientTape() as tape:
            tape.watch(actions)
            q_values = self.critic_model([inputs, actions])
            q_values = tf.squeeze(q_values)
        return tape.gradient(q_values, actions)

    def update_target_network(self):
        critic_model_w = self.critic_model.get_weights()
        target_critic_model_w = self.target_critic_model.get_weights()
        update_w = []
        for src, dest in zip(critic_model_w, target_critic_model_w):
            u_w = (src * self.tau) + (dest * (1 - self.tau))
            update_w.append(u_w)
        self.target_critic_model.set_weights(update_w)


# Taken from https://github.com/openai/baselines/blob/master/baselines/ddpg/noise.py, which is
# based on http://math.stackexchange.com/questions/1287634/implementing-ornstein-uhlenbeck-in-matlab
class OrnsteinUhlenbeckActionNoise:
    def __init__(self, mu, sigma=0.3, theta=.15, dt=1e-2, x0=None):
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
#     # TODO
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

def train(env: object, args, actor, critic, actor_noise, agent, eng) -> object:
    # Set up summary Ops
    # summary_ops, summary_vars = build_summaries()
    sub_episodes = 5
    writer = tf.summary.create_file_writer(args['summary_dir'])
    reward_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    maxq_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    max_angle_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    max_angleSpeed_result = np.zeros(int(args['max_episodes']) * sub_episodes)  # Raw code: 2500
    
     # sess.run(tf.global_variables_initializer())

    # Initialize target network weights
    actor.update_target_network()
    critic.update_target_network()

    # Initialize replay memory
    replay_buffer = ReplayBuffer(int(args['buffer_size']), int(args['random_seed']))

    # Needed to enable BatchNorm.
    # This hurts the performance on Pendulum but could be useful
    # in other environments.
    # tflearn.is_training(True)

    paths = list()

    for i in range(int(args['max_episodes'])):

        # Utilize GP from previous iteration while training current iteration
        if (agent.firstIter == 1):
            pass
        else:
            # 在第一个ep后才会将GP_model赋值与GP_model_prev
            agent.GP_model_prev = agent.GP_model.copy()
            dynamics_gp.build_GP_model(agent)

        for el in range(sub_episodes):
            # 每次结束后都会Update一次GP model(除最后一次), 但本次episodes仍使用上一ep训练的GP，不会采用本轮最新实时update的GP.
            obs, action, rewards, action_bar, action_BAR, action_RL_list = [], [], [], [], [], []
            unwrapped_state_lst = []

            s = env.reset()
            # a_temp = env.unwrapped.state
            # print(env.unwrapped.state)
            # Ensure that starting position is in "safe" region
            while (env.unwrapped.state[0] > 0.8 or env.unwrapped.state[0] < -0.8):
                s = env.reset()

            ep_reward = 0
            ep_ave_max_q = 0

            for j in range(int(args['max_episode_len'])):

                s_obs = env.unwrapped.state
                # env.render()   # Show animation

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

                # Utilize safety barrier function
                if (agent.firstIter == 1):
                    # 在第一个ep时的取值，此时未对GP_model_prev赋值.
                    [f, g, x, std] = dynamics_gp.get_GP_dynamics(agent, s, action_RL)
                else:
                    # 从第二个ep起，每个ep都将最新的GP_model赋值给GP_model_prev.
                    [f, g, x, std] = dynamics_gp.get_GP_dynamics_prev(agent, s, action_RL)
                # s:
                u_bar_ = cbf.control_barrier_qp(agent, np.squeeze(s), action_RL, f, g, x, std)
                # u_bar_ = cbf.control_barrier_sosp(agent, np.squeeze(s), action_RL, f, g, x, std, eng)

                action_ = action_RL + u_bar_  # (1, )

                s2, r, terminal, info = env.step(action_)  # r: reward
                s_next_unwrapped_state = env.unwrapped.state
                print('QP next state', s_next_unwrapped_state)
                print('QP control is:', u_bar_, '| Next state:', s_next_unwrapped_state[0], s_next_unwrapped_state[1])
                # TODO 注意维度
                # 此处存的时RL_policy原始的action, 仅对RL_Policy进行更新
                replay_buffer.add(np.reshape(s, (actor.s_dim,)), np.reshape(a, (actor.a_dim,)), r,
                                  terminal, np.reshape(s2, (actor.s_dim,)))

                # 这里计划存储的是RL_action经过u_bar与u_BAR调整的deployed action.
                # replay_buffer.add(np.reshape(s, (actor.s_dim,)), np.reshape(action_, (actor.a_dim,)), r,
                #                  terminal, np.reshape(s2, (actor.s_dim,)))

                # Keep adding experience to the memory until
                # there are at least minibatch size samples
                if replay_buffer.size() > int(args['minibatch_size']):
                    s_batch, a_batch, r_batch, t_batch, s2_batch = replay_buffer.sample_batch(
                        int(args['minibatch_size']))

                    # Calculate targets
                    target_q = critic.predict_target(s2_batch, actor.predict_target(s2_batch))

                    y_i = []
                    for k in range(int(args['minibatch_size'])):
                        if t_batch[k]:
                            y_i.append(np.array([r_batch[k]]))
                        else:
                            y_i.append(r_batch[k] + critic.gamma * target_q[k])

                    # Update the critic given the targets
                    # TODO 注意y_i 维度
                    # [1,2,..,64]

                    v_pred = critic.train(s_batch, a_batch, np.reshape(y_i, (int(args['minibatch_size']), 1)))

                    ep_ave_max_q += np.amax(v_pred)

                    # Update the actor policy using the sampled gradient
                    a_outs = actor.predict(s_batch)
                    s_grads = critic.action_gradients(s_batch, a_outs)
                    grads = np.array(s_grads).reshape((-1, actor.a_dim))

                    # grads维度
                    actor.train(s_batch, grads)

                    # Update target networks
                    actor.update_target_network()
                    critic.update_target_network()

                s = s2
                ep_reward += r

                obs.append(s)
                rewards.append(r)
                action_bar.append(u_bar_)
                action_BAR.append(u_BAR_)
                action.append(action_)
                action_RL_list.append(action_rl)
                unwrapped_state_lst.append(s_next_unwrapped_state) # theta, theta_dot

                if terminal:
                    reward_result[i * sub_episodes + el] = ep_reward
                    maxq_result[i * sub_episodes + el] = ep_ave_max_q / float(j)
                    obs_x = np.concatenate(unwrapped_state_lst).reshape((200, 2))
                    max_theta = np.max(obs_x[:, 0])
                    max_angle_result[i * sub_episodes + el] = max_theta
                    max_thetaDot = np.max(obs_x[:, 1])
                    max_angleSpeed_result[i * sub_episodes + el] = max_thetaDot

                    print('QP constrol is:', u_bar_, '| Next state:', s2[0], s2[1])
                    print('action_rl:', action_rl, '|', ' u_BAR:', u_BAR_, '|', 'u_bar:', u_bar_, '|', 'output:',
                          action_)
                    print(action[len(action) - 1], '|', action_bar[len(action_bar) - 1], '|',
                              action_RL_list[len(action_RL_list) - 1], '|', action_BAR[len(action_BAR) - 1])
                    # print(obs)
                    # print(obs.shape)
                    # breakpoint()
                    with writer.as_default():
                        tf.summary.scalar("Reward", ep_reward, step=i*sub_episodes+el)
                        tf.summary.scalar("Qmax Value",  ep_ave_max_q / float(j), step=i*sub_episodes+el)    
                        tf.summary.scalar("Max Angle", max_theta, step=i*sub_episodes+el)
                        tf.summary.scalar("Max Angle Speed", max_thetaDot, step=i*sub_episodes+el)
                    writer.flush()

                    print('| Episode: {:d} sub- {} | Reward: {:d} | Qmax: {:.4f} | MaxAngle: {} | MaxAngleSpeed: {} '.format(i, i*sub_episodes+el, int(ep_reward), (ep_ave_max_q / float(j)), max_theta, max_thetaDot))
                    print('-'*20)
                    # reward_result[i] = ep_reward
                    # maxq_result[i] = ep_ave_max_q / float(j)
                    # 200为最大max_time_steps
                    # (200, 3)
                    path = {"Observation": np.concatenate(obs).reshape(
                        (int(args['max_episode_len']), agent.observation_size)),  # (200, 3)
                            "Action": np.concatenate(action),
                            "Action_bar": np.concatenate(action_bar),
                            "Action_BAR": np.concatenate(action_BAR),
                            "Reward": np.asarray(rewards)}
                    paths.append(path)
                    cur_iter = i * sub_episodes + el + 1
                    save_dir = args['csv_dir']
                    os.makedirs(save_dir, exist_ok=True)
                    # max_angleSpeed_result
                    rec_pd = pd.DataFrame({'Reward': reward_result[:cur_iter], 
                                            'Qmax': maxq_result[:cur_iter], 
                                            'MaxAngle': max_angle_result[:cur_iter],
                                            'MaxAngleSpeed':max_angleSpeed_result[:cur_iter]},
                                            columns=['Reward', 'Qmax', 'MaxAngle', 'MaxAngleSpeed'])
                    rec_pd.to_csv(os.path.join(save_dir, 'log.csv'), index=False)     

                    # Log data in each step
                    step_data_dir = os.path.join(save_dir, 'step_data')
                    os.makedirs(step_data_dir, exist_ok=True)
                    step_data = pd.DataFrame({'Reward': rewards, 
                                              'Action':  np.concatenate(action), 
                                              'Angle': obs_x[:, 0],
                                              'AngleSpeed': obs_x[:, 1],
                                              'Action_RL':  np.concatenate(action_RL_list),
                                              'Action_ubar': np.concatenate(action_bar),
                                              'Action_uBAR':  np.concatenate(action_BAR)},
                                            columns=['Reward', 'Action', 'Angle', 'AngleSpeed', 'Action_RL', 'Action_ubar', 'Action_uBAR'])
                    step_data.to_csv(os.path.join(step_data_dir, '{}.csv'.format(cur_iter)), index=False)

                    # Trend graph
                    fig = plt.figure()
                    plt.plot(np.arange(len(maxq_result[:cur_iter])), maxq_result[:cur_iter])
                    plt.xlabel('Episode')
                    plt.ylabel('Max Q value')
                    pic_name = "Convergence of Max Q Value"
                    plt.title(pic_name)
                    img_path = os.path.join(save_dir, "MaxQ.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    fig = plt.figure()
                    plt.plot(np.arange(len(reward_result[:cur_iter])), reward_result[:cur_iter])
                    plt.xlabel('Episode')
                    plt.ylabel('Reward')
                    pic_name = "Reward vs. Episode"
                    plt.title(pic_name)
                    img_path = os.path.join(save_dir, "Reward.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    fig = plt.figure()
                    plt.plot(np.arange(len(max_angle_result[:cur_iter])), np.ones(len(max_angle_result[:cur_iter])), label='Safe Boundary', color='k', linestyle='--')
                    plt.plot(np.arange(len(max_angle_result[:cur_iter])), np.abs(max_angle_result[:cur_iter]), label='DDPG-{}'.format(args['method']), color='r')
                    plt.xlabel('Episode')
                    plt.ylabel('Max Angle (rad)')
                    pic_name = "Safety Violation"
                    plt.title(pic_name)
                    plt.legend()
                    img_path = os.path.join(save_dir, "MaxAngle.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    fig = plt.figure()
                    # max_angleSpeed_result
                    plt.plot(np.arange(len(max_angleSpeed_result[:cur_iter])), np.abs(max_angleSpeed_result[:cur_iter]), label='DDPG-{}'.format(args['method']), color='r')
                    plt.xlabel('Episode')
                    plt.ylabel('Max Angle Speed (rad/s)')
                    pic_name = "Trend of Max Angle Speed"
                    plt.title(pic_name)
                    plt.legend()
                    img_path = os.path.join(save_dir, "MaxAngleSpeed.png")
                    plt.savefig(img_path)
                    plt.close('all')

                    # First policy vs last policy

                    fig = plt.figure()
                    plt.plot(np.arange(int(args['max_episode_len'])), np.ones(int(args['max_episode_len'])), label='Safe Boundary', color='k', linestyle='--')
                    plt.plot(np.arange(int(args['max_episode_len'])), (-1) * np.ones(int(args['max_episode_len'])), color='k', linestyle='--')
                    if (i * sub_episodes + el) == 0:
                        step_policy_first = obs_x[:, 0] # 改
                    else:
                        step_policy_cur = obs_x[:, 0]
                        plt.plot(np.arange(int(args['max_episode_len'])), step_policy_cur,
                                 label='Epoch-{}'.format(i * sub_episodes + el + 1), color='b')
                    plt.plot(np.arange(int(args['max_episode_len'])), step_policy_first, label='1st Epoch', color='r')
                    plt.xlabel('Step')
                    plt.ylabel('Angle (rad)')
                    pic_name = "DDPG-{}".format(args['method'])
                    plt.title(pic_name)
                    plt.legend()
                    img_path = os.path.join(save_dir, "Step_PolicyCompare.png")
                    plt.savefig(img_path)
                    plt.close('all')
                    break

            if el <= 3:  # Update GP model each el except for the last el.
                dynamics_gp.update_GP_dynamics(agent, path)

        if (i <= 4):
            # Update the cbf_bar function.
            agent.bar_comp.get_training_rollouts(paths)
            barr_loss = agent.bar_comp.train()
        else:
            # 由于在训练后期，action通常都在safe area内，因此bar_loss都趋于0，为了加速训练，此处直接取0，且不影响结果。
            barr_loss = 0.
        agent.firstIter = 0  # 在第一个ep后变为0
    print("Training done")
    return [paths, reward_result]


def main(_argv):
    # gpu_devices = tf.config.experimental.list_physical_devices('GPU') # Search available GPU
    # tf.debugging.set_log_device_placement(True)
    # if len(gpu_devices) > 0:
    #     # os.environ['CUDA_VISIBLE_DEVICES'] = "0"
    #     tf.config.set_visible_devices(devices=gpu_devices, device_type='GPU') # Select GPU
    #     for single_gpu in gpu_devices:
    #         tf.config.experimental.set_memory_growth(single_gpu, True) # Set the graphics memory usage
    # else:
    #     print("Cannot find the GPU device.")

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    # eng = matlab.engine.connect_matlab()
    args = parse_args()

    matlab_eng_id = matlab.engine.find_matlab()
    eng = matlab.engine.connect_matlab(matlab_eng_id[0])

    # with tf.Session() as sess:

    env = gym.make(args['env'])  # unwrapped attr for raw settings in gym.
    np.random.seed(int(args['random_seed']))
    tf.random.set_seed(int(args['random_seed']))
    env.seed(int(args['random_seed']))

    # Set environment parameters for pendulum
    env.unwrapped.max_torque = 15.
    env.unwrapped.max_speed = 60.
    env.unwrapped.action_space = spaces.Box(low=-env.unwrapped.max_torque, high=env.unwrapped.max_torque, shape=(1,))
    high = np.array([1., 1., env.unwrapped.max_speed])
    env.unwrapped.observation_space = spaces.Box(low=-high, high=high)

    state_dim = env.observation_space.shape[0]  # 3
    action_dim = env.action_space.shape[0]  # 1
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

    agent = LEARNER(env)
    agent.max_episode_len = int(args['max_episode_len'])
    cbf.build_barrier(agent)  # 重复执行, 在LEARN.__init__已执行
    dynamics_gp.build_GP_model(agent)  # 重复执行, 在LEARN.__init__已执行
    agent.bar_comp = BARRIER(state_dim, action_dim, float(args['mlp_lr']),
                             args['mlp_mode'])  # 3: state_dim, 1: action_dim

    [paths, reward_result] = train(env, args, actor, critic, actor_noise, agent, eng)

    os.makedirs(args['mat_dir'], exist_ok=True)
    mat_path = os.path.join(args['mat_dir'], "data4_{}.mat".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M")))
    savemat(mat_path, dict(data=paths, reward=reward_result))
    return 0


if __name__ == '__main__':

    try:
        app.run(main)
    except SystemExit:
        pass
