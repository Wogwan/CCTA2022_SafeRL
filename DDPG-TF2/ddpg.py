""" 
Implementation of DDPG - Deep Deterministic Policy Gradient

Algorithm and hyperparameter details can be found here: 
    http://arxiv.org/pdf/1509.02971v2.pdf

The algorithm is tested on the Pendulum-v0 OpenAI gym task 
and developed with tflearn + Tensorflow

Author: Patrick Emami
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "3"
import tensorflow as tf

import numpy as np
import pandas as pd
import gym
from gym import wrappers
# import tflearn
import argparse
import pprint as pp
from scipy.io import savemat
import datetime
from absl import app
from replay_buffer import ReplayBuffer
from matplotlib import pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='provide arguments for DDPG agent')

    # agent parameters
    parser.add_argument('--actor-lr', help='actor network learning rate', default=0.0001)
    parser.add_argument('--critic-lr', help='critic network learning rate', default=0.001)
    parser.add_argument('--gamma', help='discount factor for critic updates', default=0.99)
    parser.add_argument('--tau', help='soft target update parameter', default=0.001)
    parser.add_argument('--buffer-size', help='max size of the replay buffer', default=1000000)
    parser.add_argument('--minibatch-size', help='size of minibatch for minibatch-SGD', default=64)

    # run parameters
    parser.add_argument('--env', help='choose the gym env- tested on {Pendulum-v0}', default='Pendulum-v1') # Pendulum-v0
    parser.add_argument('--random-seed', help='random seed for repeatability', default=1754)
    parser.add_argument('--max-episodes', help='max num of episodes to do while training', default=600) # 600
    parser.add_argument('--max-episode-len', help='max length of 1 episode', default=200)
    parser.add_argument('--render-env', help='render the gym env', action='store_false')
    parser.add_argument('--use-gym-monitor', help='record gym results', action='store_false')
    parser.add_argument('--monitor-dir', help='directory for storing gym results', default='./res/gym_ddpg')
    parser.add_argument('--summary-dir', help='directory for storing tensorboard info', default='./res/tf_ddpg')
    parser.add_argument('--mat-dir', help='directory for storing mat info', default='./res/mat')
    parser.add_argument('--csv-dir', help='directory for storing csv info', default='./res/csv')


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
        self.action_bound = action_bound # tf.constant(action_bound)
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
        self.num_trainable_vars = len(self.actor_model.trainable_variables) + len(self.target_actor_model.trainable_variables)

    def create_actor_network(self):
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
        # Scale output to -action_bound to action_bound
        scaled_out = tf.keras.layers.Multiply()([out, self.action_bound])
        actor_model = tf.keras.models.Model(inputs=inputs,outputs=scaled_out)
        return actor_model

    def train(self, inputs, a_gradient):
        # inputs: state
        # a_gradient: the graident of critic model loss
        with tf.GradientTape() as tape:
            unnormal_grads = tape.gradient(self.actor_model(inputs, training=True), self.actor_model.trainable_variables, -a_gradient)
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
        # inputs: states
        # action
        return self.critic_model.predict([inputs, action])

    def predict_target(self, inputs, action):
        return self.target_critic_model.predict([inputs, action])

    def action_gradients(self, inputs, actions):
        actions = tf.convert_to_tensor(actions)
        with tf.GradientTape() as tape:
            tape.watch(actions)
            q_values_x = self.critic_model([inputs, actions])
            q_values = tf.squeeze(q_values_x) # The same gradient can be obtained if remove it.

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

def train(env, args, actor, critic, actor_noise):

    # Set up summary Ops
    # summary_ops, summary_vars = build_summaries()

    writer = tf.summary.create_file_writer(args['summary_dir'])
    reward_result = np.zeros(int(args['max_episodes'])) # Raw code: 2500
    maxq_result = np.zeros(int(args['max_episodes'])) # Raw code: 2500

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
        print("Epoch: {}".format(i))
        s = env.reset()

        ep_reward = 0
        ep_ave_max_q = 0

        obs, action, rewards, action_bar, action_BAR, action_RL_list = [], [], [], [], [], []

        for j in range(int(args['max_episode_len'])):

            # env.render()


            # Added exploration noise
            #a = actor.predict(np.reshape(s, (1, 3))) + (1. / (1. + i))
            a = actor.predict(np.reshape(s, (1, actor.s_dim))) + actor_noise()

            s2, r, terminal, info = env.step(a[0])
            s_next_state = env.unwrapped.state
            print('DDPG control is:', a[0], '| Next state:', s_next_state[0], s_next_state[1])
            replay_buffer.add(np.reshape(s, (actor.s_dim,)), np.reshape(a, (actor.a_dim,)), r,
                              terminal, np.reshape(s2, (actor.s_dim,)))

            # Keep adding experience to the memory until
            # there are at least minibatch size samples
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

                # Update the critic given the targets
                # predicted_q_value是(s_t, a_t)通过critic预测得到的predicted_q_value
                predicted_q_value = critic.train(
                    s_batch, a_batch, np.reshape(y_i, (int(args['minibatch_size']), 1)))

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

            obs.append(s)
            rewards.append(r)
            action.append(a[0])

            if terminal:
                print('DDPG control is:', a[0], '| Next state:', s2[0], s2[1])
                print('output:', a[0])

                with writer.as_default():
                    # print(action[len(action) - 1], '|', action_bar[len(action_bar) - 1], '|',
                    #       action_RL_list[len(action_RL_list) - 1], '|', action_BAR[len(action_BAR) - 1])

                    tf.summary.scalar("Reward", ep_reward, step=i)
                    tf.summary.scalar("Qmax Value",  ep_ave_max_q / float(j), step=i)
                writer.flush()

                print('| Reward: {:d} | Episode: {:d} | Qmax: {:.4f}'.format(int(ep_reward), \
                        i, (ep_ave_max_q / float(j))))
                reward_result[i] = ep_reward
                maxq_result[i] = ep_ave_max_q / float(j)
                path = {"Observation":np.concatenate(obs).reshape((200,3)),
                        "Action":np.concatenate(action),
		                "Reward":np.asarray(rewards)}
                paths.append(path)

                os.makedirs(args['csv_dir'], exist_ok=True)
                rec_pd = pd.DataFrame({'Reward': reward_result[:i+1], 'Qmax': maxq_result[:i+1]}, columns=['Reward', 'Qmax'])
                rec_pd.to_csv(os.path.join(args['csv_dir'], 'log.csv'), index=False)

                # Trend graph
                fig = plt.figure(figsize=(12, 6))
                ax1 = plt.subplot(1, 2, 1)
                ax1.plot(np.arange(len(maxq_result[:i+1])), maxq_result[:i+1])
                ax1.set_xlabel('Episode')
                ax1.set_ylabel('Max Q value')
                pic_name = "Convergence of Max Q Value"
                plt.title(pic_name)

                ax2 = plt.subplot(1, 2, 2)
                ax2.plot(np.arange(len(reward_result[:i+1])), reward_result[:i+1])
                ax2.set_xlabel('Episode')
                ax2.set_ylabel('Reward')
                pic_name = "Reward vs. Episode"
                plt.title(pic_name)
                img_path = os.path.join(args['csv_dir'], "Trend.png")
                plt.savefig(img_path)
                plt.close('all')

                break

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

    args = parse_args()


    env = gym.make(args['env'])
    np.random.seed(int(args['random_seed']))
    tf.random.set_seed(int(args['random_seed']))
    env.seed(int(args['random_seed']))

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

    [paths, reward_result] = train(env, args, actor, critic, actor_noise)

    os.makedirs(args['mat_dir'], exist_ok=True)
    mat_path = os.path.join(args['mat_dir'], "data4_{}.mat".format(datetime.datetime.now().strftime("%y-%m-%d-%H-%M")))
    savemat(mat_path, dict(data=paths, reward=reward_result))

    # Save to csv
    return 0

if __name__ == '__main__':

    try:
        app.run(main)
    except SystemExit:
        pass