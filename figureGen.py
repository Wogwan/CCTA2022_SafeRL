""" 
Figure generate
"""
import os

import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
import datetime
from absl import app

class FigGen():
    
    def __init__(self):
        self.y_scale_min = None # Specify the range of y_axis for Reward graph. None is for defult value.
        cur_time = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        self.save_dir = './res/figGen/{}'.format(cur_time)
        self.config = {
            'ddpg': {
                'method': 'ONLY',
                'date': '2022-02-20-20-50-47', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                'res_dir': './DDPG-TF2/res', # The path of result
                'color': 'yellow',
                'linestyle': '--',
            },

            'ddpg_qp': {
                'method': 'QP',
                'date': '2022-02-20-19-54-20', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                'res_dir': './DDPG-CBF-TF2/res', # The path of result
                'color': 'blue',
                'linestyle': '-.',
            },

            'ddpg_sosp': {
                'method': 'SOSP',
                'date': '2022-02-20-20-38-11', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                'res_dir': './DDPG-CBF-TF2/res', # The path of result
                'color': 'red',
                'linestyle': '-',
            }
            # Reward, Qmax, MaxAngle, MaxAngleSpeed
        }
        self.algo_lst = list(self.config.keys())
        print("Algorithm list: {}".format(self.algo_lst))
        os.makedirs(self.save_dir, exist_ok=True)
        for algo in self.algo_lst:
            csv_dir = os.path.join(self.config[algo]['res_dir'], self.config[algo]['method'], 'csv', self.config[algo]['date'])
            if not os.path.exists(csv_dir):
                raise Exception("Cannot find the dir [{}]".format(csv_dir))
            self.config[algo]['csv_dir'] = csv_dir
        
    def readData(self):
        
        self.algo_eps = []
        for algo in self.algo_lst:
            csv_dir = self.config[algo]['csv_dir']
            ep_data_pd = pd.read_csv(os.path.join(csv_dir, 'log.csv'), header=0)
            eps = 137 # len(ep_data_pd['Reward'])

            self.config[algo]['Reward'] = list(ep_data_pd['Reward'])[:eps]
            self.config[algo]['Qmax'] = list(ep_data_pd['Qmax'])[:eps]
            self.config[algo]['MaxAngle'] = list(ep_data_pd['MaxAngle'])[:eps]
            self.config[algo]['MaxAngleSpeed'] = list(ep_data_pd['MaxAngleSpeed'])[:eps]
            self.config[algo]['epoch'] = eps
            self.algo_eps.append(eps)

            step_data_dir = os.path.join(csv_dir, 'step_data')
            first_ep_step_data = pd.read_csv(os.path.join(step_data_dir, '1.csv'), header=0)
            self.config[algo]['first_epoch'] = {}
            self.config[algo]['first_epoch']['Reward'] = list(first_ep_step_data['Reward'])
            self.config[algo]['first_epoch']['Action'] = list(first_ep_step_data['Action'])
            self.config[algo]['first_epoch']['Angle'] = list(first_ep_step_data['Angle'])
            self.config[algo]['first_epoch']['AngleSpeed'] = list(first_ep_step_data['AngleSpeed'])
            if algo != 'ddpg':
                self.config[algo]['first_epoch']['Action_RL'] = list(first_ep_step_data['Action_RL'])
                self.config[algo]['first_epoch']['Action_ubar'] = list(first_ep_step_data['Action_ubar'])
                self.config[algo]['first_epoch']['Action_uBAR'] = list(first_ep_step_data['Action_uBAR'])

            last_ep_step_data = pd.read_csv(os.path.join(step_data_dir, '{}.csv'.format(eps)), header=0)
            self.config[algo]['last_epoch'] = {}
            self.config[algo]['last_epoch']['Reward'] = list(last_ep_step_data['Reward'])
            self.config[algo]['last_epoch']['Action'] = list(last_ep_step_data['Action'])
            self.config[algo]['last_epoch']['Angle'] = list(last_ep_step_data['Angle'])
            self.config[algo]['last_epoch']['AngleSpeed'] = list(last_ep_step_data['AngleSpeed'])
            if algo != 'ddpg':
                self.config[algo]['last_epoch']['Action_RL'] = list(last_ep_step_data['Action_RL'])
                self.config[algo]['last_epoch']['Action_ubar'] = list(last_ep_step_data['Action_ubar'])
                self.config[algo]['last_epoch']['Action_uBAR'] = list(last_ep_step_data['Action_uBAR'])
            self.config[algo]['step_size'] = len(self.config[algo]['last_epoch']['Angle'])


        # 提示三个算法长度不一致
        flag = np.sum(np.diff(self.algo_eps))
        if flag != 0:
            print("NOTICE: The epoches are not equal.. [{}]".format(self.algo_eps))
        print("Read data done..")

    def draw_ep(self):
        ep_dir = os.path.join(self.save_dir, 'ep_fig')
        os.makedirs(ep_dir, exist_ok=True)
        ep_min = np.min(self.algo_eps) # Min. ep among algorithms.
        # Reward
        fig = plt.figure()
        for algo in self.algo_lst:
            plt.plot(np.arange(ep_min), self.config[algo]['Reward'][:ep_min], label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'])
        plt.xlabel('Episode')
        plt.ylabel('Reward')
        pic_name = "Reward vs. Episode"
        plt.title(pic_name)
        plt.legend()
        if self.y_scale_min is not None:
            plt.ylim(ymin=self.y_scale_min, ymax=0)
        img_path = os.path.join(ep_dir, "Reward.png")
        plt.savefig(img_path)
        plt.close('all')

        # Qmax
        fig = plt.figure()
        for algo in self.algo_lst:
            plt.plot(np.arange(ep_min), self.config[algo]['Qmax'][:ep_min], label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'])
        plt.xlabel('Episode')
        plt.ylabel('Max Q-value')
        pic_name = "Trend of Max Q-value"
        plt.title(pic_name)
        plt.legend()
        img_path = os.path.join(ep_dir, "MaxQ.png")
        plt.savefig(img_path)
        plt.close('all')

        # MaxAngle
        fig = plt.figure()
        plt.plot(np.arange(ep_min), np.ones(ep_min), label='Safe Boundary', color='k', linestyle='--')
        for algo in self.algo_lst:
            plt.plot(np.arange(ep_min), np.abs(self.config[algo]['MaxAngle'][:ep_min]), label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'])
        plt.xlabel('Episode')
        plt.ylabel('Max Angle (rad)')
        pic_name = "Safety Violation"
        plt.title(pic_name)
        plt.legend()
        img_path = os.path.join(ep_dir, "MaxAngle.png")
        plt.savefig(img_path)
        plt.close('all')

        # MaxAngleSpeed
        fig = plt.figure()
        # max_angleSpeed_result
        for algo in self.algo_lst:
            plt.plot(np.arange(ep_min), np.abs(self.config[algo]['MaxAngleSpeed'][:ep_min]), label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'])
        plt.xlabel('Episode')
        plt.ylabel('Max Angle Speed (rad/s)')
        pic_name = "Trend of Max Angle Speed"
        plt.title(pic_name)
        plt.legend()
        img_path = os.path.join(ep_dir, "MaxAngleSpeed.png")
        plt.savefig(img_path)
        plt.close('all')

        print("Draw(ep) done..")

    def draw_step(self):
        step_dir = os.path.join(self.save_dir, 'step_fig')
        os.makedirs(step_dir, exist_ok=True)
        # 各个算法分开画图
        for algo in self.algo_lst:
            # Reward, Action, Angle, AngleSpeed, Action_RL, Action_ubar, Action_uBAR, 
            self.data_lst = ['Reward', 'Action', 'Angle', 'AngleSpeed', 'Action_RL', 'Action_ubar', 'Action_uBAR']
            self.ylabel_lst = ['Reward', 'Torque', 'Angle (rad)', 'Angle Speed (rad/s)', 'Torque',  'Torque',  'Torque']
            for idx, data_type in enumerate(self.data_lst):
                if (algo == 'ddpg') and (data_type in ['Action_RL', 'Action_ubar', 'Action_uBAR']):
                    continue
                fig = plt.figure()
                if data_type == 'Angle':
                    plt.plot(np.arange(self.config[algo]['step_size']), np.ones(self.config[algo]['step_size']), label='Safe Boundary', color='k', linestyle='--')
                    plt.plot(np.arange(self.config[algo]['step_size']), (-1) * np.ones(self.config[algo]['step_size']), color='k', linestyle='--')

                plt.plot(np.arange(self.config[algo]['step_size']), self.config[algo]['first_epoch'][data_type], label='First Epoch', color='blue', linestyle='-')
                plt.plot(np.arange(self.config[algo]['step_size']), self.config[algo]['last_epoch'][data_type], label='Last Epoch', color='red', linestyle='--')
                plt.xlabel('Step')
                plt.ylabel(self.ylabel_lst[idx])
                pic_name = "DDPG-{}({})".format(self.config[algo]['method'], data_type)
                plt.title(pic_name)
                plt.legend()
                img_path = os.path.join(step_dir, "DDPG-{}_{}.png".format(self.config[algo]['method'], data_type))
                plt.savefig(img_path)
                plt.close('all')
            print("Draw(step, {}) done..".format(algo))
        print("Draw(step) done..")


def main(_argv):
    """
    Epoch-level: Reawrd, max_angle, max_angle_speed
    Step-level (1st - last ep): Action, angle, angle_speed, reward ->compare
    """
    fig_handler = FigGen() 
    fig_handler.readData()

    fig_handler.draw_ep()
    fig_handler.draw_step()


    return 0

if __name__ == '__main__':

    try:
        app.run(main)
    except SystemExit:
        pass
