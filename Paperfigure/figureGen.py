""" 
Figure generate
"""
import os

import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
plt.rc('font',family='Times New Roman') 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import datetime
from absl import app
from matplotlib.pyplot import MultipleLocator
from matplotlib.ticker import MaxNLocator



class FigGen():
    
    def __init__(self):
        self.y_scale_min = None # Specify the range of y_axis for Reward graph. None is for defult value.
        self.figFormat = 'eps' # 'eps'
        self.policy_ep_lst = [2, 150]
        self.step_range_display = 50
        cur_time = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        # self.save_dir = r'./res' # './res/{}'.format(cur_time)
        self.save_dir = r'C:\Users\hjhuang\Downloads\figGen\figGen\res' # './res/{}'.format(cur_time)
        
        linewidth = 2
        self.linewidth = linewidth
        self.dpi = 800
        self.config = {
            'ddpg': {
                'method': 'ONLY',
                'date': '2022-02-20-20-50-47', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                # 'res_dir': r'./data', # The path of result
                'res_dir': r'C:\Users\hjhuang\Downloads\figGen\figGen\data', # The path of result
                'color': '#F9A825', #'orange',
                'linestyle': '-.',
                'linewidth': linewidth,
            },

            'ddpg_qp': {
                'method': 'QP',
                'date': '2022-02-20-19-54-20', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                'res_dir': r'C:\Users\hjhuang\Downloads\figGen\figGen\data', # The path of result
                'color': '#388E3C', # 'green',
                'linestyle': '--',
                'linewidth': linewidth,
            },

            'ddpg_sosp': {
                'method': 'SOSP',
                'date': '2022-02-20-20-38-11', # Format: 'yyyy-mm-dd-hour-min-second', e.g., '2022-02-21 09:07:59'
                'res_dir': r'C:\Users\hjhuang\Downloads\figGen\figGen\data', # The path of result
                'color': '#2251D1', # 'blue',
                'linestyle': '-',
                'linewidth': linewidth,
            }
            # Reward, Qmax, MaxAngle, MaxAngleSpeed
        }
        self.algo_lst = list(self.config.keys())
        print("Algorithm list: {}".format(self.algo_lst))
        os.makedirs(self.save_dir, exist_ok=True)
        for algo in self.algo_lst:
            csv_dir = os.path.join(self.config[algo]['res_dir'], self.config[algo]['method'], self.config[algo]['date'])
            if not os.path.exists(csv_dir):
                raise Exception("Cannot find the dir [{}]".format(csv_dir))
            self.config[algo]['csv_dir'] = csv_dir
        
    def readData(self):
        
        self.algo_eps = []
        for algo in self.algo_lst:
            csv_dir = self.config[algo]['csv_dir']
            ep_data_pd = pd.read_csv(os.path.join(csv_dir, 'log.csv'), header=0)
            eps = 150 # len(ep_data_pd['Reward'])

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
            
            if self.step_range_display is None:
                self.config[algo]['step_size'] = len(self.config[algo]['last_epoch']['Angle'])
            else:
                self.config[algo]['step_size'] = self.step_range_display

            for mid_ep in self.policy_ep_lst:
                mid_ep_step_data = pd.read_csv(os.path.join(step_data_dir, '{}.csv'.format(mid_ep)), header=0)
                self.config[algo]['{}_epoch'.format(mid_ep)] = {}
                self.config[algo]['{}_epoch'.format(mid_ep)]['Reward'] = list(mid_ep_step_data['Reward'])
                self.config[algo]['{}_epoch'.format(mid_ep)]['Action'] = list(mid_ep_step_data['Action'])
                self.config[algo]['{}_epoch'.format(mid_ep)]['Angle'] = list(mid_ep_step_data['Angle'])
                self.config[algo]['{}_epoch'.format(mid_ep)]['AngleSpeed'] = list(mid_ep_step_data['AngleSpeed'])
                if algo != 'ddpg':
                    self.config[algo]['{}_epoch'.format(mid_ep)]['Action_RL'] = list(mid_ep_step_data['Action_RL'])
                    self.config[algo]['{}_epoch'.format(mid_ep)]['Action_ubar'] = list(mid_ep_step_data['Action_ubar'])
                    self.config[algo]['{}_epoch'.format(mid_ep)]['Action_uBAR'] = list(mid_ep_step_data['Action_uBAR'])



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
        if True:
            fig = plt.figure()
            for algo in self.algo_lst:
                # if algo == 'ddpg':
                #     continue
                plt.plot(np.arange(ep_min), self.config[algo]['Reward'][:ep_min], label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'], linewidth=self.config[algo]['linewidth'])
            font_size = 18
            plt.xlabel('Episode', fontsize=font_size)
            plt.ylabel('Reward', fontsize=font_size)
            plt.tick_params(labelsize=font_size-4) 
            plt.xlim(xmin=0, xmax=150)
            pic_name = "Reward vs. Episode"
            # plt.title(pic_name)
            plt.legend(bbox_to_anchor=(0.61, 0.78), prop={'size': font_size-5})
            x_major_locator=MultipleLocator(50)
            y_major_locator=MultipleLocator(1000)
            ax=plt.gca()
            ax.xaxis.set_major_locator(x_major_locator)
            ax.yaxis.set_major_locator(y_major_locator)
            
            # plt.ylim(ymin=-230, ymax=10)
            # if self.y_scale_min is not None:
                # plt.ylim(ymin=self.y_scale_min, ymax=0)

            # 子图
            inset_ax = fig.add_axes([0.58, 0.16, 0.3, 0.2],facecolor="white")
            for algo in self.algo_lst:
                inset_ax.plot(np.arange(ep_min), self.config[algo]['Reward'][:ep_min], label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'], linewidth=1) # self.config[algo]['linewidth'])
            # y_major_locator=MultipleLocator(1000)
            ax=plt.gca()
            # ax.yaxis.set_major_locator(y_major_locator)
            plt.tick_params(labelsize=font_size-8) # 8
            plt.ylim(ymin=-180, ymax=30)


            img_path = os.path.join(ep_dir, "Reward.{}".format(self.figFormat))
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            img_path = os.path.join(ep_dir, "Reward.png")
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            plt.close('all')

        """
        # Qmax
        fig = plt.figure()
        for algo in self.algo_lst:
            plt.plot(np.arange(ep_min), self.config[algo]['Qmax'][:ep_min], label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'])
        plt.xlabel('Episode')
        plt.ylabel('Max Q-value')
        pic_name = "Trend of Max Q-value"
        plt.title(pic_name)
        plt.legend()
        img_path = os.path.join(ep_dir, "MaxQ.{}".format(self.figFormat))
        plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
        img_path = os.path.join(ep_dir, "MaxQ.png")
        plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
        plt.close('all')
        """
        
        if True:
            # MaxAngle
            fig = plt.figure()
            font_size = 18 
            plt.plot(np.arange(ep_min), np.ones(ep_min), label='Safe Boundary', color='k', linestyle='-', linewidth=self.linewidth)
            for algo in self.algo_lst:
                plt.plot(np.arange(ep_min), np.abs(self.config[algo]['MaxAngle'][:ep_min]), label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'], linewidth=self.config[algo]['linewidth'])
            plt.xlabel('Episode', fontsize=font_size)
            plt.ylabel(r'Max |$\theta$|', fontsize=font_size)
            plt.tick_params(labelsize=font_size-4) 
            plt.yticks([0, 1, 3.14])
            plt.xlim(xmin=0, xmax=150)
            # plt.ylim(ymin=0)
            # plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True)) # y轴整数刻度
            x_major_locator=MultipleLocator(50)
            ax=plt.gca()
            ax.xaxis.set_major_locator(x_major_locator)
            

            # pic_name = "Safety Violation"
            # plt.title(pic_name)
            plt.legend(bbox_to_anchor=(0.97, 0.90), prop={'size': font_size-5})
            img_path = os.path.join(ep_dir, "MaxAngle.{}".format(self.figFormat))
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            img_path = os.path.join(ep_dir, "MaxAngle.png")
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            plt.close('all')
        
        if True:
        # MaxAngleSpeed
            fig = plt.figure()
            font_size = 18 
            # max_angleSpeed_result
            for algo in self.algo_lst:
                plt.plot(np.arange(ep_min), np.abs(self.config[algo]['MaxAngleSpeed'][:ep_min]), label='DDPG-{}'.format(self.config[algo]['method']), color=self.config[algo]['color'], linestyle=self.config[algo]['linestyle'], linewidth=self.config[algo]['linewidth'])
            plt.xlabel('Episode', fontsize=font_size)
            plt.ylabel(r'Max |$\dot{\theta}$|', fontsize=font_size)
            plt.tick_params(labelsize=font_size-4) 
            plt.xlim(xmin=0, xmax=150)

            x_major_locator=MultipleLocator(50)
            y_major_locator=MultipleLocator(2)
            ax=plt.gca()
            ax.xaxis.set_major_locator(x_major_locator)
            ax.yaxis.set_major_locator(y_major_locator)
            
            # pic_name = "Trend of Max Angle Speed"
            # plt.title(pic_name)
            plt.legend(bbox_to_anchor=(0.60, 0.75), prop={'size': font_size-5})
            img_path = os.path.join(ep_dir, "MaxAngleSpeed.{}".format(self.figFormat))
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            img_path = os.path.join(ep_dir, "MaxAngleSpeed.png")
            plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
            plt.close('all')

        print("Draw(ep) done..")

    def draw_step(self):
        step_dir = os.path.join(self.save_dir, 'step_fig')
        os.makedirs(step_dir, exist_ok=True)

        self.attr_dict = {
            'color': ['#2251D1', '#CE0202', '#F9A825', '#388E3C'],
            'linestyle': ['-', '--', '-.', ':']
        }
        font_size = 22
        # 各个算法分开画图
        for algo in self.algo_lst:
            # Reward, Action, Angle, AngleSpeed, Action_RL, Action_ubar, Action_uBAR, 
            self.data_lst = ['Reward', 'Action', 'Angle', 'AngleSpeed', 'Action_RL', 'Action_ubar', 'Action_uBAR']
            self.ylabel_lst = ['Reward', 'Torque', r'$\theta$', r'$\dot{\theta}$', 'Torque',  'Torque',  'Torque']
            for idx, data_type in enumerate(self.data_lst):
                if (algo == 'ddpg'): #and (data_type in ['Action_RL', 'Action_ubar', 'Action_uBAR']):
                    continue
                fig = plt.figure()
                if data_type == 'Angle':
                    plt.plot(np.arange(self.config[algo]['step_size'])+1, np.ones(self.config[algo]['step_size']), color='k', linestyle='-', linewidth=self.linewidth)
                    plt.plot(np.arange(self.config[algo]['step_size'])+1, (-1) * np.ones(self.config[algo]['step_size']), color='k', linestyle='-', linewidth=self.linewidth)
                else:
                    continue
                # plt.plot(np.arange(self.config[algo]['step_size']), self.config[algo]['first_epoch'][data_type][:self.config[algo]['step_size']], label='First policy', color=self.attr_dict['color'][0], linestyle=self.attr_dict['linestyle'][0], linewidth=self.linewidth)
                # plt.plot(np.arange(self.config[algo]['step_size']), self.config[algo]['last_epoch'][data_type][:self.config[algo]['step_size']], label='Last Epoch', color='red', linestyle='--')
                for mid_idx, mid_ep in enumerate(self.policy_ep_lst):
                    if mid_idx == 0:
                        cur_label = 'First policy'
                    elif mid_idx == (len(self.policy_ep_lst)-1):
                        cur_label = 'Final policy'
                    else:
                        cur_label = '{}-th policy'.format(mid_ep)
                    plt.plot(np.arange(self.config[algo]['step_size'])+1, self.config[algo]['{}_epoch'.format(mid_ep)][data_type][:self.config[algo]['step_size']], label=cur_label, color=self.attr_dict['color'][mid_idx], linestyle=self.attr_dict['linestyle'][mid_idx], linewidth=self.linewidth)
                
                plt.xlabel('Step', fontsize=font_size)
                plt.ylabel(self.ylabel_lst[idx], fontsize=font_size)
                plt.yticks([-1, 0, 1])
                plt.xticks([1, 25, 50])
                plt.tick_params(labelsize=font_size-2) 
                # plt.xlim(xmin=1)
                plt.ylim(ymin=-1.8, ymax=2.5) #-2.3, 2.3
                plt.xlim(xmin=1, xmax=50)

                pic_name = "DDPG-{}".format(self.config[algo]['method'])
                plt.title(pic_name, fontsize=font_size)
                
                if self.config[algo]['method'] == 'QP':
                    # 0.60, 0.70
                    plt.legend(bbox_to_anchor=(0.43, 0.66), prop={'size': font_size})
                elif self.config[algo]['method'] == 'SOSP':
                    # [0.97, 0.92]x起点重叠, 0.61, 0.80
                    plt.legend(bbox_to_anchor=(0.43, 0.66), prop={'size': font_size})
                    
                # plt.ylim(ymin=self.y_scale_min, ymax=0)
                img_path = os.path.join(step_dir, "DDPG-{}_{}.{}".format(self.config[algo]['method'], data_type, self.figFormat))
                plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
                img_path = os.path.join(step_dir, "DDPG-{}_{}.png".format(self.config[algo]['method'], data_type))
                plt.savefig(img_path, bbox_inches='tight', dpi=self.dpi)
                plt.close('all')
            print("Draw(step, {}) done..".format(algo))
        print("Draw(step) done..")


    def revise(self):
        for algo in self.algo_lst:
            csv_dir = self.config[algo]['csv_dir']
            # step_data
            step_data_dir = os.path.join(csv_dir, 'step_data')
            step_file_lst = os.listdir(step_data_dir)
            ep = len(step_file_lst)
            maxAngle_lst = []
            maxAngleSpeed_lst = []
            for sig_num in np.arange(1, ep + 1):
                csv_name = os.path.join(step_data_dir, '{}.csv'.format(sig_num))
                step_data = pd.read_csv(csv_name, header=0)
                step_data.rename(columns={'Angle': 'Angle_unwrap'}, inplace=True)
                step_data['Angle'] = np.arctan2(np.sin(step_data['Angle_unwrap']), np.cos(step_data['Angle_unwrap']))
                maxAngle_lst.append(np.max(np.abs(step_data['Angle'])))
                maxAngleSpeed_lst.append(np.max(np.abs(step_data['AngleSpeed'])))

                step_data.to_csv(csv_name, index=False)   
            log_name = os.path.join(csv_dir, 'log.csv')
            log_data = pd.read_csv(log_name, header=0)
            log_data.rename(columns={'MaxAngle': 'MaxAngle_raw', 'MaxAngleSpeed': 'MaxAngleSpeed_raw'}, inplace=True)
            log_data['MaxAngle'] = maxAngle_lst
            log_data['MaxAngleSpeed'] = maxAngleSpeed_lst
            log_data.to_csv(log_name, index=False)   


def main(_argv):
    """
    Epoch-level: Reawrd, max_angle, max_angle_speed
    Step-level (1st - last ep): Action, angle, angle_speed, reward ->compare
    """
    mode = "draw" # draw, revise
    fig_handler = FigGen() 
    if mode == 'draw':
        fig_handler.readData()
        fig_handler.draw_ep()
        fig_handler.draw_step()
    elif mode == 'revise':
        fig_handler.revise()
    else:
        raise Exception("Unexpected mode [{}]".format(mode))

    return 0

if __name__ == '__main__':

    try:
        app.run(main)
    except SystemExit:
        pass


"""
1.修改step data中Angle的wrap数值，先将原始数据重命名为angle_unwrap, 再生成新的angle
2.修改Max angle, max angle speed的计算方式max(abs())，重新生成log.csv
3.

"""