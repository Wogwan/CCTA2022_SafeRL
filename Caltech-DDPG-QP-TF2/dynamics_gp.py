#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:17:49 2018

@author: rcheng
"""

import numpy as np
from cvxopt import matrix
from cvxopt import solvers
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C

def build_GP_model(self):
    N = self.observation_size # 3
    GP_list = []
    noise = 0.01
    for i in range(N-1): # 0, 1
        kern = 1.0 * RBF(length_scale=2.0, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-10, 1e+1))
        #kern = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        gp = GaussianProcessRegressor(kernel=kern, alpha = noise, n_restarts_optimizer=10)
        GP_list.append(gp)
    self.GP_model = GP_list

#Get the dynamics of the system from the current time step with the RL action
def get_dynamics(self, obs, u_rl):
    dt = 0.05
    G = 10
    m = 1.4
    l = 1.4
    obs = np.squeeze(obs)
    theta = np.arctan2(obs[1], obs[0])
    theta_dot = obs[2]
    # 为什么(theta + np.pi)要+pi并在前面*(-1)
    # f[next_theta, next_theta_dot], 
    f = np.array([-3*G/(2*l)*np.sin(theta + np.pi)*dt**2 + theta_dot*dt + theta + 3/(m*l**2)*u_rl*dt**2, 
                    theta_dot - 3*G/(2*l)*np.sin(theta + np.pi)*dt + 3/(m*l**2)*u_rl*dt])
    # g(st), 即h(s)中的g(st)，与action相乘的项
    g = np.array([3/(m*l**2)*dt**2, 3/(m*l**2)*dt])
    
    x = np.array([theta, theta_dot])
    return [np.squeeze(f), np.squeeze(g), np.squeeze(x)]

#Build barrier function model
def update_GP_dynamics(self,path):
    N = self.observation_size
    X = path['Observation']
    X = X.reshape((self.max_episode_len, self.observation_size)) # (200, 3)) 200是max_time_steps, 最多200步就会terminal.
    U = path['Action'] # action_rl + u_BAR + u_bar
    L = X.shape[0]
    err = np.zeros((L-1,N-1))
    S = np.zeros((L-1,2))
    for i in range(L-1):
        [f,g,x1] = get_dynamics(self,X[i,:],U[i])
        theta_p = np.arctan2(X[i,1], X[i,0])
        theta_dot_p = X[i,2]
        theta = np.arctan2(X[i+1,1], X[i+1,0])
        theta_dot = X[i+1,2]
        S[i,:] = np.array([theta_p, theta_dot_p])
        # f为根据公式以及当前状态与动作，得到下一步计算的确定的状态变化。
        # 因此这里模拟实际环境下一步状态与计算的下一步状态的差距diff，
        # 并用这个diff作为label去拟合GP model.
        err[i,:] = np.array([theta, theta_dot]) - f 
    self.GP_model[0].fit(S,err[:,0])
    self.GP_model[1].fit(S,err[:,1])
    
def get_GP_dynamics(self, obs, u_rl):
    [f_nom, g, x] = get_dynamics(self, obs, u_rl)
    f = np.zeros(2)
    # 根据上一状态state->x(即last_theta, last_theta_dot),预测本次GP uncertainty 
    [m1, std1] = self.GP_model[0].predict(x.reshape(1,-1), return_std=True)
    [m2, std2] = self.GP_model[1].predict(x.reshape(1,-1), return_std=True)
    """
    [
       [[1,2],[3,4] ], example 1
       [[5,6],[7,8] ], example 2
    ]
    
    """
    f[0] = f_nom[0] + self.GP_model[0].predict(x.reshape(1,-1)) 
    f[1] = f_nom[1] + self.GP_model[1].predict(x.reshape(1,-1))
    return [np.squeeze(f), np.squeeze(g), np.squeeze(x), np.array([np.squeeze(std1), np.squeeze(std2)])]

def get_GP_dynamics_prev(self, obs, u_rl):
    [f_nom, g, x] = get_dynamics(self, obs, u_rl)
    f = np.zeros(2)
    [m1, std1] = self.GP_model_prev[0].predict(x.reshape(1,-1), return_std=True)
    [m2, std2] = self.GP_model_prev[1].predict(x.reshape(1,-1), return_std=True)
    f[0] = f_nom[0] + self.GP_model_prev[0].predict(x.reshape(1,-1))
    f[1] = f_nom[1] + self.GP_model_prev[1].predict(x.reshape(1,-1))
    return [np.squeeze(f), np.squeeze(g), np.squeeze(x), np.array([np.squeeze(std1), np.squeeze(std2)])]

