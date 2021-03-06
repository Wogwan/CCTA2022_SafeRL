#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:17:49 2018

@author: rcheng
"""

import numpy as np
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
from cvxopt import matrix
from cvxopt import solvers
import dynamics_gp
import matlab

#Build barrier function model
def build_barrier(self):
    N = self.action_size
    #self.P = matrix(np.eye(N), tc='d')
    self.P = matrix(np.diag([1., 1e24]), tc='d')
    self.q = matrix(np.zeros(N+1))
    self.H1 = np.array([1, 0.05]) # y1 = H1[0]s[0] + H1[1]s[1]
    self.H2 = np.array([1, -0.05])
    self.H3 = np.array([-1, 0.05])
    self.H4 = np.array([-1, -0.05])
    self.F = 1

#Get compensatory action based on satisfaction of barrier function
def control_barrier_sosp(self, obs, u_rl, f, g, x, std, eng):
    """
    Input: 
    obs: observation: (3, ), current state
    u_rl: action: (1, ), a_RL + u_BAR
    f: (2, ), [next_theta, next_theta_dot]
    g: (1, ), g(s_t)
    x: (2, ), [original_theta, original_theta_dot]
    std: (2, ), [GP_Model_1_std, GP_Model_2_std]

    Output:
    u_bar: (1, )
    """

    '''
        Solves a SOS program

            minimize    c_1*u_1 + c_2*u_2 + ... + c_n*u_n = u
            subject to  Δh - L_1*h >= 0.
                        
    '''
    # eng.cd(r'C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\matlab_code', nargout=0)
    x2 = x
    x1 = x[0].tolist()
    x2 = x[1].tolist()
    u_rl = u_rl[0].tolist()
    u_sos = eng.sos_compute_nosim(u_rl,x1,x2,nargout=1)
    return np.array([u_sos])


def control_barrier_qp(self, obs, u_rl, f, g, x, std): # Original name control_barrier()
    """
    Input: 
    obs: observation: (3, ), current state
    u_rl: action: (1, ), a_RL + u_BAR
    f: (2, ), [next_theta, next_theta_dot]
    g: (1, ), g(s_t)
    x: (2, ), [original_theta, original_theta_dot]
    std: (2, ), [GP_Model_1_std, GP_Model_2_std]

    Output:
    u_bar: (1, )
    """
    #Define gamma for the barrier function
    gamma_b = 0.5
    
    #Set up Quadratic Program to satisfy the Control Barrier Function
    kd = 1.5
    u_a = 0

    G = np.array([[-np.dot(self.H1,g), 
                    -np.dot(self.H2,g), 
                    -np.dot(self.H3,g), 
                    -np.dot(self.H4,g), 
                    1, 
                    -1, 
                    g[1], 
                    -g[1]], 
                    [-1, -1, -1, -1, 0, 0, 0, 0]])
    G = np.transpose(G)

    '''
        Solves a quadratic program

            minimize    (1/2)*x'*P*x + [q'*x]
            subject to  G*x <= h.
                        [A*x = b]
    '''

    h = np.array([gamma_b*self.F + np.dot(self.H1,f) + np.dot(self.H1,g)*u_a - (1-gamma_b)*np.dot(self.H1,x) - kd*np.dot(np.abs(self.H1),std),
                # gamma_b*self.F + np.dot(self.H1,f) + np.dot(self.H1,g)*u_a - (1-gamma_b)*np.dot(self.H1,x) - kd*np.dot(np.abs(self.H1),std)
                # f = f(s_{t+1}) = f(s_t) + g(s_t)a_t + \mu(s_t) - k\sigma(s_t)
                # - (1-gamma_b)
                  gamma_b*self.F + np.dot(self.H2,f) + np.dot(self.H2,g)*u_a - (1-gamma_b)*np.dot(self.H2,x) - kd*np.dot(np.abs(self.H2),std),
                  gamma_b*self.F + np.dot(self.H3,f) + np.dot(self.H3,g)*u_a - (1-gamma_b)*np.dot(self.H3,x) - kd*np.dot(np.abs(self.H3),std),
                  gamma_b*self.F + np.dot(self.H4,f) + np.dot(self.H4,g)*u_a - (1-gamma_b)*np.dot(self.H4,x) - kd*np.dot(np.abs(self.H4),std),
                  -u_rl + self.torque_bound, # -u + b > 0, u<b
                  u_rl + self.torque_bound, # u + b >0 , u > -b
                  -f[1] - g[1]*u_rl + self.max_speed,
                  f[1] + g[1]*u_rl + self.max_speed])
    h = np.squeeze(h).astype(np.double) #(1,1,3,5,1) => (3,5), (8, ) => (8, )
    
    #Convert numpy arrays to cvx matrices to set up QP
    G = matrix(G,tc='d')
    h = matrix(h,tc='d')

    solvers.options['show_progress'] = False
    sol = solvers.qp(self.P, self.q, G, h)
    u_bar = sol['x']
    #if np.abs(u_bar[1]) > 0.001:
        #print("Violation of Safety: ")
        #print(u_bar[1])

    if (np.add(np.squeeze(u_rl), np.squeeze(u_bar[0])) - 0.001 >= self.torque_bound):
        u_bar[0] = self.torque_bound - u_rl
        print("Error in QP")
    elif (np.add(np.squeeze(u_rl), np.squeeze(u_bar[0])) + 0.001 <= -self.torque_bound):
        u_bar[0] = -self.torque_bound - u_rl
        print("Error in QP")
    else:
        pass

    return np.expand_dims(np.array(u_bar[0]), 0)
