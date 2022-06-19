import numpy as np
import cbf
import dynamics_gp

class LEARNER():
    def __init__(self, env):
        self.firstIter = 1
        self.count = 1
        self.env = env
        self.torque_bound = 15.
        self.max_speed = 60.

        self.observation_space = env.observation_space
        self.action_space = env.action_space
        print('Observation space: {}'.format(self.observation_space))
        print('Action space: {}'.format(self.action_space))

        self.observation_size = self.env.observation_space.shape[0]
        self.action_size = self.action_space.shape[0]
        
        # Build barrier function model
        cbf.build_barrier(self)
        
        # Build GP model of dynamics
        dynamics_gp.build_GP_model(self)


