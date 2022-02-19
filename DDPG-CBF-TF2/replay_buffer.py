""" 
Data structure for implementing experience replay

Author: Patrick Emami
"""
from collections import deque
import random
import numpy as np

class ReplayBuffer(object):

    def __init__(self, buffer_size, random_seed=123):
        """
        The right side of the deque contains the most recent experiences 
        """
        self.buffer_size = buffer_size
        self.buffer = deque(maxlen=self.buffer_size)
        random.seed(random_seed)

    def add(self, s, a, r, t, s2):
        # s: state, a: action, r: reward, t: terminal, s2: next state
        self.buffer.append([s, a, r, t, s2])
        
        # experience = (s, a, r, t, s2)
        # if self.count < self.buffer_size: 
        #     self.buffer.append(experience)
        #     self.count += 1
        # else:
        #     self.buffer.popleft()
        #     self.buffer.append(experience)

    def size(self):
        return len(self.buffer)

    def sample_batch(self, batch_size):
        batch = []

        if len(self.buffer) < batch_size:
            batch = random.sample(self.buffer, len(self.buffer))
        else:
            batch = random.sample(self.buffer, batch_size)
        
        s_batch, a_batch, r_batch, t_batch, s2_batch = map(np.asarray, zip(*batch))
        """
        [ [s1, a1, r1, t1, ns1], [s2, a2, r2, .., ns2], ..]
        After line 43:
        [ s_batch: [s1, s2, s3, s4, .., s64], a_batch: [a1, a2, a3, ..., a64] , [r1, .., r64] ]
        """
        # s_batch = np.array([_[0] for _ in batch])
        # a_batch = np.array([_[1] for _ in batch])
        # r_batch = np.array([_[2] for _ in batch])
        # t_batch = np.array([_[3] for _ in batch])
        # s2_batch = np.array([_[4] for _ in batch])

        return s_batch, a_batch, r_batch, t_batch, s2_batch

    def clear(self):
        self.buffer.clear()
