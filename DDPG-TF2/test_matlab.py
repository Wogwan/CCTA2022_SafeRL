import matlab
import matlab.engine
import os
import gym
import sys
import numpy as np

print("Current working directory: {0}".format(os.getcwd()))


def build_barrier(eng):
    eng.cd(r'C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\matlab_code', nargout=0)
    aa = 1
    c, d, e, f = eng.main_function(aa, nargout=4)
    return c, d, e, f

def validate_safe_state(x1, x2, sys_bar, region, sys_d2, eng):
    eng.cd(r'C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\matlab_code', nargout=0)
    c = eng.check_state(x1, x2, sys_bar, region, sys_d2, nargout=1)
    return c


if __name__ == '__main__':
    matlab_eng_id = matlab.engine.find_matlab()
    eng = matlab.engine.connect_matlab(matlab_eng_id[0])
    env = gym.make('Pendulum-v1')
    s = env.reset()
    record = env.unwrapped.state
    coe, sys_bar, region, sys_d2 = build_barrier(eng)
    kk = record.tolist()
    s1 = kk[0]
    s2 = kk[1]
    k = validate_safe_state(s1, s2, sys_bar, region, sys_d2, eng)
    print(k)

# MATLAB Usage
# PATH = engine.sosaddpath
# print(PATH)
# eng = matlab.engine.start_matlab()
# print("Current working directory: {0}".format(os.getcwd()))
# k = 1
# a = engine.mat_py_test(k,nargin=1)
# a = engine.mat_py_test
# a = engine.workspace

# eng.cd(r'/home/wogwan/Documents/safeRL/ASCC2022_SafeRL-dev/DDPG-TF2/matlab_code',nargout=0)
# Important to specify the output
# c = eng.eval('main_function',nargout = 2)
# print(c)
# eng.quit()

# engine.mat_py_test(nargout=0)
# engine.quit()

