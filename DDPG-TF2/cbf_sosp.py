import numpy
import os

import numpy as np
import scipy.io
import matlab
import matlab.engine
# import os
# import sys
# sys.path.append('/home/wogwan/Documents/safeRL/ASCC2022_SafeRL-dev/DDPG-TF2/matlab_code')
# os.chdir(r'/home/wogwan/Documents/safeRL/ASCC2022_SafeRL-dev/DDPG-TF2/matlab_code')


def build_barrier(self, eng):
    # Start running matlab
    # eng = matlab.engine.start_matlab()
    # Save some observations
    # scipy.io.savemat(self.obs_file_path,{'data':self.obs})
    # Connect to an established engine
    # eng_id = matlab.engine.find_matlab()
    # eng = matlab.engine.connect_matlab(eng_id[0])
    eng.cd(r'C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\matlab_code', nargout=0)
    # Important to specify the output
    c, d = eng.eval('main_function', nargout=2)
    print(c)
    self.a_bar = np.array(c)
    eng.quit()

# def control_barrier_function(self):
#
#
#
# def build_optimal_controller(self):


