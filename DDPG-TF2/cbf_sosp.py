import numpy
import os

import numpy as np
import scipy.io
import matlab
import matlab.engine

def build_barrier(self, eng):
    eng.cd(r'.\matlab_code', nargout=0)
    c, d = eng.eval('main_function', nargout=2)
    print(c)
    self.a_bar = np.array(c)
    eng.quit()


