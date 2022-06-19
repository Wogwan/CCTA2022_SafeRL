#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last update: 2022-06-17
Name: cbf_sosp.py
Description: API of programming controllers
"""

import numpy as np

def build_barrier(self, eng):
    eng.cd(r'.\matlab_code', nargout=0)
    c, d = eng.eval('main_function', nargout=2)
    print(c)
    self.a_bar = np.array(c)
    eng.quit()


