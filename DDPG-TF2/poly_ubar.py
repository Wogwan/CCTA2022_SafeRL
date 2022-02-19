import numpy as np
from sympy import *


def poly_ubar_match(b, s_input):
    x1 = np.arcsin(s_input[0, 0])
    x2 = s_input[0, 3]
    fac = [x1**2, x1*x2, x1, x2**2, x2, 1]
    fun = 0
    for k in range(len(b)):
        fun = fac[k]*b[k]+fun

# def calculate_ubar(b, s_input, a_input):
#     k = poly_ubar_match(b)

