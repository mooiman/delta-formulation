#
# Programmer: Jan Mooiman
# delta_formulation_content.pdf
# equation B.18: 1/24 (kdx)^2 < 0.01 (i.e 1 percent)
#

import os  # file exists
import sys  # system
from datetime import datetime, date, time, timedelta
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PyQt6

def main():
    nx = 4
    acc = np.zeros(nx, dtype=np.float64)
    c_psi = np.zeros(nx, dtype=np.float64)
    k2dx2 = np.zeros(nx, dtype=np.float64)  # (kdx)^2
    N = np.zeros(nx, dtype=np.float64)  # (kdx)^2
    c_error = np.zeros(nx, dtype=np.float64)  # (kdx)^2

    acc[0] = 0.05
    acc[1] = 0.02
    acc[2] = 0.01
    acc[3] = 0.005

    # equation B.18: 1/24 (kdx)^2 < 0.01 (i.e 1 percent)
    for i in range(0,nx):
        k2dx2[i] = 24.0 * acc[i]
        c_psi[i] = 1./k2dx2[i]
        N[i] = (2.0 * np.pi)/np.sqrt(k2dx2[i])
        c_error[i] = c_psi[i] * np.pi/2. * 1./k2dx2[i]


    for i in range(0,nx):
        print('acc [%]; c_psi; c_error; N; (kdx)^2')
        print(f"{acc[i]*100.:.1f}")
        print(f"{c_psi[i]:.4f}")
        print(f"{c_error[i]:.4f}")
        print(f"{N[i]:.4f}")
        print(f"{k2dx2[i]:.4f}")
        print()

    return acc
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    start_time = datetime.now()
    print('\nStart: %s' % start_time)

    rtn_value = main()

    print('\nStart: %s' % start_time)
    print('End  : %s' % datetime.now())

    sys.exit(0)
