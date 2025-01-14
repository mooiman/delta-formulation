#
# From M. Borsboom
# Applied Numerical Mathematics 26 (1998)
# Borsboom_developerrorminadaptgridmethod_ApplNumerMath1998.pdf
#
# Programmer: Jan Mooiman
#

import os  # file exists
import sys  # system
from datetime import datetime, date, time, timedelta
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PyQt6

matplotlib.use('QtAgg')


def cm2inch(cm):
    return cm / 2.54


def thomas_algorithm(a, b, c, d):
    n = len(d)
    c_star = np.zeros(n, float)
    d_star = np.zeros(n, float)
    f = np.zeros(n, float)

    c_star[0] = c[0] / b[0]
    d_star[0] = d[0] / b[0]

    for i in range(1, n):
        m = 1.0 / (b[i] - a[i] * c_star[i - 1])
        c_star[i] = c[i] * m
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m

    f[n - 1] = d_star[n - 1]
    for i in range(n - 2, -1, -1):
        f[i] = d_star[i] - c_star[i] * f[i + 1]

    return f


def main(length=6000., dx=10):
    #length = 100. * dx
    nx = int(length / dx + 1)
    bathymetry = 1

    x = np.zeros(nx, dtype=np.float64)  # x-coordinate
    y = np.zeros(nx, dtype=np.float64)  # y-coordinate
    s = np.zeros(nx, dtype=np.float64)  # initial water level
    u = np.zeros(nx, dtype=np.float64)  # initial water level
    v = np.zeros(nx, dtype=np.float64)  # initial water level

    a0 = 0.001  # Amplitude at boundary
    for i in range(0, nx):
        x[i] = i * dx - length/2.
        s[i] = 2.0 * a0 * np.exp(-x[i] * x[i] / (500. * 500.))


    # -------------------------------------------------------------------
    # plot the series
    fig, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(1, 1))
    # plt.ylim(-2, 2)
    # plt.ylim(0, 4)
    # ax1[0].set_title('$u_t + u_x = 0$')
    ax1[0].set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
    ax1[0].set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

    major_grid_color = 'c'
    ax1[0].xaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
    ax1[0].yaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
    ax1[0].set_ylim([-0.0005, 0.0025])

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1[0].yaxis.set_major_formatter(y_formatter)

    tekst1 = 'Initial water level'
    ax1[0].plot(x, s, '-', color='b', label=tekst1)

    handles, labels = ax1[0].get_legend_handles_labels()
    ax1[0].legend(handles, labels, bbox_to_anchor=(0.85, 0.80), loc='lower left')

    fig.set_size_inches(cm2inch(35.0), cm2inch(17.5))
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    if not os.path.exists('figures'):
        os.mkdir('figures')
    if not os.path.exists('data'):
        os.mkdir('data')

    plt.show()

    with open("data/initial_water_level_delft3dflow.ini", "w") as logfile:
        # Delft3D-FLOW format
        # water level
        for j in range(0,3):
            for i in range(0, nx):
                tekst = ('%.15f  ' % (s[i]))
                logfile.write( tekst )
            tekst = ('%.15f\n' % (s[nx-1]))
            logfile.write(tekst)
        # u-velocity
        for j in range(0,3):
            for i in range(0, nx):
                tekst = ('%.15f  ' % (u[i]))
                logfile.write( tekst )
            tekst = ('%.15f\n' % (u[nx-1]))
            logfile.write(tekst)
        # v-velocity
        for j in range(0,3):
            for i in range(0, nx):
                tekst = ('%.15f  ' % (v[i]))
                logfile.write( tekst )
            tekst = ('%.15f\n' % (v[nx-1]))
            logfile.write(tekst)
        logfile.close()
    with open("data/initial_water_level_dflowfm.xyz", "w") as logfile:
        # D-FLow FM format
        for i in range(0, nx):
            tekst = ('%.4f %.4f %.15f\n' % (x[i], y[i], s[i]))
            logfile.write( tekst )
        logfile.close()

    return 0


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    start_time = datetime.now()
    print('\nStart: %s' % start_time)

    rtn_value = main()

    print('\nStart: %s' % start_time)
    print('End  : %s' % datetime.now())

    if rtn_value != 0:
        sys.exit(1)

    sys.exit(0)
