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


def main(Lx=6000., dx=10):
    #length = 100. * dx
    nx = int(Lx / dx) + 1 + 2  # two extra virtual points

    x = np.zeros(nx, dtype=np.float64)  # x-coordinate
    y = np.zeros(nx, dtype=np.float64)  # y-coordinate
    s = np.zeros(nx, dtype=np.float64)  # initial water level

    sigma = 350.0
    a0 = 0.001  # Amplitude at boundary
    for i in range(0, nx):
        x[i] = i * dx - Lx/2 - dx
        s[i] = 2.0 * a0 * np.exp(-x[i] * x[i] / (2. * sigma * sigma))

    srange = abs(max(s) - min(s))
    if srange < 0.01:
        srange = 0.01
    # -------------------------------------------------------------------
    # plot the series
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(1, 1))
    fig.set_size_inches(cm2inch(35.0), cm2inch(12.5))

    ax1.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
    ax1.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

    major_grid_color = 'c'
    ax1.xaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
    ax1.yaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
    ax1.set_ylim(min(s) - 0.01*srange, max(s) + 0.01 * srange)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)

    tekst1 = 'Initial water level'
    ax1.plot(x, s, '-', color='b', label=tekst1)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, bbox_to_anchor=(0.80, 0.85), loc='lower left')

    #figManager = plt.get_current_fig_manager()
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
