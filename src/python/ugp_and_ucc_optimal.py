#
# Initialize a quadratic function  for the FVE method
#
# Programmer: Jan Mooiman
#

import os  # file exists
import sys  # system
from datetime import datetime, date, time, timedelta
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
import PyQt6

matplotlib.use('QtAgg')

def u(x, A, B, C, D, x1):
    return (
        A
        + C * x
        + ((3 * (B - A) - (2 * C + D) * x1) / x1**2) * x**2
        + (((C + D) * x1 - 2 * (B - A)) / x1**3) * x**3
    )

def fgp(xi, c0, c1, c2):
    # parabola through three grid points
    return (
        c0 * (1. - xi) +
        c1 * xi +
        0.5 * ( c0 - 2. * c1 + c2) * (xi - 1.) * xi
    )
def fcc(xi, c0, c1, c2):
    # parabola through cell centres
    return (
        c0 * (1. - xi) +
        c1 * xi +
        0.5 * ( c0 - 2. * c1 + c2) * (xi - 0.5)  * (xi - 0.5)
    )

def cm2inch(cm):
    return cm / 2.54

def main(Lx=12., dx=2.0):
    nx = int(Lx / dx) + 1 + 2  # two extra virtual points

    mass = np.zeros(3, dtype=np.float64)
    alpha = 0.125
    mass[0] = alpha
    mass[1] = 1. - 2. * alpha
    mass[2] = alpha

    refine = 1024
    x_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    u_gp = np.full(refine*(nx-1) + 1, np.nan, dtype=np.float64)
    u_cc = np.full(refine*(nx-1) + 1, np.nan, dtype=np.float64)
    u_opt = np.full(refine*(nx-1) + 1, np.nan, dtype=np.float64)
    x = np.zeros(nx, dtype=np.float64)
    ubar = np.zeros(nx, dtype=np.float64)

    for i in range(0, refine*(nx-1)+1):
        x_ana[i] = float(i) * dx/float(refine) - dx
    for i in range(0, nx):
        x[i] = float(i - 1) * dx  # needed for ubar

    ubar[0] = 5.0
    ubar[1] = 4.5
    ubar[2] = 7.0
    ubar[3] = 3.0
    ubar[4] = 4.0
    if (nx >= 6): ubar[5] = 0.0
    if (nx >= 7): ubar[6] = 5.0
    if (nx >= 8): ubar[7] = 3.0
    if (nx >= 9): ubar[8] = 2.0
    if (nx >= 12):
        ubar[9] = 2.0
        ubar[10] = 3.0
        ubar[11] = 5.0
        ubar[12] = 3.0
        ubar[13] = 5.0
        ubar[14] = 3.0
    #
    # Draw a parabola through 3 grid points over a control volume
    #
    for i in range(0, nx - 2):
        c0 = ubar[i]
        c1 = ubar[i+1]
        c2 = ubar[i+2]
        for j in range(0, refine):
            k = (i) * refine + j + int(refine/2)
            xi = 0.5 + float(j)/float(refine)
            utmp = fgp(xi, c0, c1, c2)
            u_gp[k] = utmp

    # for each control volume (ie from cell centre to cell centre) compute the spline with given
    # values and gradient at the control volume faces    for i in range(0, nx - 2):
    for i in range(0, nx - 2):
        c0 = ubar[i]
        c1 = ubar[i+1]
        c2 = ubar[i+2]
        for j in range(0, refine):
            k = (i) * refine + j + int(refine/2)
            xi = 0.5 + float(j)/float(refine)
            utmp = fcc(xi, c0, c1, c2)
            u_cc[k] = utmp

    for i in range(0, refine * (nx - 1) + 1):
        gamma = 2./3.
        u_opt[i] = (1. - gamma) * u_gp[i] + gamma * u_cc[i]
    # ==================================================================
    # plot the series
    # fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))

    #fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))
    # plt.ylim(-2, 2)
    # plt.ylim(0, 4)
    # ax[0].set_title('$u_t + u_x = 0$')
    ax1.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
    ax1.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

    major_grid_color = 'black'
    ax1.xaxis.set_ticks(np.arange(-dx, Lx+2.* dx, dx))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.xaxis.grid(True, 'major', linestyle='-', linewidth=0.0, color=major_grid_color)
    ax1.xaxis.grid(True, 'minor', linestyle='-', linewidth=1.0, color='g')
    ax1.yaxis.grid(True, 'major', linestyle='-', linewidth=0.15, color=major_grid_color)

    infty = 1e12
    x_min = infty
    x_max = -infty
    u_min = infty
    u_max = - infty
    for i in range(0, nx):
        x_min = min(x_min, x[i])
        x_max = max(x_max, x[i])
        u_min = min(u_min, ubar[i])
        u_max = max(u_max, ubar[i])
    x_threshold = (x_max - x_min)* 0.025
    u_threshold = (u_max - u_min)* 0.075
    ax1.set_xlim([x_min - x_threshold, x_max + x_threshold])
    ax1.set_ylim([u_min - u_threshold, u_max + u_threshold])
    ax1.axvspan(x_min - x_threshold, x_min + 0.5*dx, color='gray', alpha=0.5)
    ax1.axvspan(x_max + x_threshold, x_max - 0.5*dx, color='gray', alpha=0.5)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)

    tekst2 = '$\overline{u}(x)$'
    tekst1 = '$u_{gp}, \gamma = 0$'
    tekst3 = '$u_{cc}, \gamma = 1$'
    tekst4 = '$u_{opt}, \gamma = 2/3$'

    ax1.plot(x, ubar, '-', color='b', marker = 'o', label=tekst2)
    ax1.plot(x_ana, u_gp, '--', color='0', markersize=2, label=tekst1)
    ax1.plot(x_ana, u_cc, ':', color='0', markersize=2, label=tekst3)
    ax1.plot(x_ana, u_opt, '-', color='0', linewidth=2, label=tekst4)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='lower left')

    fig1.set_size_inches(cm2inch(25.0), cm2inch(12.5))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()

    if not os.path.exists('figures'):
        os.mkdir('figures')
    text = ('figures/ugp_and_ucc_optimal_lx=%s_dx=%s.pdf' % (str(Lx), str(dx)))
    plt.savefig(text, format='pdf')
    fig1.canvas.manager.set_window_title('Optimal second order function approximation')
    plt.show()

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
