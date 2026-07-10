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

def cm2inch(cm):
    return cm / 2.54

def main(Lx=12., dx=3.0):
    nx = int(Lx / dx) + 1 + 2  # two extra virtual points

    mass = np.zeros(3, dtype=np.float64)
    alpha = 0.125
    mass[0] = alpha
    mass[1] = 1. - 2. * alpha
    mass[2] = alpha

    refine = 16
    x_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    u_g0 = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    u_g1 = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    x = np.zeros(nx, dtype=np.float64)
    ubar = np.zeros(nx, dtype=np.float64)

    for i in range(0, refine*(nx-1)+1):
        x_ana[i] = float(i) * dx/float(refine) - dx
    for i in range(0, nx):
        x[i] = float(i - 1) * dx  # needed for ubar

    mid = Lx/2.
    quart = Lx/4.
    a_coef = 1.
    for i in range(0, refine*(nx-1)+1):
        if x_ana[i] < mid:
            u_g0[i] = a_coef*(x_ana[i] - mid + quart) ** 2.
        else:
            u_g0[i] = -a_coef*(x_ana[i] - mid - quart) ** 2. + 2. * a_coef*(x_ana[refine] - mid + quart)**2.
    for i in range(0, refine*(nx-1)+1):
        u_g1[i] = np.nan
    for i in range(0, nx):
        if x[i] < mid:
            ubar[i] = a_coef*(x[i] - mid + quart) **2.
            #ubar[i] = x[i]
        else:
            ubar[i] = -a_coef*(x[i] - mid - quart) ** 2. + 2. * a_coef*(x[1] - mid + quart) **2.
            #ubar[i] = x[i]
    # for each control volume (ie from cell centre to cell centre) compute the spline with given
    # values and gradient at the control volume faces
    for i in range(1, nx-1):
        a = 0.5 * (ubar[i] + ubar[i - 1])
        b = 0.5 * (ubar[i] + ubar[i + 1])
        da = (ubar[i] - ubar[i-1])/dx
        db = (ubar[i+1] - ubar[i])/dx
        xa = 0.5 * (x[i] + x[i-1])
        xb = 0.5 * (x[i] + x[i+1])
        # u(x)=(2t3−3t2+1)a+(−2t3+3t2)b+(t3−2t2+t)(2 da)+(t3−t2)(2 db)
        for j in range(0, refine):
            k = (i - 1) * refine + j + int(refine/2)
            t = float(j)/float(refine) * dx
            utmp = u(t, a, b, da, db, dx)
            u_g1[k] = utmp
    # -------------------------------------------------------------------
    # plot the series
    fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))
    # plt.ylim(-2, 2)
    # plt.ylim(0, 4)
    # ax[0].set_title('$u_t + u_x = 0$')
    ax2.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
    ax2.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

    major_grid_color = 'black'
    ax2.xaxis.set_ticks(np.arange(-dx, Lx+2.* dx, dx))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.xaxis.grid(True, 'major', linestyle='-', linewidth=0.0, color=major_grid_color)
    ax2.xaxis.grid(True, 'minor', linestyle='-', linewidth=1.0, color='g')
    ax2.yaxis.grid(True, 'major', linestyle='-', linewidth=0.15, color=major_grid_color)

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
    u_threshold = (u_max - u_min)* 0.05
    ax2.set_xlim([x_min - x_threshold, x_max + x_threshold])
    ax2.set_ylim([u_min - u_threshold, u_max + u_threshold])
    ax2.axvspan(x_min - x_threshold, x_min + 0.5*dx, color='gray', alpha=0.5)
    ax2.axvspan(x_max + x_threshold, x_max - 0.5*dx, color='gray', alpha=0.5)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.yaxis.set_major_formatter(y_formatter)

    tekst2 = '$\overline{u}(x)$'
    tekst1 = '$\widehat{u}(x), \gamma = 0$'
    tekst3 = '$\widehat{u}(x), \gamma = 1$'

    ax2.plot(x_ana, u_g0, '--', color='0', markersize=2, label=tekst1)
    ax2.plot(x_ana, u_g1, ':', color='0', markersize=2, label=tekst3)
    ax2.plot(x, ubar, '-', color='b', marker = 'o', label=tekst2)

    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='lower left')

    fig.set_size_inches(cm2inch(25.0), cm2inch(12.5))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()

    if not os.path.exists('figures'):
        os.mkdir('figures')
    text = ('figures/ugp_and_ucc_lx=%s_dx=%s.pdf' % (str(Lx), str(dx)))
    plt.savefig(text, format='pdf')
    fig.canvas.manager.set_window_title('Compatible initialization')
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
