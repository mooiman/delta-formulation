import os  # file exists
import sys  # system
from datetime import datetime, date, time, timedelta
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
from scipy.interpolate import CubicSpline
import PyQt6

def cm2inch(cm):
    return cm / 2.54

def main(Lx=12., dx=2.0):
    nx = int(Lx / dx) + 1 + 2  # two extra virtual points
    x_gp = np.zeros(nx, dtype=np.float64)
    u_gp = np.zeros(nx, dtype=np.float64)
    x_cc = np.zeros(nx-1, dtype=np.float64)
    u_cc = np.zeros(nx-1, dtype=np.float64)
    u_opt = np.zeros(500, dtype=np.float64)

    u_gp[0] = 5.0
    u_gp[1] = 4.5
    u_gp[2] = 7.0
    u_gp[3] = 3.0
    u_gp[4] = 4.0
    if (nx >= 6): u_gp[5] = 0.0
    if (nx >= 7): u_gp[6] = 5.0
    if (nx >= 8): u_gp[7] = 3.0
    if (nx >= 9): u_gp[8] = 2.0
    if (nx >= 12):
        u_gp[9] = 2.0
        u_gp[10] = 3.0
        u_gp[11] = 5.0
        u_gp[12] = 3.0
        u_gp[13] = 5.0
        u_gp[14] = 3.0

    for i in range(0,nx-1):
        u_cc[i] = (u_gp[i+1] + u_gp[i])/2.

    # x coordinates
    for i in range(0, nx):
        x_gp[i] = float(i - 1) * dx  # needed for u_gp
    for i in range(0, nx - 1):
        x_cc[i] = (x_gp[i+1] + x_gp[i])/2.

    # Create cubic spline
    cs_gp = CubicSpline(x_gp, u_gp)
    cs_cc = CubicSpline(x_cc, u_cc)

    # Evaluate on a finer grid
    x_gp_fine = np.linspace(x_gp.min(), x_gp.max(), 500)
    u_gp_fine = cs_gp(x_gp_fine)
    x_cc_fine = np.linspace(x_cc.min(), x_cc.max(), 500)
    u_cc_fine = cs_cc(x_cc_fine)

    x_opt_fine = np.linspace(x_gp.min(), x_gp.max(), 500)
    u_opt_fine = cs_cc(x_opt_fine)

    for i in range(0,500):
        u_opt[i] = 1./3. * u_gp_fine[i] + 2./3. * u_opt_fine[i]
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
        x_min = min(x_min, x_gp[i])
        x_max = max(x_max, x_gp[i])
    for i in range(0, len(x_gp_fine)):
        u_min = min(u_min, u_gp_fine[i])
        u_max = max(u_max, u_gp_fine[i])
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

    ax1.plot(x_gp, u_gp, '-', color='b', marker = 'o', label=tekst2)
    ax1.plot(x_gp_fine, u_gp_fine, '--', color='0', markersize=2, label=tekst1)
    ax1.plot(x_cc_fine, u_cc_fine, ':', color='0', markersize=2, label=tekst3)
    ax1.plot(x_gp_fine, u_opt, '-', color='0', linewidth=2, label=tekst3)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='lower left')

    fig1.set_size_inches(cm2inch(25.0), cm2inch(12.5))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()

    if not os.path.exists('figures'):
        os.mkdir('figures')
    text = ('figures/ugp_and_ucc_with_splines_lx=%s_dx=%s.pdf' % (str(Lx), str(dx)))
    plt.savefig(text, format='pdf')
    fig1.canvas.manager.set_window_title('Optimal second order function approximation')
    plt.show()


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    main()