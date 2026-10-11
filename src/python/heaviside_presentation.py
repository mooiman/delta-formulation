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
from matplotlib.ticker import MultipleLocator
import numpy as np
import regularization_function as regf

matplotlib.use('QtAgg')


def cm2inch(cm):
    return cm / 2.54


def main(bath_in = 5, Lx_in=1000., dx_in=100.0, c_psi_in= 2., left_in = 0.0, right_in = 1.0):  # c_psi paragraph after eq. 10 of article
    bathymetry = int(bath_in)
    Lx = float(Lx_in)
    dx = float(dx_in)
    c_psi = float(c_psi_in)
    step_left = float(left_in)
    step_right = float(right_in)

    if (int(Lx/dx) != Lx/dx):
        print("Grid space is not a integer divider of domain length: Lx/dx = ", Lx/dx)
        return(1)

    nx = int(Lx / dx) + 1 + 2  # two extra virtual points
    # bathymetry:
    # 5: Step function step_left and step_right

    refine = 256
    x_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    ugiv_ana = np.zeros(refine * (nx - 1) + 1, dtype=np.float64)
    ubar_ana = np.zeros(refine * (nx - 1) + 1, dtype=np.float64)

    x = np.zeros(nx, dtype=np.float64)
    ugrid = np.zeros(nx, dtype=np.float64)
    ubar = np.zeros(nx, dtype=np.float64)
    ugiv = np.zeros(nx, dtype=np.float64)

    for i in range(0, refine*(nx-1)+1):
        x_ana[i] = i/(refine*(nx-3)) * Lx - dx
    for i in range(0, nx):
        x[i] = i/(nx-3) * Lx - dx

    bathymetry_desc = "Not (yet) specified"
    if bathymetry == 5:
        bathymetry_desc = "Single step (half way)"
        step_height = step_right - step_left

        xc = 1. / 2. * Lx
        jc = int(xc / dx) + 1
        step_at = 0.0
        xc = x[jc] + step_at * dx

        for i in range(0, nx):
            ugiv[i] = step_left
            if x[i] > xc:
                ugiv[i] = step_right
        for i in range(0, refine * (nx - 1) + 1):
            ugiv_ana[i] = step_left
            if x_ana[i] > xc:
                ugiv_ana[i] = step_right
    else:
        print("No valid bathymetry option defined, value '%s' is not supported." % bathymetry)
        return(1)

    #---------------------------------------------------------------------------------------
    #utilde = regf.get_ftilde(ugiv, ugiv_ana, c_psi, dx, nx, refine)  # taken from funcfitxxx.mw
    ubar, Err, psi, diff_max, _ = regf.compute_regularization(c_psi, ugiv, dx, nx, ugiv_ana, refine)
    #---------------------------------------------------------------------------------------

    if diff_max < 1e-5:
        for j in range(0, nx - 1):
            for i in range(0, refine):
                k = j * refine + i
                fac = float(i) / (float(refine))
                ubar_ana[k] = (1.0 - fac) * ubar[j] + fac * ubar[j + 1]
        ubar_ana[refine * (nx - 1)] = ubar[nx - 1]

        # ==================================================================
        # plot the series
        fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(1, 1))
        fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(1, 1))
        ax1.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
        ax1.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines
        ax2.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
        ax2.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

        if nx < 500:
            ax1.xaxis.set_major_locator(MultipleLocator(dx))
            ax1.xaxis.set_minor_locator(MultipleLocator(dx/2))
            ax2.xaxis.set_major_locator(MultipleLocator(dx))
            ax2.xaxis.set_minor_locator(MultipleLocator(dx/2))

        minor_grid_color = 'g'
        ax1.xaxis.grid(True, 'minor', linestyle='-', linewidth=0.99, color=minor_grid_color)
        ax2.xaxis.grid(True, 'minor', linestyle='-', linewidth=0.99, color=minor_grid_color)

        x_offset = 10. * dx
        if bathymetry == 5:
            ax1.set_xlim([0.5 * Lx - x_offset, 0.5 * Lx + x_offset])
            ax2.set_xlim([0.5 * Lx - x_offset, 0.5 * Lx + x_offset])
            ax1.set_xlim([0., 1000.])
            ax2.set_xlim([0., 1000.])
            ax1.set_ylim([step_left - 0.05 * step_height, step_right + 0.05 * step_height])
            ax2.set_ylim([step_left - 0.05 * step_height, step_right + 0.05 * step_height])

        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax1.yaxis.set_major_formatter(y_formatter)
        ax2.yaxis.set_major_formatter(y_formatter)

        tekst1 = 'Heaviside function'
        tekst5 = 'Grid points'
        ax1.plot(x_ana, ugiv_ana, '-', color='b', label=tekst1)
        ax1.plot(x, ugrid, '.', color='black', marker='o', markerfacecolor='none', label=tekst5)

        tekst2 = 'Regularized function'
        ax2.plot(x_ana, ubar_ana, '-', color='b', label=tekst2)
        ax2.plot(x, ubar, '.', color='b', )
        ax2.plot(x, ugrid, '.', color='black', marker='o', markerfacecolor='none', label=tekst5)

        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(handles1, labels1, prop={"size": 12}, loc='upper left')
        ax2.legend(handles2, labels2, prop={"size": 12}, loc='upper left')

        fig1.set_size_inches(cm2inch(35.0), cm2inch(12.5))
        fig2.set_size_inches(cm2inch(35.0), cm2inch(12.5))



        figManager = plt.get_current_fig_manager()
        # figManager.window.showMaximized()
        if not os.path.exists('figures'):
            os.mkdir('figures')

        text = ('figures/heaviside_step_at_%s.pdf' % (str(step_at)))
        fig1.savefig(text, format='pdf')
        text = ('figures/heaviside_step_at_%s_regularized.pdf' % (str(step_at)))
        fig2.savefig(text, format='pdf')

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
