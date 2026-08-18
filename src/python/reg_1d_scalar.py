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


def main(bath_in = 11, Lx_in=1000., dx_in=50.0, c_psi_in= 4.0, left_in = 0.0, right_in = 1.0):  # c_psi paragraph after eq. 10 of article
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
    # 0: tanh + step, as in the literature
    # 1: tanh with deeper step
    # 2: Experiment of Frank Platzek
    # 3: shoal: -10 [m] to -2.5 [m]
    # 4: Weir (Borsboom_FVmethoddesigned4erroranalysis_ProcFVCAIII2002.pdf)
    # 5: Step function step_left and step_right
    # 6: Constant value of 1.0
    # 7: Step at right boundary of 1.0
    # 8: Interface
    # 9: Summer-winterbed
    # 10: Wiggle
    # 11: Smooth step
    #

    Psi = c_psi * dx *dx

    refine = 256
    x_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    ugiv_ana = np.zeros(refine * (nx - 1) + 1, dtype=np.float64)
    ubar_ana = np.zeros(refine * (nx - 1) + 1, dtype=np.float64)

    x = np.zeros(nx, dtype=np.float64)
    utilde = np.zeros(nx, dtype=np.float64)
    ubar = np.zeros(nx, dtype=np.float64)
    ugiv = np.zeros(nx, dtype=np.float64)
    ugrid = np.zeros(nx, dtype=np.float64)
    ubar_xixi = np.zeros(nx, dtype=np.float64)
    ubar_xx = np.zeros(nx, dtype=np.float64)
    ugiv_xixi = np.zeros(nx, dtype=np.float64)
    ugiv_xx = np.zeros(nx, dtype=np.float64)
    diff_abs = np.zeros(nx, dtype=np.float64)
    Err = np.zeros(nx, dtype=np.float64)
    psi = np.zeros(nx, dtype=np.float64)

    for i in range(0, refine*(nx-1)+1):
        x_ana[i] = i/(refine*(nx-3)) * Lx - dx
    for i in range(0, nx):
        x[i] = i/(nx-3) * Lx - dx

    bathymetry_desc = "Not (yet) specified"
    if bathymetry == -1:
        bathymetry_desc = "0 to 2x with step at x = 0.5"
        for i in range(0, nx):
            ugiv[i] = 0.0
            if x[i] > 0.5 * Lx:
                ugiv[i] =  2. #  2. * x[i]
        for i in range(0, refine*(nx-1)+1):
            ugiv_ana[i] = 0.0
            if x_ana[i] > 0.5 * Lx:
                ugiv_ana[i] = 2.  # 2. * x_ana[i] / Lx
    elif bathymetry == 0:
        bathymetry_desc = "Standard tanh with step"
        for i in range(0, nx):
            ugiv[i] = 1.0
            if x[i] < 0.65 * Lx:
                ugiv[i] =  0.5 - 0.5 * np.tanh(20. * x[i] / Lx - 6.)
        for i in range(0, refine*(nx-1)+1):
            ugiv_ana[i] = 1.0
            if x_ana[i] < 0.65 * Lx:
                ugiv_ana[i] = 0.5 - 0.5 * np.tanh(20. * x_ana[i] / Lx - 6.)
        step_left = 0.
        step_right = 1.
        step_height = step_right - step_left
    elif bathymetry == 1:
        bathymetry_desc = "Higher standard tanh with step"
        for i in range(0, nx):
            ugiv[i] = 10.
            if x[i] < 0.65 * Lx:
                ugiv[i] = 10. * (0.5 - 0.5 * np.tanh(20. * x[i] / Lx - 6.))
        for i in range(0, refine * (nx - 1) + 1):
            ugiv_ana[i] = 10.0
            if x_ana[i] <= 0.65 * Lx:
                ugiv_ana[i] = 10.0 * (0.5 - 0.5 * np.tanh(20. * x_ana[i] / Lx - 6.))
        step_left = 0.
        step_right = 15.
        step_height = step_right - step_left
    elif bathymetry == 2:
        bathymetry_desc = "Experiment of Frank"
        step_left = 0.
        step_right = 25.
        step_height = step_right - step_left
        for i in range(0, nx):
            if x[i] < 0.3 * Lx:
                ugiv[i] = 20.0
            elif x[i] < 0.5 * Lx:
                ugiv[i] = 5.0
            elif x[i] < 0.7 * Lx:
                ugiv[i] = 0.0
            else:
                ugiv[i] = 20.0
        for i in range(0, refine * (nx - 1) + 1):
            if x_ana[i] < 0.3 * Lx - 0.5 * dx:
                ugiv_ana[i] = 20.0
            elif x_ana[i] < 0.5 * Lx - 0.5 * dx:
                ugiv_ana[i] = 5.0
            elif x_ana[i] < 0.7 * Lx - 0.5 * dx:
                ugiv_ana[i] = 0.0
            else:
                ugiv_ana[i] = 20.0
    elif bathymetry == 3:
        bathymetry_desc = "Shoaling bedlevel from -10 [m] to -2.5 [m]"
        slope_begin = 250.
        slope_end = 350.
        for i in range(0, nx):
            if x[i] < slope_begin:
                ugiv[i] = -10.0
            elif x[i] < slope_end:
                ugiv[i] = -10. + (x[i] - slope_begin) / (slope_end - slope_begin) * 7.5
            else:
                ugiv[i] = -2.5
        for i in range(0, refine * (nx - 1) + 1):
            if x_ana[i] < slope_begin:
                ugiv_ana[i] = -10.0
            elif x_ana[i] < slope_end:
                ugiv_ana[i] = -10. + (x_ana[i] - slope_begin) / (slope_end - slope_begin) * 7.5
            else:
                ugiv_ana[i] = -2.5
        step_left = -10.
        step_right = -2.5
        step_height = step_right - step_left
    elif bathymetry == 4:
        bathymetry_desc = "Weir: from -12 [m] to -5 [m] and from -5 [m] to -10 [m]"
        slope_up_begin = 200.
        slope_up_end = 250.
        slope_down_begin = 350.
        slope_down_end = 450.
        zb_begin = -12.0
        zb_weir = -5.0
        zb_end = -10.0
        step_left = -12.
        step_right = -5.
        step_height = step_right - step_left
        for i in range(0, nx):
            if x[i] < slope_up_begin:
                ugiv[i] = zb_begin
            elif x[i] < slope_up_end:
                ugiv[i] = zb_begin + (x[i] - slope_up_begin) / (slope_up_end - slope_up_begin) * (zb_weir - zb_begin)
            elif x[i] < slope_down_begin:
                ugiv[i] = zb_weir
            elif x[i] < slope_down_end:
                ugiv[i] = zb_weir + (x[i] - slope_down_begin) / (slope_down_end - slope_down_begin) * (zb_end - zb_weir)
            elif x[i] >= slope_down_end:
                ugiv[i] = zb_end
        for i in range(0, refine * (nx - 1) + 1):
            if x_ana[i] < slope_up_begin:
                ugiv_ana[i] = zb_begin
            elif x_ana[i] < slope_up_end:
                ugiv_ana[i] = zb_begin + (x_ana[i] - slope_up_begin) / (slope_up_end - slope_up_begin) * (
                            zb_weir - zb_begin)
            elif x_ana[i] < slope_down_begin:
                ugiv_ana[i] = zb_weir
            elif x_ana[i] < slope_down_end:
                ugiv_ana[i] = zb_weir + (x_ana[i] - slope_down_begin) / (slope_down_end - slope_down_begin) * (
                            zb_end - zb_weir)
            elif x_ana[i] >= slope_down_end - 0.5 * dx:
                ugiv_ana[i] = zb_end
    elif bathymetry == 5:
        bathymetry_desc = "Single step (half way)"
        step_height = step_right - step_left

        xc = 1. / 2. * Lx  # - 0.50 * dx
        jc = int(xc / dx)
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
    elif bathymetry == 6:
        bathymetry_desc = "Constant value"
        step_left = 0.
        step_right = 2.
        step_height = step_right - step_left
        for i in range(0, nx):
            ugiv[i] = 1.5
        for i in range(0, refine * (nx - 1) + 1):
            ugiv_ana[i] = 1.5
    elif bathymetry == 7:
        bathymetry_desc = "Single step (at right)"
        step_left = 0.0
        step_right = 1.0
        step_height = step_right - step_left
        for i in range(0, nx):
            if x[i] < 1.0 * Lx:
                ugiv[i] = step_left
            else:
                ugiv[i] = step_right
            ugiv[0] = 1.0
            ugiv[1] = 1.0
        for i in range(0, refine * (nx - 1) + 1):
            if x_ana[i] < 1.0 * Lx:
                ugiv_ana[i] = step_left
            else:
                ugiv_ana[i] = step_right
            for j in range(0, refine*(2-1) + 1):
                ugiv_ana[j] = 1.0
    elif bathymetry == 8:
        bathymetry_desc = "f(x)=0 x<0.5, f(x)=2x x>0.5: Interface problem"
        for i in range(0, nx):
            ugiv[i] = 0.0
            if x[i] > 0.5:
                ugiv[i] = 2. * x[i]
        for i in range(0, refine * (nx - 1) + 1):
            ugiv_ana[i] = 0.0
            if x_ana[i] > 0.5:
                ugiv_ana[i] = 2. * x_ana[i]
    elif bathymetry == 9:
        bathymetry_desc = "0.25 * Lx > x < 0.75 * Lx; f(x) = Winterbed"
        for i in range(0, nx):
            ugiv[i] = step_right  # summerbed
            if x[i] > 0.25 * Lx and x[i] < 0.75 * Lx:
                ugiv[i] = step_left   # winterbed
        for i in range(0, refine * (nx - 1) + 1):
            ugiv_ana[i] = step_right  # summerbed
            if x_ana[i] > 0.25 * Lx and x_ana[i] < 0.75 * Lx:
                ugiv_ana[i] = step_left   # winterbed
    elif bathymetry == 10:
        bathymetry_desc = "Wiggle (2 Dx)"
        for i in range(0, nx):
            ugiv[i] = (-1)**float(i)
        for i in range(1, nx):
            for j in range(0, refine):
                k = j + (i - 1) * (refine)
                alpha = j/(refine)
                ugiv_ana[k] = -(-1)**float(i) *np.cos(np.pi * float(j)/float(refine))
        ugiv_ana[refine * (nx - 1)] = ugiv[nx-1]
    elif bathymetry == 11:
        bathymetry_desc = "Smooth step"
        for i in range(1, nx):
            x_sc = (x[i] + dx)/Lx
            ugiv[i] = step_left + (step_right - step_left) * np.exp(-1./x_sc) / (np.exp(-1./x_sc) + np.exp(-1./((x[nx-1] + dx)/Lx - x_sc)))
        ugiv[0] = step_left
        for i in range(1, refine * (nx - 1) + 1):
            x_sc = (x_ana[i]+dx) / Lx
            ugiv_ana[i] = step_left + (step_right - step_left) * np.exp(-1./x_sc) / (np.exp(-1./x_sc) + np.exp(-1./((x[nx-1] + dx)/Lx - x_sc)))
        ugiv_ana[0] = step_left
    else:
        print("No valid bathymetry option defined, value '%s' is not supported." % bathymetry)
        return(1)

    #---------------------------------------------------------------------------------------
    utilde = regf.get_ftilde(ugiv, ugiv_ana, c_psi, dx, nx, refine)  # taken from funcfitxxx.mw
    ubar, Err, psi, diff_max, ugiv_xx_max = regf.compute_regularization(c_psi, ugiv, dx, nx, ugiv_ana, refine)
    #---------------------------------------------------------------------------------------
    if diff_max < 1e-5:
        for j in range(0, nx - 1):
            for i in range(0, refine):
                k = j * refine + i
                fac = float(i) / (float(refine))
                ubar_ana[k] = (1.0 - fac) * ubar[j] + fac * ubar[j + 1]
        ubar_ana[refine * (nx - 1)] = ubar[nx - 1]

        for i in range(1, nx - 1):
            ubar_xixi[i] = (ubar[i + 1] - 2. * ubar[i] + ubar[i - 1])
            ubar_xx[i] = (ubar[i + 1] - 2. * ubar[i] + ubar[i - 1])/(dx*dx)
            ugiv_xixi[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1])
            ugiv_xx[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1])/(dx*dx)

        max_ubar_xixi = 0.0
        max_ugiv_xixi = 0.0
        max_ubar_ugiv = 0.0
        for i in range(0, nx):
            diff_abs[i] = np.abs(ubar[i] - ugiv[i])
            max_ubar_xixi = max(max_ubar_xixi, ubar_xixi[i])
            max_ugiv_xixi = max(max_ugiv_xixi, ugiv_xixi[i])
            max_ubar_ugiv = max(max_ubar_ugiv, diff_abs[i])

        norm_ubar_ugiv = 0.0
        for i in range(0, refine * (nx - 1) + 1):
            norm_ubar_ugiv += abs(ubar_ana[i] - ugiv_ana[i])
        delta_ubar_ugiv = norm_ubar_ugiv / float(nx)

        # ==================================================================
        # plot the series
        fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(1, 1))
        fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(1, 1))
        # plt.ylim(-2, 2)
        # plt.ylim(0, 4)
        # ax1.set_title('$u_t + u_x = 0$')
        ax1.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
        ax1.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

        ubarrange = abs(max(ubar) - min(ubar))
        if ubarrange < 0.01:
            ubarrange = 0.01
            if bathymetry == 1:
                ubarrange = 20.
        ubar_xx_range = abs(max(ubar_xx) - min(ubar_xx))
        if ubar_xx_range < 0.01:
            ubar_xx_range = 0.01
        if nx < 500:
            ax1.xaxis.set_major_locator(MultipleLocator(dx))
            ax1.xaxis.set_minor_locator(MultipleLocator(dx/2))
            #ax1.yaxis.set_major_locator(MultipleLocator(ubarrange/10))

        minor_grid_color = 'g'
        major_grid_color = 'c'
        #ax1.xaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
        ax1.xaxis.grid(True, 'minor', linestyle='-', linewidth=0.99, color=minor_grid_color)
        #ax1.yaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)

        x_offset = 10. * dx
        if bathymetry == 0:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ax1.set_xlim([0, Lx])
            ax2.set_xlim([0, Lx])
            ax1.set_ylim([-0.1, 1.1])
        if bathymetry == 1:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ax1.set_xlim([0, Lx])
            ax2.set_xlim([0, Lx])
            ax1.set_ylim([-1, 11])
        if bathymetry == 2:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ax1.set_ylim([-1, 21.0])
        if bathymetry == 3:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ax1.set_ylim([-11., 1.0])
        if bathymetry == 4:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ax1.set_ylim([-12.5, 0.5])
        if bathymetry == 5:
            ax1.set_xlim([0.5 * Lx - x_offset, 0.5 * Lx + x_offset])
            ax2.set_xlim([0.5 * Lx - x_offset, 0.5 * Lx + x_offset])
            ax1.set_xlim([0., 1000.])
            ax2.set_xlim([0., 1000.])
            ax1.set_ylim([step_left - 0.05 * step_height, step_right + 0.05 * step_height])
        if bathymetry == 6:
            ax1.set_ylim([max(ubar) + ubarrange*0.05, min(ubar) - ubarrange*0.05])
            ax2.set_ylim([max(ubar_xx) + ubar_xx_range * 0.05, min(ubar_xx) - ubar_xx_range * 0.05])
        if bathymetry == 7:
            ax1.set_xlim([-dx, Lx+dx])
            ax2.set_xlim([-dx, Lx+dx])
            ymin1 = min(ubar)
            ymin2 = min(ugiv)
            ymax1 = max(ubar)
            ymax2 = max(ugiv)
            ax1.set_ylim([min(ymin1, ymin2) - ubarrange*0.05, max(ymax1, ymax2) + ubarrange*0.05])
            #ax2.set_ylim([max(ubar_xx) + ubar_xx_range * 0.05, min(ubar_xx) - ubar_xx_range * 0.05])
        if bathymetry == 10:
            ax1.set_xlim([0, Lx])
            ax2.set_xlim([0, Lx])
            # ax1.set_xlim([-dx, Lx+dx])
            # ax2.set_xlim([-dx, Lx+dx])
            ax1.set_ylim([-0.2, 0.2])
        if bathymetry == 11:
            ax1.set_xlim([0, Lx])
            ax2.set_xlim([0, Lx])
            ax1.set_ylim(step_left - 0.05*np.abs(step_left - step_right), step_right + 0.05 * np.abs(step_left - step_right))

        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax1.yaxis.set_major_formatter(y_formatter)

        tekst1 = 'Given function (analytic)'
        tekst2 = 'Regularized $\overline{u}(x)$'
        tekst3 = 'abs($\overline{u}-u_{giv}$)'
        tekst4 = '$\widetilde{u}(x_i)$'
#            if bathymetry == 0:
#                i_step = math.floor(0.65 * Lx/dx)
#                ax1.plot(x[0:i_step], ugiv[0:i_step], '-', color='b', label=tekst1)
#                ax1.step(x[i_step-1:], ugiv[i_step-1:], '-', color='b')
#            else:
#                 ax1.plot(x, ugiv, '-', color='b', label=tekst1)
        ax1.plot(x_ana, ugiv_ana, '-', color='b', label=tekst1)
        tekst1 = 'Given function (on grid)'
        ax1.plot(x, ugiv, '-', color='c', label=tekst1, marker = 'o')
        ax1.plot(x_ana, utilde, '.', color='black', markersize=1, label=tekst4)
        ax1.plot(x, ubar, '-', color='r', label=tekst2, marker = 'o')
        # if bathymetry == 1 or bathymetry == 2:
        #     ax1.plot(x, diff_abs, '--', color='b', label=tekst3)
        # else:
        #     ax1.plot(x, diff_abs, '--', color='b', label=tekst3)
        tekst5 = 'Grid points'
        ax1.plot(x, ugrid, '.', color='black', marker='o', markerfacecolor='none', label=tekst5)

        handles, labels = ax1.get_legend_handles_labels()
        if bathymetry == 0 or bathymetry == 1 or bathymetry == 11:
            ax1.legend(handles, labels, prop={"size": 12}, loc='center right')
        elif bathymetry == 9:
            ax1.legend(handles, labels, prop={"size": 12}, loc='lower right')
        else:
            ax1.legend(handles, labels, prop={"size": 12}, loc='upper right')

        ax2.set_ylabel('??? $\\longrightarrow$')  # after definition of spines
        ax2.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines
        ax2.xaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
        ax2.yaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
        # ax2.set_ylim([step_left - 0.05 * step_height, step_right + 0.05 * step_height])
        text1 = 'psi'
        text2 = 'ugiv_xixi (giv)'
        text3 = 'ubar_xixi (reg)'
        text4 = 'Error'
        # ax2.plot(x, psi, '-', color='magenta', label=text1)
        ax2.plot(x, np.abs(ugiv_xixi), '-', color='c', label=text2, marker = 'o')
        ax2.plot(x, np.abs(ubar_xixi), '-', color='r', label=text3, marker = 'o')
        tekst5 = 'Grid points'
        ax2.plot(x, ugrid, '.', color='black', marker='o', markerfacecolor='none', label=tekst5)

        handles, labels = ax2.get_legend_handles_labels()
        ax2.legend(handles, labels, prop={"size": 12}, loc='center right')

        tekst = ('$c_{\Psi}$ = %.1f; $\Delta x$ = %.3f [m]; ${\Psi = c_{\Psi}\Delta x^2}$ = %.3f' % (c_psi, dx, Psi))
        fig1.text(0.125, 0.90, tekst, fontsize=10)
        tekst = ('max|$\overline{u}$ - $u_{giv}$|= %.5e; Sum|$\overline{u}$ - $u_{giv}$|= %.5e' % (max_ubar_ugiv, delta_ubar_ugiv))
        fig1.text(0.125, 0.95, tekst, fontsize=10)

        tekst = ('max($ugiv_{xixi}$) = %.4f' % (max_ugiv_xixi))
        fig2.text(0.75, 0.93, tekst, fontsize=10)
        tekst = ('max($\overline{u}_{xixi}$) = %.4f' % (max_ubar_xixi))
        fig2.text(0.75, 0.90, tekst, fontsize=10)

        uxx_max = 0.
        for i in range(1, nx - 1):
            uxx_max = max(uxx_max, np.abs(ubar_xixi[i]))

        #tekst = ('max($u_{xixi}$) = %.5f' % (max_ubar_xixi))
        #fig1.text(0.1275, 0.35, tekst, fontsize=10)

        fig1.set_size_inches(cm2inch(35.0), cm2inch(12.5))
        fig2.set_size_inches(cm2inch(35.0), cm2inch(12.5))
        figManager = plt.get_current_fig_manager()
        # figManager.window.showMaximized()
        if not os.path.exists('figures'):
            os.mkdir('figures')
        if not os.path.exists('data'):
            os.mkdir('data')

        with open("data/bed_level_regularized.tek", "w") as logfile:
            logfile.write('* %s\n' % bathymetry_desc)
            logfile.write('* Regularized\n')
            logfile.write('* column 1: x-coordinate\n')
            logfile.write('* column 2: z_bed\n')
            logfile.write('bedlevel\n')
            tekst = ('%.d %.d\n' % (nx, 2))
            logfile.write(tekst)
            for i in range(0, nx):
                tekst = ('%.8f %.8f\n' % (x[i], ubar[i]))
                logfile.write( tekst )
            logfile.close()
        with open("data/bed_level_ugiv.tek", "w") as logfile:
            logfile.write('* %s\n' % bathymetry_desc)
            logfile.write('* u_given\n')
            logfile.write('* column 1: x-coordinate\n')
            logfile.write('* column 2: z_bed\n')
            logfile.write('bedlevel\n')
            tekst = ('%.d %.d\n' % (nx, 2))
            logfile.write(tekst)
            for i in range(0, nx):
                tekst = ('%.8f %.8f\n' % (x[i], ugiv[i]))
                logfile.write( tekst )
            logfile.close()

        text = ('figures/regul_1_1d_scalar_dx%s_cpsi%s.pdf' % (str(dx), str(c_psi)))
        fig1.savefig(text, format='pdf')
        text = ('figures/regul_2_1d_scalar_dx%s_cpsi%s.pdf' % (str(dx), str(c_psi)))
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
