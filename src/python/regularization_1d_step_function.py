#
# Based on document from M. Borsboom (update September 2022)
# To smooth or not to smooth
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
import PyQt6

matplotlib.use('QtAgg')


def cm2inch(cm):
    return cm / 2.54


def thomas_algorithm_5(a, b, c, d, e, f):
    n = len(f)

    a[0] = 0.0
    b[0] = 0.0
    d[0] = d[0] / c[0]
    e[0] = e[0] / c[0]
    f[0] = f[0] / c[0]
    c[0] = c[0] / c[0]

    for i in range(1, n-1):
        d[i  ] = (d[i]-b[i]*e[i-1])/(c[i]-b[i]*d[i-1])
        e[i  ] = (e[i]            )/(c[i]-b[i]*d[i-1])
        f[i  ] = (f[i]-b[i]*f[i-1])/(c[i]-b[i]*d[i-1])
        c[i  ] =  c[i]-b[i]*d[i-1]

        b[i+1] = b[i+1]-a[i+1]*d[i-1]
        c[i+1] = c[i+1]-a[i+1]*e[i-1]
        f[i+1] = f[i+1]-a[i+1]*f[i-1]

    i    = n-1
    f[i] = (f[i]-b[i]*f[i-1])/(c[i]-b[i]*d[i-1])
    c[i] =  c[i]-b[i]*d[i-1]
#-----------------------------------------------------------------------
#     back sweep
#-----------------------------------------------------------------------
    i=n-2
    f[i]=f[i]-d[i]*f[i+1]
    for i in range(n - 3, -1, -1):
        f[i] = f[i]-d[i]*f[i+1]-e[i]*f[i+2]

    return f


def main(Lx=2000., dx=100., c_psi = 0.125):
    # paragraph after eq. 10 of article
    #Lx = 100. * dx
    Lx += dx
    imax = int(Lx / dx) + 1
    bathymetry = 4
    Psi = c_psi * dx * dx

    step_left = 0.
    step_right = 1.
    #step_right = 2.0
    step_height = step_right - step_left

    refine  = 100
    x_ana = np.zeros(refine*(imax-1)+1, dtype=np.float64)
    ugiv_ana = np.zeros(refine*(imax-1)+1, dtype=np.float64)
    utilde_ana = np.zeros(refine*(imax-1)+1, dtype=np.float64)
    ubar_ana = np.zeros(refine*(imax-1)+1, dtype=np.float64)
    x = np.zeros(imax, dtype=np.float64)
    ugiv = np.zeros(imax, dtype=np.float64)
    ubar = np.zeros(imax, dtype=np.float64)
    ugrid = np.zeros(imax, dtype=np.float64)  # needed to draw grid points
    u0 = np.zeros(imax, dtype=np.float64)
    u0_xx = np.zeros(imax, dtype=np.float64)
    u1_xx = np.zeros(imax, dtype=np.float64)
    diff_abs = np.zeros(imax, dtype=np.float64)
    psi = np.zeros(imax, dtype=np.float64)
    a = np.zeros(imax, dtype=np.float64)
    b = np.ones(imax, dtype=np.float64)
    c = np.zeros(imax, dtype=np.float64)
    d = np.zeros(imax, dtype=np.float64)
    e = np.zeros(imax, dtype=np.float64)
    f = np.zeros(imax, dtype=np.float64)
    err = np.zeros(imax, dtype=np.float64)

    for i in range(0, refine*(imax-1) + 1):
        x_ana[i] = i/(refine*(imax-1)) * Lx
    for i in range(0, imax):
        x[i] = i/(imax-1) * Lx

    xc = 1./2. * Lx  # - 0.50 * dx
    jc = int(xc/dx)
    step_at = 0.00
    xc = x[jc] + step_at * dx

    for i in range(0, refine*(imax-1) + 1):
        if bathymetry == 4:
            ugiv_ana[i] = step_left
            utilde_ana[i] = step_left + 0.5 * step_height * np.exp((x_ana[i] - xc)/np.sqrt(Psi))
            if x_ana[i] > xc:
                ugiv_ana[i] = step_right
                utilde_ana[i] = step_right - 0.5 * step_height * np.exp((xc - x_ana[i])/np.sqrt(Psi))
    for i in range(0, imax):
        if bathymetry == 4:
            ugiv[i] = step_left
            ubar[i] = step_left + 0.5 * step_height * np.exp((x[i] - xc)/np.sqrt(Psi))
            if x[i] > xc:
                ugiv[i] = step_right
                ubar[i] = step_right - 0.5 * step_height * np.exp((xc - x[i])/np.sqrt(Psi))
        u0[i] = ugiv[i]
    for j in range(0, imax-1):
        for i in range(0, refine):
            k = j*refine + i
            fac = float(i)/(float(refine))
            ubar_ana[k] = (1.0 - fac) * ubar[j] + fac * ubar[j+1]
    ubar_ana[refine*(imax-1)] = ubar[imax-1]

    iter_max = 1000
    diff_max0 = 0
    for it in range(0, iter_max):
        for i in range(1, imax - 1):
            u0_xx[i] = (u0[i + 1] - 2. * u0[i] + u0[i - 1]) # / (dx * dx)
        u0_xx[0] = u0_xx[1]
        u0_xx[imax - 1] = u0_xx[imax - 2]

        c_error = c_psi

        for i in range(1, imax - 1):
            a[i] = 0.0
            a[i] = 1./8. - c_error
            b[i] = 6./8. + 2. * c_error
            c[i] = 1./8. - c_error
            e[i] = 0.0
            f[i] = np.abs(u0_xx[i])
        i = 0
        a[i] = 0.0
        b[i] = 0.0
        c[i] = 1.0
        d[i] = -1.0
        e[i] = 0.0
        f[0] = 0.0
        i = imax - 1
        a[i] = 0.0
        b[i] = -1.0
        c[i] = 1.0
        d[i] = 0.0
        e[i] = 0.0
        f[i] = 0.0
        Err = thomas_algorithm_5(a, b, c, d, e, f)

        max_Err = 0.0
        for i in range(0,imax):
            max_Err = max(max_Err, Err[i])

        for i in range(0, imax):
            psi[i] = c_psi * dx * dx * Err[i] # /(max_Err + 1e-12)

        for i in range(1, imax - 1):
            psi_1 = 0.5 * (psi[i - 1] + psi[i])
            psi_2 = 0.5 * (psi[i] + psi[i + 1])
            a[i] = 0.0
            b[i] = 1./8. * dx - psi_1 / dx
            c[i] = 6./8. * dx + psi_1 / dx + psi_2 / dx
            d[i] = 1./8. * dx - psi_2 / dx
            e[i] = 0.0
            f[i] = 1./8. * ugiv[i - 1] + 6./8. * ugiv[i] + 1./8. * ugiv[i + 1]
        i = 0
        a[i] = 0.0
        b[i] = 0.0
        c[i] = 1.0
        d[i] = -1.0
        e[i] = 0.0
        f[0] = 0.0
        i = imax - 1
        a[i] = 0.0
        b[i] = -1.0
        c[i] = 1.0
        d[i] = 0.0
        e[i] = 0.0
        f[i] = 0.0
        u1 = thomas_algorithm_5(a, b, c, d, e, f)

        u1_xx_max = 0.
        for i in range(1, imax - 1):
            u1_xx[i] = (u1[i + 1] - 2. * u1[i] + u1[i - 1]) / (dx * dx)
            u1_xx_max = max(u1_xx_max, np.abs(u1_xx[i]))
        diff_max1 = 0.0
        for i in range(0, imax):
            diff_abs[i] = np.abs(u1[i] - ugiv[i]) - 25.
            diff_max1 = max(diff_max1, np.abs(u1[i] - ugiv[i]))
        # print("u1 - ugiv:", diff_max1)
        diff_max10 = abs(diff_max1 - diff_max0)
        print("diff:", diff_max10)
        diff_max0 = diff_max1
        u0 = u1

        if it+1 == iter_max and diff_max10 >= 1e-8:
            print("Not converged: ", it+1, diff_max10)

        if diff_max10 < 1e-8 or it+1==iter_max:

            # -------------------------------------------------------------------
            # plot the series
            fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))
            # plt.ylim(-2, 2)
            # plt.ylim(0, 4)
            # ax1.set_title('$u_t + u_x = 0$')
            ax1.set_ylabel('Amplitude [m]  $\\longrightarrow$')  # after definition of spines
            ax1.set_xlabel('x [m] $\\longrightarrow$')  # after definition of spines

            ax1.xaxis.set_major_locator(MultipleLocator(dx))
            ax1.xaxis.set_minor_locator(MultipleLocator(dx/2))
            ax1.yaxis.set_major_locator(MultipleLocator(step_height/10))

            major_grid_color = 'c'
            minor_grid_color = 'g'
            # ax1.xaxis.grid(True, 'major', linestyle='-', linewidth=0.25, color=major_grid_color)
            ax1.xaxis.grid(True, 'minor', linestyle='-', linewidth=0.99, color=minor_grid_color)
            #x_offset = min(5.0 * dx, Lx/4)
            #x_offset = min(10.0 * dx, Lx/2)
            #ax1.set_xlim([0.5 * Lx - x_offset, 0.5 * Lx + x_offset])
            #ax1.set_xlim([-Lx/2, Lx/2])
            ax1.set_xlim([500, 1500])
            ax1.set_ylim([step_left - 0.05*step_height, step_right + 0.05*step_height])

            y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ax1.yaxis.set_major_formatter(y_formatter)

            tekst1a = 'Given function (analytic)'
            tekst1b = 'Given function (on grid)'
            tekst2 = 'Smoothed $\widetilde{u}(x)$'
            tekst3 = 'Piecewise lin. $\overline{u}(x)$'
            tekst4 = 'Smoothed with error estimation'
            ax1.step(x_ana, ugiv_ana, '-', color='b', label=tekst1a, linewidth=2.0)
            ax1.plot(x_ana, utilde_ana, ':', color='b', label=tekst2)
            # ax1.plot(x, ugiv, '-', color='c', marker = 'o', label=tekst1b)
            ax1.plot(x, ubar, '-.', color='r', marker = 'o', label=tekst3)
            # ax1.plot(x, u0, '-.', color='g', marker = 'o', label=tekst4)
            tekst5 = 'Grid points'
            ax1.plot(x, ugrid, '.', color='black', marker = 'o', markerfacecolor='none', label=tekst5)

            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels, loc='upper left')

            tekst_color = '#000001'
            tekst = ('To smooth or not to smooth')
            #fig.text(0.125, 0.95, tekst, fontsize=10, color=tekst_color)
            tekst = ('$c_{\Psi}$ = %.3f; $\Delta x$ = %.1f [m]; ${\Psi = c_{\Psi}\Delta x^2}$ = %.1f' % (c_psi, dx, Psi))
            tekst_color = 'b'
            if it+1 == iter_max:
                tekst_color = 'r'
            fig.text(0.125, 0.90, tekst, fontsize=10, color=tekst_color)

            # tekst = ('Max $|\overline{u}_{xx}|$ = %.5e' % (u1_xx_max))
            # fig.text(0.75, 0.90, tekst, fontsize=10, color=tekst_color)

            fig.set_size_inches(cm2inch(35.0), cm2inch(12.5))
            figManager = plt.get_current_fig_manager()
            # figManager.window.showMaximized()
            if not os.path.exists('figures'):
                os.mkdir('figures')
            if not os.path.exists('data'):
                os.mkdir('data')

            text = ('figures/regul_1d_step_function_dx%s_cpsi%s.pdf' % (str(dx), str(c_psi)))
            plt.savefig(text, format='pdf')
            fig.canvas.manager.set_window_title('Regularization_1d_step_function')

            text = ('data/regul_1d_step_function_dx%s_cpsi%s.tek' % (str(dx), str(c_psi)))
            with open(text, "w") as logfile:
                #for j in range(0, len(c_psi)):
                    logfile.write('*\n')
                    logfile.write('* Step function\n')
                    logfile.write('*\n')
                    logfile.write('* column 1: x-coordinate\n')
                    logfile.write('* column 2: u_giv_ana\n')
                    logfile.write('* column 3: u_tilde_ana\n')
                    logfile.write('* column 4: u_bar_ana\n')
                    logfile.write('* column 5: u_tilde_ana - u_giv_ana\n')
                    logfile.write('* column 6: u_bar_ana - u_giv_ana\n')
                    logfile.write('* column 7: u_bar_ana - u_tilde_ana\n')
                    logfile.write('*\n')
                    text = ('* cpsi_%s\n' % str(c_psi))
                    logfile.write(text)
                    logfile.write('*\n')
                    logfile.write('analytic_solutions\n')
                    tekst = ('%.d %.d\n' % (refine * (imax - 1) + 1, 7))
                    logfile.write(tekst)
                    for i in range(0, refine * (imax - 1) + 1):
                        tekst = ('%8.8f  %8.8f  %8.8f  %8.8f  %8.8f  %8.8f  %8.8f\n' % (
                            x_ana[i], ugiv_ana[i], utilde_ana[i], ubar_ana[i],
                            abs(utilde_ana[i] - ugiv_ana[i]),
                            abs(ubar_ana[i] - ugiv_ana[i]),
                            abs(ubar_ana[i] - utilde_ana[i] )))
                        logfile.write(tekst)
                    logfile.close()

            plt.show()
            break

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
