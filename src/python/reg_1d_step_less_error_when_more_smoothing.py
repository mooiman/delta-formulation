#
# Based on document from M. Borsboom (update September 2022)
# To smooth or not to smooth
#
# More smoothing less difference between u-tilde and u-bar, also less dependent on location of the discontinuity
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
# import pyqt
from regularization_function import compute_regularization

matplotlib.use('QtAgg')


def cm2inch(cm):
    return cm / 2.54


def main(Lx=1000., dx=20., c_psi = 2.):
    # paragraph after eq. 10 of article
    #Lx = 100. * dx
    imax = int(Lx / dx + 1)
    step_left = 0.
    step_right = 1.0
    step_height = step_right - step_left

    refine = 100
    x_ana = np.zeros(refine*(imax-1) + 1, dtype=np.float64)
    ugiv_ana = np.zeros(refine*(imax-1) + 1, dtype=np.float64)
    utilde_ana = np.zeros(refine*(imax-1) + 1, dtype=np.float64)
    ubar_ana = np.zeros(refine*(imax-1) + 1, dtype=np.float64)
    x = np.zeros(imax, dtype=np.float64)
    ugiv = np.zeros(imax, dtype=np.float64)
    ubar = np.zeros(imax, dtype=np.float64)
    utilde = np.zeros(imax, dtype=np.float64)
    u = np.zeros(imax, dtype=np.float64)
    u0 = np.zeros(imax, dtype=np.float64)
    u1 = np.zeros(imax, dtype=np.float64)
    u0x = np.zeros(imax, dtype=np.float64)

    for i in range(0, refine*(imax-1)+1):
        x_ana[i] = i/(refine*(imax-1)) * Lx
    for i in range(0, imax):
        x[i] = i/(imax-1) * Lx


    # c_psis = [0.01, 0.0675, 0.0925, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1., 1.5, 2., 3., 4., 5., 7., 10.]
    # c_psis = [0.09, 0.1, 0.11]
    c_psis = np.linspace(0.01,10,1000)
    delta_utilde_ugiv = np.zeros(len(c_psis))
    delta_ubar_ugiv = np.zeros(len(c_psis))
    delta_ubar_utilde = np.zeros(len(c_psis))
    delta_i = -1
    for c_psi in c_psis:
        delta_i += 1
        Psi = c_psi * dx * dx

        # detect in which dx the step is
        xc = 1. / 2. * Lx
        jc = int(xc/dx)
        step_at = 0.0
        xc = x[jc] + step_at * dx

        for i in range(0, refine*(imax-1)+1):
            ugiv_ana[i] = 0.
            utilde_ana[i] = 0.5 * np.exp((x_ana[i] - xc)/np.sqrt(Psi))
            if x_ana[i] > xc:
                ugiv_ana[i] = 1.
                utilde_ana[i] = 1. - 0.5 * np.exp((xc - x_ana[i])/np.sqrt(Psi))
        for i in range(0, imax):
            ugiv[i] = 0.
            utilde[i] = 0.5 * np.exp((x[i] - xc)/np.sqrt(Psi))
            if x[i] > xc:
                ugiv[i] = 1.
                utilde[i] = 1. - 0.5 * np.exp((xc - x[i])/np.sqrt(Psi))
            u0[i] = ugiv[i]

        for i in range(0, imax):
            u[i] = 0.
            if x[i] >= xc:
                u[i] = 1.

        delta_c = (xc - x[jc])/dx
        alpha = c_psi
        beta  = 0.0
        if alpha != 1./8.:
            gamma = 0.5 /(alpha-1./8.)
            if alpha < 1./8.:
                beta = 1 + gamma + np.sqrt( (1. + gamma)**2 - 1)
            elif alpha > 1./8.:
                beta = 1 + gamma - np.sqrt( (1. + gamma)**2. - 1)

        u[jc] = 0.5 - delta_c * (1. - beta)/(1. + beta)
        for j in range(0, imax):
            if j < jc:
                u[j] = beta ** (jc-j) * u[jc]
                print('%d %.8f %.8f %.8f' % (j, x[j], beta, u[j]))
            elif j == jc:
                u[j] = 0.5 - delta_c * (1. - beta)/(1. + beta)
                print(' --- %d %.8f %.8f %.8f ---' % (j, x[j], beta, u[j]))
            else:
                u[j] = 1. - beta ** (j - jc) * (1.0 - u[jc])
                print('%d %.8f %.8f %.8f' % (j, x[j], beta, -u[j]+1.))
        for j in range(0, imax-1):
            for i in range(0, refine):
                k = j*refine + i
                fac = float(i)/(float(refine))
                ubar_ana[k] = (1.0 - fac) * u[j] + fac * u[j+1]
        ubar_ana[refine*(imax-1)] = u[imax-1]

        # delta_utilde[delta_i], delta_ubar[delta_i] = compute_regularization(c_psi, ugiv, u0, dx, imax)

        norm_utilde_ugiv = 0
        norm_ubar_ugiv = 0
        norm_ubar_utilde = 0
        for i in range(0, refine*(imax-1)+1):
            norm_utilde_ugiv += abs(utilde_ana[i] - ugiv_ana[i])
            norm_ubar_ugiv += abs(ubar_ana[i] - ugiv_ana[i])
            norm_ubar_utilde += abs(ubar_ana[i] - utilde_ana[i])
        delta_utilde_ugiv[delta_i] = norm_utilde_ugiv/refine
        delta_ubar_ugiv[delta_i] = norm_ubar_ugiv/refine
        delta_ubar_utilde[delta_i] = norm_ubar_utilde/refine
        # delta_ubar_utilde[delta_i] = abs(delta_ubar_ugiv[delta_i]) - abs(delta_utilde_ugiv[delta_i])
        a = 1

    # -------------------------------------------------------------------
    # plot the series
    fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(2, 1))

    ax2.set_ylabel('$L^1$-norm $\\longrightarrow$')  # after definition of spines
    ax2.set_xlabel('$c_{\psi} \\longrightarrow$')  # after definition of spines
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax2.set_xlim([0.01, 10.])
    ax2.set_ylim([0.01, 10.])

    ax2.grid(which='major', linestyle='-', linewidth='1', color = 'gray')
    ax2.grid(which='minor', linestyle='--', linewidth='0.5', color = 'red')

    tekst1 = '$||\widetilde{u} - u_{giv}||$'
    tekst2 = '$||\overline{u} - u_{giv}||$'
    tekst3 = '$||\overline{u} - \widetilde{u}||$'
    ax2.plot(c_psis, delta_utilde_ugiv, '-', color='b', label=tekst1)
    ax2.plot(c_psis, delta_ubar_ugiv, ':', color='b', label=tekst2)
    ax2.plot(c_psis, (delta_ubar_utilde), ':', color='r', label=tekst3)

    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, prop={"size": 12}, loc='upper left')

    # tekst = ('$c_{\Psi}$ = %.3f; $\Delta x$ = %.2f [m]; ${\Psi = c_{\Psi}\Delta x^2}$ = %.2f' % (c_psi, dx, Psi))
    # fig.text(0.125, 0.90, tekst, fontsize=10)
    tekst = ('Step offset from cell centre: %.2f $\Delta x$' % (step_at))
    fig.text(0.125, 0.94, tekst, fontsize=10)

    fig.set_size_inches(cm2inch(17.5), cm2inch(12.5))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()
    if not os.path.exists('figures'):
        os.mkdir('figures')
    if not os.path.exists('data'):
        os.mkdir('data')

    text = ('figures/regul_1d_step_at=%s_dx%s.pdf' % (
        str(step_at), str(dx)))
    plt.savefig(text, format='pdf')
    fig.canvas.manager.set_window_title('Regularization_1d_step_function')

    text = ('data/regul_1d_step_at=%s_dx%s_cpsi%s_to_%s.tek' % (
        str(step_at), str(dx), str(c_psis[0]), str(c_psis[-1])))
    with open(text, "w") as logfile:
        for j in range(0, len(c_psis)):
            logfile.write('*\n')
            logfile.write('* Step function\n')
            logfile.write('*\n')
            logfile.write('* column 1: x-coordinate\n')
            logfile.write('* column 2: u_giv_ana\n')
            logfile.write('* column 3: u_tilde_ana\n')
            logfile.write('* column 4: u_bar_ana\n')
            logfile.write('* column 5: abs(u_tilde_ana - u_giv_ana)\n')
            logfile.write('* column 6: abs(u_bar_ana - u_giv_ana)\n')
            logfile.write('* column 7: abs(u_bar_ana - u_tilde_ana)\n')
            logfile.write('*\n')
            text = ('* cpsi_%s\n' % str(c_psis[j]))
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
                abs(ubar_ana[i] - utilde_ana[i])))
                logfile.write(tekst)
        logfile.close()
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
