#
# From M. Borsboom
# funcfit&smoothing_Fouriermodeanalyse.pdf
#
# Programmer: Jan Mooiman
#

import os  # file exists
import sys  # system
from datetime import datetime, date, time, timedelta
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

import numpy as np
import PyQt6

matplotlib.use('QtAgg')


def format_func(value, tick_number):
    return r"${0}$ $\pi$".format(value)


def cm2inch(cm):
    return cm / 2.54


def rtilde(c_psi=0, kdx=0):
    tel = 1.
    noem = 1. + c_psi * kdx * kdx
    return tel / noem


def rbar(c_psi=0, kdx=0):
    if kdx ==0:
        tel = 0.
    else:
        # leaving out sqrt(-1)
        tel = ( np.exp(0.5 * kdx) - np.exp(-0.5 * kdx))/kdx
    noem = 3./4. + 2. * c_psi + (1./8. - c_psi) * 2. * math.cos(kdx)
    return np.abs(tel / noem)


def main():
    imax = int(1000)
    x = np.linspace(0.0 , 1. , imax)
    alpha = 1. / 8.
    tilde2 = np.zeros(imax, dtype=np.float64)
    tilde4 = np.zeros(imax, dtype=np.float64)
    tilde8 = np.zeros(imax, dtype=np.float64)
    bar2 = np.zeros(imax, dtype=np.float64)
    bar4 = np.zeros(imax, dtype=np.float64)
    bar8 = np.zeros(imax, dtype=np.float64)
    d = 9999
    for i in range(0, imax):
        kdx = x[i] * np.pi
        # ratio[i] = (1. - 2. * alpha) + 2. * d + (alpha - d) * 2. * math.cos(x[i])
        tilde2[i] = rtilde(2, kdx)
        tilde4[i] = rtilde(4, kdx)
        tilde8[i] = rtilde(8, kdx)
        bar2[i] = rbar(2, kdx)
        bar4[i] = rbar(4, kdx)
        bar8[i] = rbar(8, kdx)


    # -------------------------------------------------------------------
    # plot the series
    fig, ax = plt.subplots(1)
    # plt.ylim(-2, 2)
    # plt.ylim(0, 4)
    # ax1[0].set_title('$u_t + u_x = 0$')
    ax.set_ylabel('$\widetilde{r} = | \widetilde{u}/u_{giv} |$ or $\overline{r} = | \overline{u}/u_{giv} |$ $\\longrightarrow$')  # after definition of spines
    ax.set_xlabel('k $\Delta x$ $\\longrightarrow$')  # after definition of spines

    major_grid_color = 'c'
#    plt.xaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
#    plt.yaxis.grid(True, 'major', linestyle='-', linewidth=0.5, color=major_grid_color)
    ax.set_xlim([1.e-2, 1.0])
    ax.set_ylim([1.e-3, 2.])
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    ax.grid(which='major', linestyle='-', linewidth='1', color = 'gray')
    ax.grid(which='minor', linestyle='--', linewidth='0.5', color = 'red')

    tekst1 = '$\widetilde{r}$ ($c_{\Psi}=2$)'
    tekst2 = '$\widetilde{r}$ ($c_{\Psi}=4$)'
    tekst3 = '$\widetilde{r}$ ($c_{\Psi}=8$)'
    tekst4 = '$\overline{r}$ ($c_{\Psi}=2$)'
    tekst5 = '$\overline{r}$ ($c_{\Psi}=4$)'
    tekst6 = '$\overline{r}$ ($c_{\Psi}=8$)'
    ax.plot(x, tilde2, '-', color='r', label=tekst1)
    ax.plot(x, tilde4, '--', color='r', label=tekst2)
    ax.plot(x, tilde8, ':', color='r', label=tekst3)
    ax.plot(x, bar2, '-', color='b', label=tekst4)
    ax.plot(x, bar4, '--', color='b', label=tekst5)
    ax.plot(x, bar8, ':', color='b', label=tekst6)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={"size": 12}, loc='lower left')

    # tekst = ('ratio =  1/(3/4 + 2 d + (1/8 - d) 2 cos(kdx))')
    # fig.text(0.01, 0.91, tekst, fontsize=10)
    #tekst = ('$d = \Psi/\Delta x^2$')
    #fig.text(0.125, 0.89, tekst, fontsize=10)

    fig.set_size_inches(cm2inch(35.0), cm2inch(17.5))
    figManager = plt.get_current_fig_manager()
    ##figManager.window.showMaximized()
    if not os.path.exists('figures'):
        os.mkdir('figures')
    if not os.path.exists('data'):
        os.mkdir('data')

    text = ('figures/ratio_u_u_giv.pdf')
    plt.savefig(text, format='pdf')
    plt.show()

    with open("data/ratio_u_u_giv.tkl", "w") as logfile:
        logfile.write('* Ratio abs(u_giv/u)\n')
        logfile.write('* column 1: x-coordinate contains kdx\n')
        logfile.write('* column 2: ratio rtilde d = Psi/dx^2 = 2\n')
        logfile.write('* column 3: ratio rbar d = Psi/dx^2 = 2\n')
        logfile.write('ratio\n')
        tekst = ('%.d %.d\n' % (imax, 2))
        logfile.write(tekst)
        for i in range(0, imax):
            tekst = ('%.8f %.8f %.8f\n' % (x[i], tilde2[i], bar2[i]))
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
