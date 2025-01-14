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


def rtilde(d=0, kdx=0):
    tel = 1.
    noem = 1. + d * kdx * kdx
    return tel / noem


def rbar(d=0, kdx=0):
    if kdx == 0:
        tel = 0. + 0.j
    else:
        tel = 0. - 1.j * ( np.exp(0.5 * kdx) - np.exp(-0.5 * kdx))/kdx
    noem = 3./4. + 2. * d + (1./8. - d) * 2. * math.cos(kdx) + 0.j
    return  np.abs(tel / noem)


def main():
    imax = int(1000)
    x = np.linspace(0.0 , 1. , imax)
    # alpha = 1. / 8.
    errorFVEgridpoint = np.zeros(imax, dtype=np.float64)
    errorFVEcellcentre = np.zeros(imax, dtype=np.float64)
    errorFVEoptimal = np.zeros(imax, dtype=np.float64)
    for i in range(1, imax):
        kdx = x[i] * np.pi
        # ratio[i] = (1. - 2. * alpha) + 2. * d + (alpha - d) * 2. * math.cos(x[i])
        # errorFVEgridpoint[i] = np.abs(  rbar(0, kdx) - 1. )
        # errorFVEcellcentre[i] = np.abs( rbar(0, kdx) * np.cos(0.5 * kdx) - 1.)
        # errorFVEoptimal[i] = np.abs( rbar(0, kdx) * (1./3. + 2./3. * np.cos(0.5 * kdx)) - 1.)
        errorFVEgridpoint[i] = np.abs( 1. - (8. * np.sin(0.5 * kdx))/(kdx * (3. + np.cos(kdx))) )
        errorFVEcellcentre[i] = np.abs( 1. - (8. * np.sin(0.5 * kdx) * np.cos(0.5 * kdx))/(kdx * (3. + np.cos(kdx))) )
        errorFVEoptimal[i] = np.abs( 1. - (8. * np.sin(0.5 * kdx) * (1./3. + 2./3. * np.cos(0.5 * kdx)))/(kdx * (3. + np.cos(kdx))))

    # -------------------------------------------------------------------
    # plot the series
    fig, ax = plt.subplots(1)
    # plt.ylim(-2, 2)
    # plt.ylim(0, 4)
    # ax1[0].set_title('$u_t + u_x = 0$')
    ax.set_ylabel('Discretization error $\\longrightarrow$')  # after definition of spines
    ax.set_xlabel('k $\Delta x$  $\\longrightarrow$')  # after definition of spines
    ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
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

    tekst1 = 'errorFVEgridpoint'
    tekst2 = 'errorFVEcellcentre'
    tekst3 = 'errorFVEoptimal'
    ax.plot(x, errorFVEgridpoint, '-', color='b', label=tekst1)
    ax.plot(x, errorFVEcellcentre, '--', color='b', label=tekst2)
    # ax.plot(x, errorFVEoptimal, ':', color='b', label=tekst3)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={"size": 12}, loc='upper left')

    # tekst = ('ratio =  1/(3/4 + 2 d + (1/8 - d) 2 cos(kdx))')
    # fig.text(0.01, 0.91, tekst, fontsize=10)

    fig.set_size_inches(cm2inch(35.0), cm2inch(17.5))
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    if not os.path.exists('figures'):
        os.mkdir('figures')
    if not os.path.exists('data'):
        os.mkdir('data')

    text = ('figures/discr_error_u_u_giv.pdf')
    plt.savefig(text, format='pdf')
    plt.show()

    with open("data/disc_error_u_u_giv.tkl", "w") as logfile:
        logfile.write('* column 1: x-coordinate contains kdx\n')
        logfile.write('* column 2: errorFVEgridpoint\n')
        logfile.write('* column 3: errorFVEcellcentre \n')
        logfile.write('error\n')
        tekst = ('%.d %.d\n' % (imax, 83))
        logfile.write(tekst)
        for i in range(0, imax):
            tekst = ('%.8f %.8f %.8f\n' % (x[i], errorFVEgridpoint[i], errorFVEcellcentre[i]))
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
