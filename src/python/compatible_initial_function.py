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


def main(Lx=12., dx=2.):
    nx = int(Lx / dx) + 1 + 2  # two extra virtual points

    refine = 32
    x_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    u_ana = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    x = np.zeros(nx, dtype=np.float64)
    u0 = np.zeros(nx, dtype=np.float64)
    cv_ana = np.zeros(nx, dtype=np.float64)
    cv_coarse = np.zeros(nx, dtype=np.float64)
    cv_compatible = np.zeros(nx, dtype=np.float64)

    u1 = np.zeros(nx, dtype=np.float64)
    a = np.zeros(nx, dtype=np.float64)
    b = np.ones(nx, dtype=np.float64)
    c = np.zeros(nx, dtype=np.float64)
    d = np.zeros(nx, dtype=np.float64)
    e = np.zeros(nx, dtype=np.float64)
    f = np.zeros(nx, dtype=np.float64)

    for i in range(0, refine*(nx-1)+1):
        x_ana[i] = float(i) * dx/float(refine) - dx
    for i in range(0, nx):
        x[i] = float(i - 1) * dx

    mid = Lx/2.
    quart = Lx/4.
    a_coef = 1.
    for i in range(0, refine*(nx-1)+1):
        if x_ana[i] < mid:
            u_ana[i] = a_coef*(x_ana[i] - mid + quart) ** 2.
            #u_ana[i] = a_coef*(x_ana[i] - mid + quart) ** 4.
            #u_ana[i] = x_ana[i]
        else:
            u_ana[i] = -a_coef*(x_ana[i] - mid - quart) ** 2. + 2. * a_coef*(x_ana[refine] - mid + quart)**2.
            #u_ana[i] = -a_coef*(x_ana[i] - mid - quart) ** 4. + 2. * a_coef*(x_ana[0] - mid + quart)**4.
            #u_ana[i] = x_ana[i]
        #u_ana[i] = np.sin(1.* math.pi * x_ana[i]/Lx) +2.
    for i in range(0, nx):
        if x[i] < mid:
            u0[i] = a_coef*(x[i] - mid + quart) **2.
            #u0[i] = a_coef*(x[i] - mid + quart) **4.
            #u0[i] = x[i]
        else:
            u0[i] = -a_coef*(x[i] - mid - quart) ** 2. + 2. * a_coef*(x[1] - mid + quart) **2.
            #u0[i] = -a_coef*(x[i] - mid - quart) ** 4. + 2. * a_coef*(x[0] - mid + quart) **4.
            #u0[i] = x[i]
        #u0[i] = np.sin(1. * math.pi * x[i] / Lx) +2.

    # integral of the control volumes (xi-1/2, xi+1/2)
    for i in range(1, nx-1):
        cv_coarse[i] =  dx/2. * 0.25 * (3. * u0[i] +u0[i-1]) + dx/2. * 0.25 * (3. * u0[i] + u0[i+1])
    i = 0
    cv_coarse[i] = dx/2. * 0.25 * (3. * u0[i] + u0[i+1])
    i = nx-1
    cv_coarse[i] = dx/2. * 0.25 * (3. * u0[i] + u0[i-1])

    for i in range(1, nx-1):
        cv_ana[i] = 0.
        for j in range(0, refine):
            k = i * refine + j - int(refine/2)
            cv_ana[i] += dx/float(refine) * 0.5 * (u_ana[k] + u_ana[k+1])
    i = 0
    cv_ana[i] = 0.
    for j in range(0, int(refine/2)):
        k = i * refine + j
        cv_ana[i] += dx/float(refine) * 0.5 * (u_ana[k] + u_ana[k+1])
    i = nx-1
    cv_ana[i] = 0.
    for j in range(0, int(refine/2)):
        k = i * refine + j - int(refine/2)
        cv_ana[i] += dx/float(refine) * 0.5 * (u_ana[k] + u_ana[k+1])


    alpha = 0.125
    c_error = 0.0
    a[0] = 0.
    b[0] = 0.
    c[0] = 1./12.
    d[0] = 10./12.
    e[0] = 1./12.
    f[0] = u0[1]
    a[nx - 1] = 1./12.
    b[nx - 1] = 10./12.
    c[nx - 1] = 1./12.
    d[nx - 1] = 0.
    e[nx - 1] = 0.
    f[nx - 1] = u0[(nx-1)-1]
    for i in range(1, nx - 1):
        b[i] = dx * alpha - c_error
        c[i] = dx * (1 - 2. * alpha )+ 2. * c_error
        d[i] = dx * alpha - c_error
        f[i] = cv_ana[i]
    #u1 = thomas_algorithm_3(b, c, d, f)
    u1 = thomas_algorithm_5(a, b, c, d, e, f)

    for i in range(1, nx-1):
        cv_compatible[i] =  dx/2. * 0.25 * (3. * u1[i] +u1[i-1]) + dx/2. * 0.25 * (3. * u1[i] + u1[i+1])
    i = 0
    cv_compatible[i] = dx/2. * 0.25 * (3. * u1[i] + u1[i+1])
    i = nx-1
    cv_compatible[i] = dx/2. * 0.25 * (3. * u1[i] + u1[i-1])


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
        u_min = min(u_min, u1[i])
        u_max = max(u_max, u1[i])
    x_threshold = (x_max - x_min)* 0.025
    u_threshold = (u_max - u_min)* 0.05
    ax2.set_xlim([x_min - x_threshold, x_max + x_threshold])
    ax2.set_ylim([u_min - u_threshold, u_max + u_threshold])
    ax2.axvspan(x_min - x_threshold, x_min + 0.5*dx, color='gray', alpha=0.5)
    ax2.axvspan(x_max + x_threshold, x_max - 0.5*dx, color='gray', alpha=0.5)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax2.yaxis.set_major_formatter(y_formatter)

    tekst1 = 'Given function $u(x)$'
    tekst2 = '$\overline{u}(x)$'
    tekst3 = 'Compatible $\overline{u}(x)$'

    ax2.plot(x_ana, u_ana, '-', color='b', label=tekst1)
    ax2.plot(x, u0, '-.', color='g', marker = 'o', label=tekst2)
    ax2.plot(x, u1, '-', color='r', marker = 'o', label=tekst3)

    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='lower left')

    fig.set_size_inches(cm2inch(35.0), cm2inch(12.5))
    figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()
    if not os.path.exists('data'):
        os.mkdir('data')
    data_file =('data/compatible_initialization_lx=%s_dx=%s.tek' % (str(Lx), str(dx)))
    with open(data_file, "w") as logfile:
        logfile.write('*\n')
        logfile.write('* column 1: x-coordinate\n')
        logfile.write('* column 2: u0\n')
        logfile.write('* column 3: u1\n')
        logfile.write('* column 4: cv_analytic\n')
        logfile.write('* column 5: cv_coarse\n')
        logfile.write('* column 6: cv_compatible\n')
        logfile.write('data\n')
        tekst = ('%.d %.d\n' % (nx, 6))
        logfile.write(tekst)
        for i in range(0, nx):
            tekst = ('%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n' % (x[i], u0[i], u1[i], cv_ana[i], cv_coarse[i], cv_compatible[i]))
            logfile.write(tekst)
        logfile.close()

    if not os.path.exists('figures'):
        os.mkdir('figures')
    text = ('figures/compatible_initialization_lx=%s_dx=%s.pdf' % (str(Lx), str(dx)))
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
