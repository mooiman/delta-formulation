import numpy as np


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

def compute_regularization(c_psi, ugiv, dx, nx, ugiv_ana, refine):
    iter_max = 200
    diff_max0 = 0
    ugiv_xx_max = 0.

    u0 = np.zeros(nx, dtype=np.float64)
    ugiv_xx = np.zeros(nx, dtype=np.float64)
    ugiv_xixi = np.zeros(nx, dtype=np.float64)
    u0_xx = np.zeros(nx, dtype=np.float64)
    u0_xixi = np.zeros(nx, dtype=np.float64)
    a = np.zeros(nx, dtype=np.float64)
    b = np.ones(nx, dtype=np.float64)
    c = np.zeros(nx, dtype=np.float64)
    d = np.zeros(nx, dtype=np.float64)
    e = np.zeros(nx, dtype=np.float64)
    f = np.zeros(nx, dtype=np.float64)
    Err = np.zeros(nx, dtype=np.float64)
    psi = np.zeros(nx, dtype=np.float64)
    diff_abs = np.zeros(nx, dtype=np.float64)

    for i in range(0, nx):
        u0[i] = ugiv[i]

    for it in range(0, iter_max):
        u0_xx_max = 0.
        u0_xixi_max = 0.
        ugiv_xx_max = 0.
        ugiv_xixi_max = 0.
        for i in range(1, nx - 1):
            u0_xx[i] = (u0[i + 1] - 2. * u0[i] + u0[i - 1]) / (dx * dx)
            u0_xixi[i] = (u0[i + 1] - 2. * u0[i] + u0[i - 1])
            ugiv_xx[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1]) / (dx * dx)
            ugiv_xixi[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1])
            u0_xx_max = max(u0_xx_max, np.abs(u0_xx[i]))
            u0_xixi_max = max(u0_xixi_max, np.abs(u0_xixi[i]))
            ugiv_xx_max = max(ugiv_xx_max, np.abs(ugiv_xx[i]))
            ugiv_xixi_max = max(ugiv_xixi_max, np.abs(ugiv_xixi[i]))
        u0_xixi[0] = u0_xixi[1]
        u0_xixi[nx - 1] = u0_xixi[nx - 2]
        u0_xx[0] = u0_xx[1]
        u0_xx[nx - 1] = u0_xx[nx - 2]

        c_error = c_psi

        # eq8
        alpha = 0.125
        for i in range(1, nx - 1):
            a[i] = 0.
            b[i] = alpha - c_error
            c[i] = 1 - 2. * alpha + 2. * c_error
            d[i] = alpha - c_error
            e[i] = 0.
            f[i] = np.abs(u0_xixi[i])
        i = 0
        a[i] = 0.
        b[i] = 1.
        c[i] = -2.
        d[i] = 1.
        e[i] = 0.
        f[i] = 0.
        i = 1
        a[i] = 0.
        # b[i] = 1./12.
        # c[i] = 10./12.
        # d[i] = 1./12.
        b[i] = 11./24.
        c[i] = 14./24.
        d[i] = -1./24.
        e[i] = 0.
        f[i] = f[i]
        i = nx - 2
        a[i] = 0.
        b[i] = -1./24.
        c[i] = 14./24.
        d[i] = 11./24.
        e[i] = 0.
        f[i] = f[i]
        i = nx - 1
        a[i] = 1.
        b[i] = -2.
        c[i] = 1.
        d[i] = 0.
        e[i] = 0.
        f[i] = 0.
        
        # solve eq8
        Err = thomas_algorithm_5(a, b, c, d, e, f)

        #eq6
        for i in range(0, nx):
            psi[i] = c_psi * dx * dx * Err[i]

        # eq7 
        for i in range(1, nx - 1):
            psi_im12 = psi[i]  # 0.5 * (psi[i - 1] + psi[i])
            psi_ip12 = psi[i]  # 0.5 * (psi[i] + psi[i + 1])
            a[i] = 0.
            b[i] = alpha * dx - psi_im12 / dx
            c[i] = 0.5 * ((1. - 2. * alpha) * dx + (1. - 2. * alpha) * dx) + psi_im12 / dx + psi_ip12 / dx
            d[i] = alpha * dx - psi_ip12 / dx
            e[i] = 0.
            f[i] = 0.
            if refine == 1:
                f[i] = dx * (alpha * ugiv[i - 1] + (1. - 2. * alpha) * ugiv[i] + alpha * ugiv[i + 1])
            else:
                for j in range (0, refine):
                    k = i * refine + j - int(refine/2)+1
                    f[i] += 0.5 * (ugiv_ana[k] + ugiv_ana[k+1]) * dx / refine

        i = 0
        a[i] = 0.
        b[i] = 0.
        c[i] = 1.
        d[i] = -2.
        e[i] = 1.
        f[i] = 0.
        i = 1
        a[i] = 0.
        b[i] = 0.
        c[i] = 1.
        d[i] = 0.
        e[i] = 0.
        f[i] = ugiv[i]

        i = nx - 2
        a[i] = 0.
        b[i] = 0.
        c[i] = 1.
        d[i] = 0.
        e[i] = 0.
        f[i] = ugiv[i]
        i = nx - 1
        a[i] = 1.
        b[i] = -2.
        c[i] = 1.
        d[i] = 0.
        e[i] = 0.
        f[i] = 0.0
        # solve eq7
        u0 = thomas_algorithm_5(a, b, c, d, e, f)

        diff_max1 = 0.0
        sum_max1 = 0.0
        for i in range(0, nx):
            diff_max1 = max(diff_max1, np.abs(u0[i] - ugiv[i]))
            sum_max1 += np.abs(u0[i] - ugiv[i])
        # print("u0 - ugiv:", diff_max1)
        diff_max10 = abs(diff_max1 - diff_max0)
        print("diff:", it, diff_max10)
        diff_max0 = diff_max1

        if diff_max10 < 1e-8:
            break

    return u0, Err, psi, diff_max10, ugiv_xx_max

