import numpy as np

def thomas_algorithm(a, b, c, d):
    n = len(d)
    c_star = np.zeros(n, float)
    d_star = np.zeros(n, float)
    f = np.zeros(n, float)

    c_star[0] = c[0] / b[0]
    d_star[0] = d[0] / b[0]

    for i in range(1, n):
        m = 1.0 / (b[i] - a[i] * c_star[i - 1])
        c_star[i] = c[i] * m
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m

    f[n - 1] = d_star[n - 1]
    for i in range(n - 2, -1, -1):
        f[i] = d_star[i] - c_star[i] * f[i + 1]

    return f

def compute_regularization(c_psi, ugiv, dx, imax, ugiv_ana, refine):
    iter_max = 200
    diff_max0 = 0
    ugiv_xx_max = 0.

    u0 = np.zeros(imax, dtype=np.float64)
    ugiv_xx = np.zeros(imax, dtype=np.float64)
    ugiv_xixi = np.zeros(imax, dtype=np.float64)
    u0_xx = np.zeros(imax, dtype=np.float64)
    u0_xixi = np.zeros(imax, dtype=np.float64)
    a = np.zeros(imax, dtype=np.float64)
    b = np.ones(imax, dtype=np.float64)
    c = np.zeros(imax, dtype=np.float64)
    d = np.zeros(imax, dtype=np.float64)
    Err = np.zeros(imax, dtype=np.float64)
    psi = np.zeros(imax, dtype=np.float64)
    diff_abs = np.zeros(imax, dtype=np.float64)

    for i in range(0, imax):
        u0[i] = ugiv[i]

    for it in range(0, iter_max):
        u_xx_max = 0.
        u_xixi_max = 0.
        ugiv_xx_max = 0.
        ugiv_xixi_max = 0.
        for i in range(1, imax - 1):
            u0_xx[i] = (u0[i + 1] - 2. * u0[i] + u0[i - 1]) / (dx * dx)
            u0_xixi[i] = (u0[i + 1] - 2. * u0[i] + u0[i - 1])
            ugiv_xx[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1]) / (dx * dx)
            ugiv_xixi[i] = (ugiv[i + 1] - 2. * ugiv[i] + ugiv[i - 1])
            u_xx_max = max(u_xx_max, np.abs(u0_xx[i]))
            u_xixi_max = max(u_xixi_max, np.abs(u0_xixi[i]))
            ugiv_xx_max = max(ugiv_xx_max, np.abs(ugiv_xx[i]))
            ugiv_xixi_max = max(ugiv_xixi_max, np.abs(ugiv_xixi[i]))
        u0_xixi[0] = u0_xixi[1]
        u0_xixi[imax - 1] = u0_xixi[imax - 2]
        u0_xx[0] = u0_xx[1]
        u0_xx[imax - 1] = u0_xx[imax - 2]

        c_error = c_psi

        alpha = 0.125
        for i in range(1, imax - 1):
            a[i] = alpha - c_error
            b[i] = 1 - 2. * alpha + 2. * c_error
            c[i] = alpha - c_error
            d[i] = np.abs(u0_xixi[i])
        a[0] = 0.
        b[0] = 2.
        c[0] = -1.
        d[0] = d[1]
        a[imax - 1] = -1.
        b[imax - 1] = 2.
        c[imax - 1] = 0.
        d[imax - 1] = d[imax-2]
        Err = thomas_algorithm(a, b, c, d)

        for i in range(0, imax):
            psi[i] = c_psi * dx * dx * Err[i]

        for i in range(1, imax - 1):
            psi_1 = 0.5 * (psi[i - 1] + psi[i])
            psi_2 = 0.5 * (psi[i] + psi[i + 1])
            a[i] = alpha * dx - psi_1 / dx
            b[i] = 0.5 * ((1. - 2. * alpha) * dx + (1. - 2. * alpha) * dx) + psi_1 / dx + psi_2 / dx
            c[i] = alpha * dx - psi_2 / dx
            d[i] = 0.
            for j in range (0, refine):
                k = i * refine + j - int(refine/2)+1
                d[i] += 0.5 * (ugiv_ana[k] + ugiv_ana[k+1]) * dx / refine
            if refine == 1:
                d[i] = dx * (alpha * ugiv[i - 1] + (1. - 2. * alpha) * ugiv[i] + alpha * ugiv[i + 1])

        i = 0
        a[i] = 0.
        b[i] = 1.
        c[i] = 0.
        d[i] = ugiv[i]
        i = imax - 1
        a[i] = 0.
        b[i] = 1.
        c[i] = 0.
        d[i] = ugiv[i]
        u1 = thomas_algorithm(a, b, c, d)

        diff_max1 = 0.0
        sum_max1 = 0.0
        for i in range(0, imax):
            diff_abs[i] = np.abs(u1[i] - ugiv[i]) - 25.
            diff_max1 = max(diff_max1, np.abs(u1[i] - ugiv[i]))
            sum_max1 += np.abs(u1[i] - ugiv[i])
        # print("u1 - ugiv:", diff_max1)
        diff_max10 = abs(diff_max1 - diff_max0)
        print("diff:", it, diff_max10)
        diff_max0 = diff_max1
        u0 = u1

        if diff_max10 < 1e-8:
            break

    return u0, Err, psi, diff_max10, ugiv_xx_max

