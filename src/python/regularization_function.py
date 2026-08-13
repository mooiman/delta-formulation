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

def get_ftilde(ugiv, ugiv_ana, c_psi, dx, nx, refine):
    # function_fitxxx.mw. See getftilde := proc( csm )
    a = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    b = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    c = np.ones(refine*(nx-1) + 1, dtype=np.float64)
    d = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    e = np.zeros(refine*(nx-1) + 1, dtype=np.float64)
    f = np.zeros(refine*(nx-1) + 1, dtype=np.float64)

    c_error = c_psi * dx * dx
    dx_small = dx / refine
    refine_bnd = int(refine/2)
    refine_bnd = refine

    for i in range(0, refine_bnd):
        a[i] = 0.
        b[i] = 0.0
        c[i] = 1.0  # 11. / 24.  # 1./12.  # 11. / 24.
        d[i] = 0.0
        e[i] = 0.0
        f[i] = ugiv[0]

    i = refine_bnd
    a[i] = 0.
    b[i] = 11. / 24.  # 1./12.  # 11. / 24.0.0
    c[i] = 14. / 24.  # 10./12. # 14. / 24.11. / 24.  # 1./12.  # 11. / 24.
    d[i] = -1. / 24.  # 1./12.  # -1. / 24.14. / 24.  # 10./12. # 14. / 24.
    e[i] = 0.0  # 1./12.  # -1. / 24.
    f[i] = ugiv[0]

    i = refine*(nx-1)+1 - 1 - refine_bnd
    a[i] = -1. / 24.  # 1./12.  # -1. / 24.
    b[i] = 14. / 24.  # 10./12. # 14. / 24.
    c[i] = 11. / 24.  # 1./12.  # 11. / 24.
    d[i] = 0.
    e[i] = 0.
    f[i] = ugiv[nx - 1]

    for i in range(refine*(nx-1)+1 - 1 - refine_bnd, refine*(nx-1)+1 - 1):
        a[i] = 0.
        b[i] = 0.0
        c[i] = 1.0  # 11. / 24.  # 1./12.  # 11. / 24.
        d[i] = 0.0
        e[i] = 0.0
        f[i] = ugiv[nx-1]

    # eq7
    for i in range(refine_bnd + 1, refine*(nx-1)+1 - 1 - refine_bnd):
        a[i] = 0.
        b[i] = 1. / 8. * dx_small - c_error / dx_small
        c[i] = 6. / 8. * dx_small + 2. * c_error / dx_small
        d[i] = 1. / 8. * dx_small - c_error / dx_small
        e[i] = 0.
        f[i] = (1./4. * ugiv_ana[i-1] + 3./4. * ugiv_ana[i]) * dx_small/2. + (3./4. * ugiv_ana[i] + 1./4. * ugiv_ana[i+1]) * dx_small/2.

    ftilde = thomas_algorithm_5(a, b, c, d, e, f)

    return ftilde

def compute_regularization(c_psi, ugiv, dx, nx, ugiv_ana, refine):
    iter_max = 500  # from funcfitxxx.mw
    csmscale = 1.0  # from funcfitxxx.mw

    diff_max0 = 0
    ugiv_xx_max = 0.
    rhside = np.zeros(nx, dtype=np.float64)

    ugiv_xx = np.zeros(nx, dtype=np.float64)
    ugiv_xixi = np.zeros(nx, dtype=np.float64)
    u0_xx = np.zeros(nx, dtype=np.float64)
    u0_xixi = np.zeros(nx, dtype=np.float64)
    fbar = np.zeros(nx, dtype=np.float64)
    ddfbar = np.zeros(nx, dtype=np.float64)
    a = np.zeros(nx, dtype=np.float64)
    b = np.ones(nx, dtype=np.float64)
    c = np.zeros(nx, dtype=np.float64)
    d = np.zeros(nx, dtype=np.float64)
    e = np.zeros(nx, dtype=np.float64)
    f = np.zeros(nx, dtype=np.float64)
    Err = np.zeros(nx, dtype=np.float64)
    csmold = np.zeros(nx, dtype=np.float64)
    csmnew = np.zeros(nx, dtype=np.float64)
    diff_abs = np.zeros(nx, dtype=np.float64)

    dx_small = dx / refine
    for i in range(1, nx - 1):
        csmnew[i] = dx * dx
    csmnew[0] = 0.0
    csmnew[nx - 1] = 0.0

    for i in range(1, nx - 1):
        if refine == 1:
            rhside[i] = dx * (1. / 8. * ugiv_ana[i - 1] + 6. / 8. * ugiv_ana[i] + 1. / 8. * ugiv_ana[i + 1])
        else:
            rhside[i] = 0.0
            for j in range(0, refine):
                k = i * refine + j - int(refine / 2) + 1
                # f[i] += 0.5 * (ugiv_ana[k] + ugiv_ana[k + 1]) * dx / refine
                rhside[i] += 0.5 * (ugiv_ana[k-1] + ugiv_ana[k]) * dx_small

    for it in range(0, iter_max):
        #=======================================================================
        # eq7 
        for i in range(1, nx - 1):
            csmleft = 0.5 * (csmnew[i - 1] + csmnew[i])
            csmright = 0.5 * (csmnew[i] + csmnew[i + 1])
            a[i] = 0.
            b[i] = 1./8. * dx - csmleft / dx
            c[i] = 6./8. * dx + (csmleft + csmright) / dx
            d[i] = 1./8. * dx - csmright / dx
            e[i] = 0.
            f[i] = rhside[i]

        i = 0
        a[i] = 0.
        b[i] = 0.0
        c[i] = 1./12.  # 11./24. #
        d[i] = 10./12. # 14./24. #
        e[i] = 1./12.  # -1./24. #
        f[i] = ugiv[i]

        i = nx - 1
        a[i] = 1./12.  #  -1./24. #
        b[i] = 10./12. #  14./24. #
        c[i] = 1./12.  #  11./24. #
        d[i] = 0.
        e[i] = 0.
        f[i] = ugiv[i]

        # solve eq7
        f = thomas_algorithm_5(a, b, c, d, e, f)
        for i in range(0, nx):
            fbar[i] = f[i]
        # ============================================================
        # Construct and solve system for smoothing coefficient
        # ============================================================
        for i in range(0, nx):
            csmold[i] = csmnew[i]

        # compute: dx^2 * u0_xixi
        for i in range(1, nx-1):
            ddfbar[i] = dx * dx * (fbar[i-1] - 2. * fbar[i]  + fbar[i+1])
        # Constant extrapolation at boundaries
        ddfbar[0] = ddfbar[1]
        ddfbar[nx-1] = ddfbar[nx-2]

        # eq 8
        c_error = c_psi
        for i in range(1, nx - 1):
            a[i] = 0.
            b[i] = 1./8. - c_error
            c[i] = 6./8. + 2. * c_error
            d[i] = 1./8. - c_error
            e[i] = 0.
            f[i] = np.abs(ddfbar[i])

        i = 0
        a[i] = 0.
        b[i] = 0.
        c[i] = 1.
        d[i] = -1.
        e[i] = 0.
        f[i] = 0.0  # == csmnew[0]
        i = nx - 1
        a[i] = 0.
        b[i] = -1.
        c[i] = 1.
        d[i] = 0.
        e[i] = 0.
        f[i] = 0.0  # == csmnew[nx - 1]

        # solve eq8
        f = thomas_algorithm_5(a, b, c, d, e, f)
        for i in range(0, nx):
            csmnew[i] = f[i]


        # ============================================================
        # Convergence error
        # ============================================================

        dcsmmax = (
            np.max(np.abs(csmnew - csmold))
            / np.max(csmold)
        )

        # ============================================================
        # Under-relaxation
        # ============================================================
        omega = 1.0
        for i in range(0, nx):
            csmnew[i] = omega * csmnew[i] + (1 - omega) * csmold[i]

        diff_max1 = 0.0
        for i in range(0, nx):
            diff_max1 = max(diff_max1, np.abs(csmnew[i] - csmold[i]))
        # print("csmnew - csmold:", diff_max1)
        diff_max10 = abs(diff_max1 - diff_max0)
        print("diff:", it, diff_max10)
        diff_max0 = diff_max1

        if diff_max10 < 1e-8:
            break

    return fbar, Err, csmnew, diff_max10, ugiv_xx_max
