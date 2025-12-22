//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
//    Copyright (C) 2025 Jan Mooiman
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------------

#include "viscosity.h"
#include "interpolations.h"
#include "jacobians.h"

//   nw - - - - - - - n - - - - - - - ne        nw - - - - - - - n - - - - - - - ne
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    | - - - - - x - | - x - - - - - |          | - - - - - 6 - | - 5 - - - - - |
//    |       |       |       |       |          |       |       |       |       |
//    |       x       |       x       |          |       7 scv_3 | scv_2 4       |
//    |       |       |       |       |          |       |       |       |       |
//    w - - - - - - - 0 - - - - - - - e          w - - - - - - - 0 - - - - - - - e
//    |       |       |       |       |          |       |       |       |       |
//    |       x       |       x       |          |       0 scv_0 | scv_1 3       |
//    |       |       |       |       |          |       |       |       |       |
//    | - - - - - x - | - x - - - - - |          | - - - - - 1 - | - 2 - - - - - |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//   sw - - - - - - - s - - - - - - - se        sw - - - - - - - s - - - - - - - se

int viscosity_matrix_and_rhs(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& visc, double theta, size_t nx, size_t ny)
{
    double n_xi = 0.0;  // xi-component of the outward normal vector
    double n_eta = 0.0;  // eta-component of the outward normal vector
    double scvf_fac;

    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    int p_0 = c_eq/(3*27);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (p_0 % ny == 0) { return 1; }  // south boundary
    if ((p_0 + 1) % ny == 0) { return 2; }  // north boundary

    size_t p_sw = p_0 - ny - 1;
    size_t p_w  = p_0 - ny;
    size_t p_nw = p_0 - ny + 1;
    size_t p_s  = p_0 - 1; 
    size_t p_n  = p_0 + 1;
    size_t p_se = p_0 + ny - 1;
    size_t p_e  = p_0 + ny;
    size_t p_ne = p_0 + ny + 1;

    std::vector<double> x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
    std::vector<double> y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
    double scv_area_0 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
    y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
    double scv_area_1 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
    y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
    double scv_area_2 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
    y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
    double scv_area_3 = polygon_area(x_pol, y_pol);

    //double cv_area_0 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    //double cv_area_1 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    //double cv_area_2 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    //double cv_area_3 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;

    //==========================================================================
    // c-equation
    // 
    // No contribution to the continuity equation
    //
    //==========================================================================
    // q-momentum equation
    //==========================================================================
    bool q_mom = true;
    bool r_mom = true;
    if (q_mom)
    {
    size_t col_sw = q_eq;
    size_t col_w  = q_eq + 3;
    size_t col_nw = q_eq + 6;
    size_t col_s  = q_eq + 9;
    size_t col_0  = q_eq + 12;
    size_t col_n  = q_eq + 15;
    size_t col_se = q_eq + 18;
    size_t col_e  = q_eq + 21;
    size_t col_ne = q_eq + 24;

    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta = 0.0;

    double visc_0   = scvf_xi(visc[p_0], visc[p_w], visc[p_s], visc[p_sw]);
    double htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    double qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
    double rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
    double dhdxi_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    double dqdxi_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
    double drdxi_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
    double dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    double dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
    double drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    double dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
    double dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

    double aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    double bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdxi_0;
    double cc_q = visc_0 / scv_area_0 * 1.0;
    double dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_sw + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
          -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
        ) * n_xi
        );

    //
    // sub control volume 0 ============================================
    // scv_0 face_1
    n_xi = 0.0;
    n_eta = 0.0;  -1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
    dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
    dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
    dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    double aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    double bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
    double cc_r = visc_0 / scv_area_0 * 1.0;
    double dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dx_dxi_0 * dy_deta_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_sw + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );

    scvf_fac = 0.5 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_sw + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            + dx_dxi_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_eta
        );
    rhs[row + 2] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_eta
        );

    //
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi = 0.0;
    n_eta = 0.0;  -1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_e]);
    dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_e]);
    dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
    dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);

    aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_1 * 1.0;
    dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_1 * 1.0;
    dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dx_dxi_0 * dy_deta_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );   
    add_value(values, col_e     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_se    , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );
    add_value(values, col_e  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_se + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );

    scvf_fac = 0.5 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_se    , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_se + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            + dx_dxi_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_eta
        );
    rhs[row + 2] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_eta
        );

    //
    // sub control volume 1 ============================================
    // scv_1 face_3
    n_xi = 1.0;
    n_eta = 0.0;

    visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_se], visc[p_s]);
    htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
    dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
    dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
    dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

    aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_1 * 1.0;
    dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_se    , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_se + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
          -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
        ) * n_xi
        );

    //
    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi = 1.0;
    n_eta = 0.0;

    visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
    dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);

    aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_2 * 1.0;
    dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_ne + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
          -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
        ) * n_xi
        );

    //
    // sub control volume 2 ============================================
    // scv_2 face_5
    n_xi = 0.0;
    n_eta = 0.0;  1.0;

    visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_ne], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
    dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);

    aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_2 * 1.0;
    dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_2 * 1.0;
    dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dx_dxi_0 * dy_deta_0 * n_eta;
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_e     , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );

    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );
    add_value(values, col_ne + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_e  + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );

    scvf_fac = 0.5 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_ne + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            + dx_dxi_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_eta
        );
    rhs[row + 2] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_eta
        );

    //
    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi = 0.0;
    n_eta = 0.0;  1.0;

    visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_nw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
    dhdeta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    dqdeta_0 = dcdy_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    drdeta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
    dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);

    aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_3 * 1.0;
    dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_3 * 1.0;
    dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dx_dxi_0 * dy_deta_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );   
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 1./2. + dd_r * 1./2.) );
    add_value(values, col_nw    , theta * scvf_fac * (aa_r * 1./2. - dd_r * 1./2.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );
    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 1./2.  + cc_r * 1./2.) );
    add_value(values, col_nw + 2, theta * scvf_fac * (bb_r * 1./2.  - cc_r * 1./2.) );

    scvf_fac = 0.5 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_nw    , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_nw + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            + dx_dxi_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_eta
        );
    rhs[row + 2] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_eta
        );

    //
    // sub control volume 3 ============================================
    // scv_3 face_7
    n_xi = -1.0;
    n_eta = 0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_n], visc[p_nw]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
    dhdxi_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    dqdxi_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
    drdxi_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);
    dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);

    aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_3 * 1.0;
    dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_q * 3./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 3./8. - dd_q * 3./4.) );   
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 1./8. + dd_q * 1./4.) );
    add_value(values, col_nw    , theta * scvf_fac * (aa_q * 1./8. - dd_q * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_q * 3./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_q * 3./4.) );
    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_q * 1./4.) );
    add_value(values, col_nw + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_q * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
          -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
        ) * n_xi
        );
    }
    //==========================================================================
    // r-momentum equation
    //==========================================================================
    // 
    if (r_mom)
    {
    size_t col_sw = r_eq;
    size_t col_w  = r_eq + 3;
    size_t col_nw = r_eq + 6;
    size_t col_s  = r_eq + 9;
    size_t col_0  = r_eq + 12;
    size_t col_n  = r_eq + 15;
    size_t col_se = r_eq + 18;
    size_t col_e  = r_eq + 21;
    size_t col_ne = r_eq + 24;

    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    double n_xi = -1.0;
    double n_eta = 0.0;

    double visc_0   = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    double htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    double qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    double rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    double dhdxi_0  = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    double dqdxi_0  = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
    double drdxi_0  = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
    double dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    double dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
    double drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    double dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
    double dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

    double aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    double bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
    double cc_q = visc_0 / scv_area_0 * 1.0;
    double dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    double aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    double bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
    double cc_r = visc_0 / scv_area_0 * 1.0;
    double dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_sw + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    scvf_fac = 0.5 * dy_deta_0 * dx_dxi_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_sw + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );

    rhs[row + 1] += - (
        0.5 * (
            + dy_deta_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_xi
        );
    rhs[row + 2] += - (
        0.5 * (
            - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_xi
        );

    //
    // sub control volume 0 ============================================
    // scv_0 face_1
    n_xi = 0.0;
    n_eta = -1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
    dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
    dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
    dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);

    aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
    bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdeta_0;
    cc_r = visc_0 / scv_area_0 * 1.0;
    dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_sw + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    rhs[row + 2] += - (
        0.5 * (
          -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
        ) * n_eta
        );

    //
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi = 0.0;
    n_eta = -1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_s]);
    dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_s]);
    dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
    dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);

    aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
    bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdeta_0;
    cc_r = visc_0 / scv_area_1 * 1.0;
    dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_sw + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    rhs[row + 2] += - (
        0.5 * (
          -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
        ) * n_eta
        );

    //
    // sub control volume 1 ============================================
    // scv_1 face_3
    n_xi = 1.0;
    n_eta = 0.0;

    visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_se], visc[p_s]);
    htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
    drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
    dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
    dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);
    dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);

    aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_1 * 1.0;
    dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_1 * 1.0;
    dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_e     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_se    , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_e  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_se + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_s  + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    scvf_fac = 0.5 * dy_deta_0 * dx_dxi_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );   
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_se    , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_se + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );

    rhs[row + 1] += - (
        0.5 * (
            + dy_deta_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_xi
        );
    rhs[row + 2] += - (
        0.5 * (
            - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_xi
        );
    
    //
    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi = 1.0;
    n_eta = 0.0;

    visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
    dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);

    aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_2 * 1.0;
    dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_2 * 1.0;
    dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_e     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_e  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_ne + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    scvf_fac = 0.5 * dy_deta_0 * dx_dxi_0 * n_xi;
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );
    add_value(values, col_ne + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );

    rhs[row + 1] += - (
        0.5 * (
            + dy_deta_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_xi
        );
    rhs[row + 2] += - (
        0.5 * (
            - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_xi
        );

    //
    // sub control volume 2 ============================================
    // scv_2 face_5
    n_xi = 0.0;
    n_eta = 1.0;

    visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_ne], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
    dhdeta_0 = dcdx_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    dqdeta_0 = dcdx_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
    drdeta_0 = dcdx_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
    dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);

    aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
    bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdeta_0;
    cc_r = visc_0 / scv_area_2 * 1.0;
    dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_e     , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_ne + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_e  + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    rhs[row + 2] += - (
        0.5 * (
          -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
        ) * n_eta
        );

    //
    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi = 0.0;
    n_eta = 1.0;

    visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_nw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);
    dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
    drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
    dhdeta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    dqdeta_0 = dcdy_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    drdeta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
    dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);

    aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
    bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdeta_0;
    cc_r = visc_0 / scv_area_3 * 1.0;
    dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

    scvf_fac = - 0.5 * 2.0 * dx_dxi_0 * dx_dxi_0 * n_eta;
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_nw    , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_nw + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    rhs[row + 2] += - (
        0.5 * (
          -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
        ) * n_eta
        );
    
    //
    // sub control volume 3 ============================================
    // scv_3 face_7
    n_xi = -1.0;
    n_eta = 0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    dhdxi_0  = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    dqdxi_0  = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
    drdxi_0  = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
    dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
    drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);
    dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);

    aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_3 * 1.0;
    dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_3 * 1.0;
    dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

    scvf_fac = -0.5 * dy_deta_0 * dy_deta_0 * n_xi;
    add_value(values, col_0     , theta * scvf_fac * (aa_r * 3./8. + dd_r * 3./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_r * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_n     , theta * scvf_fac * (aa_r * 1./8. + dd_r * 1./4.) );
    add_value(values, col_nw    , theta * scvf_fac * (aa_r * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 2, theta * scvf_fac * (bb_r * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_w  + 2, theta * scvf_fac * (bb_r * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_n  + 2, theta * scvf_fac * (bb_r * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_nw + 2, theta * scvf_fac * (bb_r * 1./8.  - cc_r * 1./4.) );

    scvf_fac = 0.5 * dy_deta_0 * dx_dxi_0 * n_xi;
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );   
    add_value(values, col_nw    , theta * scvf_fac * (aa_q * 1./2. + dd_q * 1./2.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./2. - dd_q * 1./2.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );
    add_value(values, col_nw + 1, theta * scvf_fac * (bb_q * 1./2.  + cc_q * 1./2.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./2.  - cc_q * 1./2.) );

    rhs[row + 1] += - (
        0.5 * (
            + dy_deta_0 * dx_dxi_0  * ( cc_q * dqdeta_0 + dd_q * dhdeta_0)
        ) * n_xi
        );
    rhs[row + 2] += - (
        0.5 * (
            - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
        ) * n_xi
        );
    }
//------------------------------------------------------------------------------

    return 0;
}
int viscosity_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& visc, size_t nx, size_t ny)                          // RHS vector [h, q, r]^{n}
{
    // Viscosity for post processing; WITHOUT integration over the control volumes, just the line-integral.

    std::fill_n(rhs_q.data(), rhs_q.size(), 0.0);
    std::fill_n(rhs_r.data(), rhs_r.size(), 0.0);

    for (size_t i = 1; i < nx - 1; ++i)
    {
        for (size_t j = 1; j < ny - 1; ++j)
        {
            size_t p_0  = viscosity_idx(i    , j    , ny); // central point of control volume
            size_t p_sw = viscosity_idx(i - 1, j - 1, ny);  
            size_t p_s  = viscosity_idx(i    , j - 1, ny);  
            size_t p_se = viscosity_idx(i + 1, j - 1, ny);  
            size_t p_w  = viscosity_idx(i - 1, j    , ny);  
            size_t p_e  = viscosity_idx(i + 1, j    , ny);  
            size_t p_nw = viscosity_idx(i - 1, j + 1, ny);  
            size_t p_n  = viscosity_idx(i    , j + 1, ny);  
            size_t p_ne = viscosity_idx(i + 1, j + 1, ny);  

            std::vector<double> x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
            std::vector<double> y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
            double scv_area_0 = polygon_area(x_pol, y_pol);

            x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
            y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
            double scv_area_1 = polygon_area(x_pol, y_pol);

            x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
            y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
            double scv_area_2 = polygon_area(x_pol, y_pol);

            x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
            y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
            double scv_area_3 = polygon_area(x_pol, y_pol);

            //==================================================================
            // q-momentum equation
            //==================================================================
            // scv_0 face_0
            double n_xi = -1.0;
            double n_eta = 0.0;

            double visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
            double htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
            double qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
            double rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
            double dhdxi_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
            double dqdxi_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
            double drdxi_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
            double dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
            double dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
            double drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

            double dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
            double dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
            double dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
            double dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

            double aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
            double bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdxi_0;
            double cc_q = visc_0 / scv_area_0 * 1.0;
            double dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
                ) * n_xi;

            // scv_0 face_1
            n_xi = 0.0;
            n_eta = -1.0;

            visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
            htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
            qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
            rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

            dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
            dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
            dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
            dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);

            aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_0 * 1.0;
            dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

            double aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            double bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
            double cc_r = visc_0 / scv_area_0 * 1.0;
            double dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                    - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_eta;

            // scv_1 face_2
            n_xi = 0.0;
            n_eta = -1.0;

            visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
            htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
            qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
            rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

            dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_e]);
            dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_e]);
            dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
            dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);

            aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_1 * 1.0;
            dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_1 * 1.0;
            dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

            rhs_q[p_0] += 
            0.5 * (
                - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
            ) * n_eta;

            // scv_1 face_3
            n_xi = 1.0;
            n_eta = 0.0;

            visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_se], visc[p_s]);
            htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
            qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
            rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

            dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
            dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
            dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
            dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

            aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
            bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdxi_0;
            cc_q = visc_0 / scv_area_1 * 1.0;
            dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
                ) * n_xi;

            // scv_2 face_4
            n_xi = 1.0;
            n_eta = 0.0;

            visc_0   = scvf_xi(visc[p_e], visc[p_0], visc[p_ne], visc[p_n]);
            htheta_0 = scvf_xi(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
            qtheta_0 = scvf_xi(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
            rtheta_0 = scvf_xi(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

            dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
            dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
            dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
            dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);

            aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
            bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdxi_0;
            cc_q = visc_0 / scv_area_2 * 1.0;
            dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
                ) * n_xi;

            // scv_2 face_5
            n_xi = 0.0;
            n_eta = 1.0;

            visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_ne], visc[p_e]);
            htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

            dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
            dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
            dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
            dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);

            aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_2 * 1.0;
            dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_2 * 1.0;
            dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                    - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_eta;

            // scv_3 face_6
            n_xi = 0.0;
            n_eta = 1.0;

            visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_nw], visc[p_w]);
            htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
            dhdeta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            dqdeta_0 = dcdy_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            drdeta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

            dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
            dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
            dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
            dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);

            aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_3 * 1.0;
            dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_3 * 1.0;
            dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                    - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_eta;

            // scv_3 face_7
            n_xi = -1.0;
            n_eta = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_n], visc[p_nw]);
            htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
            qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
            rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
            dhdxi_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
            dqdxi_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
            drdxi_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

            dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
            dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
            dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);
            dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);

            aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
            bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdxi_0;
            cc_q = visc_0 / scv_area_3 * 1.0;
            dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

            rhs_q[p_0] += 
                0.5 * (
                -2. * dy_deta_0 * dy_deta_0 * ( cc_q * dqdxi_0  + dd_q * dhdxi_0 )
                ) * n_xi;
            //==================================================================
            // r-momentum equation
            //==================================================================
            // scv_0 face_0
            n_xi = -1.0;
            n_eta = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
            htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
            qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
            rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

            dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_w], x[p_sw]);
            dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_w], y[p_sw]);
            dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
            dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);

            aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_0 * 1.0;
            dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_0 * 1.0;
            dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                    - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dy_deta_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_xi;

            // scv_0 face_1
            n_xi = 0.0;
            n_eta = -1.0;

            visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
            htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
            qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
            rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_s], rtheta[p_sw]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_w], qtheta[p_sw]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

            dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
            dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
            dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
            dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);

            aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
            bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdeta_0;
            cc_r = visc_0 / scv_area_0 * 1.0;
            dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                  -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
                ) * n_eta;

            // scv_1 face_2
            n_xi = 0.0;
            n_eta = -1.0;

            visc_0   = scvf_eta(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
            htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
            qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
            rtheta_0 = scvf_eta(qtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

            dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_e]);
            dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_e]);
            dx_deta_0 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
            dy_deta_0 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);

            aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
            bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdeta_0;
            cc_r = visc_0 / scv_area_1 * 1.0;
            dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                  -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
                ) * n_eta;

            // scv_1 face_3
            n_xi = 0.0;
            n_eta = -1.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
            htheta_0 = scvf_xi(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
            qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
            rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_se], rtheta[p_s]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_0], qtheta[p_s], qtheta[p_e], qtheta[p_se]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

            dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
            dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
            dx_deta_0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);
            dy_deta_0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);

            aa_q = visc_0 / scv_area_1 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_1 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_1 * 1.0;
            dd_q = visc_0 / scv_area_1 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_1 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_1 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_1 * 1.0;
            dd_r = visc_0 / scv_area_1 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                    - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dy_deta_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_xi;

            // scv_2 face_4
            n_xi = 1.0;
            n_eta = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
            htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
            qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
            rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

            dx_dxi_0  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
            dy_dxi_0  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
            dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
            dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);

            aa_q = visc_0 / scv_area_2 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_2 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_2 * 1.0;
            dd_q = visc_0 / scv_area_2 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_2 * 1.0;
            dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                    - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dy_deta_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_xi;

            // scv_2 face_5
            n_xi = 0.0;
            n_eta = 1.0;

            visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_ne], visc[p_e]);
            htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_e], rtheta[p_0], rtheta[p_ne], rtheta[p_n]);
            dhdeta_0 = dcdx_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
            dqdeta_0 = dcdx_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_ne], qtheta[p_e]);
            drdeta_0 = dcdx_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

            dx_dxi_0  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
            dy_dxi_0  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
            dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
            dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);

            aa_r = visc_0 / scv_area_2 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
            bb_r = visc_0 / scv_area_2 * -1.0 / htheta_0 * dhdeta_0;
            cc_r = visc_0 / scv_area_2 * 1.0;
            dd_r = visc_0 / scv_area_2 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                  -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
                ) * n_eta;

            // scv_3 face_6
            n_xi = 0.0;
            n_eta = 1.0;

            visc_0   = scvf_eta(visc[p_n], visc[p_0], visc[p_nw], visc[p_w]);
            htheta_0 = scvf_eta(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            qtheta_0 = scvf_eta(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            rtheta_0 = scvf_eta(qtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);
            dhdxi_0  = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
            dqdxi_0  = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
            drdxi_0  = dcdx_scvf_t(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
            dhdeta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            dqdeta_0 = dcdy_scvf_n(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            drdeta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

            dx_dxi_0  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
            dy_dxi_0  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
            dx_deta_0 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
            dy_deta_0 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);

            aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
            bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdeta_0;
            cc_r = visc_0 / scv_area_3 * 1.0;
            dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                  -2. * dx_dxi_0 * dx_dxi_0 * ( cc_r * drdeta_0  + dd_r * dhdeta_0 )
                ) * n_eta;

            // scv_3 face_7
            n_xi = -1.0;
            n_eta = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
            htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
            qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
            rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
            dhdxi_0  = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
            dqdxi_0  = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);
            drdxi_0  = dcdx_scvf_n(rtheta[p_0], rtheta[p_w], rtheta[p_n], rtheta[p_nw]);
            dhdeta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
            dqdeta_0 = dcdy_scvf_t(qtheta[p_n], qtheta[p_0], qtheta[p_nw], qtheta[p_w]);
            drdeta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

            dx_dxi_0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
            dy_dxi_0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
            dx_deta_0 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);
            dy_deta_0 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);

            aa_q = visc_0 / scv_area_3 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
            bb_q = visc_0 / scv_area_3 * -1./ htheta_0 * dhdeta_0;
            cc_q = visc_0 / scv_area_3 * 1.0;
            dd_q = visc_0 / scv_area_3 * - qtheta_0 / htheta_0;

            aa_r = visc_0 / scv_area_3 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
            bb_r = visc_0 / scv_area_3 * -1.0 / htheta_0 * dhdxi_0;
            cc_r = visc_0 / scv_area_3 * 1.0;
            dd_r = visc_0 / scv_area_3 * - rtheta_0 / htheta_0;

            rhs_r[p_0] += 
                0.5 * (
                    - dy_deta_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
                    + dy_deta_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
                ) * n_xi;
        }
    }
    return 0;
}
inline size_t viscosity_idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}

inline void add_value(double * values, size_t col, double data){ 
    values[col] += data; 
}


