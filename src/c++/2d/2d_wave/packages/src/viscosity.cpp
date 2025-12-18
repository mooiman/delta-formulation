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

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

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

    double cv_area_0 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_1 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_2 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_3 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // No contribution to the continuity equation
    //
    //--------------------------------------------------------------------------
    // q-momentum equation
    // 
    int col_sw = q_eq;
    int col_w  = q_eq + 3;
    int col_nw = q_eq + 6;
    int col_s  = q_eq + 9;
    int col_0  = q_eq + 12;
    int col_n  = q_eq + 15;
    int col_se = q_eq + 18;
    int col_e  = q_eq + 21;
    int col_ne = q_eq + 24;

    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta =  0.0;

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



    // scv_0 face_1
    n_xi =  0.0;
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
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_r * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./8. + dd_r * 1./4.) );
    add_value(values, col_sw    , theta * scvf_fac * (aa_q * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_sw + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_r * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
            + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
        ) * n_eta
        );

    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdeta_0;
    bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0;
    cc_r = visc_0 / scv_area_0 * 1.0;
    dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;


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
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. + dd_r * 3./4.) );
    add_value(values, col_s     , theta * scvf_fac * (aa_q * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./8. + dd_r * 1./4.) );
    add_value(values, col_se    , theta * scvf_fac * (aa_q * 1./8. - dd_r * 1./4.) );

    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_s  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_se + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_r * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
            + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
        ) * n_eta
        );

    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

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


    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

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


    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_0 * 1.0;
    dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

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
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_ne    , theta * scvf_fac * (aa_q * 1./8. + dd_r * 1./4.) );
    add_value(values, col_e     , theta * scvf_fac * (aa_q * 1./8. - dd_r * 1./4.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_ne + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_e  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_r * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
            + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
        ) * n_eta
        );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdeta_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdeta_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

    aa_r = visc_0 / scv_area_0 * rtheta_0 / (htheta_0 * htheta_0) * dhdxi_0;
    bb_r = visc_0 / scv_area_0 * -1.0 / htheta_0 * dhdxi_0;
    cc_r = visc_0 / scv_area_0 * 1.0;
    dd_r = visc_0 / scv_area_0 * - rtheta_0 / htheta_0;

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
    add_value(values, col_n     , theta * scvf_fac * (aa_q * 3./8. + dd_r * 3./4.) );
    add_value(values, col_0     , theta * scvf_fac * (aa_q * 3./8. - dd_r * 3./4.) );   
    add_value(values, col_nw    , theta * scvf_fac * (aa_q * 1./8. + dd_r * 1./4.) );
    add_value(values, col_w     , theta * scvf_fac * (aa_q * 1./8. - dd_r * 1./4.) );

    add_value(values, col_n  + 1, theta * scvf_fac * (bb_q * 3./8.  + cc_r * 3./4.) );
    add_value(values, col_0  + 1, theta * scvf_fac * (bb_q * 3./8.  - cc_r * 3./4.) );
    add_value(values, col_nw + 1, theta * scvf_fac * (bb_q * 1./8.  + cc_r * 1./4.) );
    add_value(values, col_w  + 1, theta * scvf_fac * (bb_q * 1./8.  - cc_r * 1./4.) );

    rhs[row + 1] += - (
        0.5 * (
            - dx_dxi_0 * dy_deta_0 * ( cc_r * drdxi_0  + dd_r * dhdxi_0)
            + dx_dxi_0 * dx_dxi_0  * ( cc_r * dqdeta_0 + dd_r * dhdeta_0)
        ) * n_eta
        );


    // scv_3 face_7
    n_xi = -1.0;
    n_eta =  0.0;

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

    aa_q = visc_0 / scv_area_0 * qtheta_0 /(htheta_0 * htheta_0) * dhdxi_0;
    bb_q = visc_0 / scv_area_0 * -1./ htheta_0 * dhdxi_0;
    cc_q = visc_0 / scv_area_0 * 1.0;
    dd_q = visc_0 / scv_area_0 * - qtheta_0 / htheta_0;

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

    //==========================================================================
    // r-momentum equation
    // 
    col_sw = r_eq;
    col_w  = r_eq + 3;
    col_nw = r_eq + 6;
    col_s  = r_eq + 9;
    col_0  = r_eq + 12;
    col_n  = r_eq + 15;
    col_se = r_eq + 18;
    col_e  = r_eq + 21;
    col_ne = r_eq + 24;
    {
    double visc_0;
    double htheta_0;
    double rtheta_0;
    double dhtheta_0;
    double drtheta_0;

    double aa_4;
    double bb_4;
    double cc_4;
    double dd_4;

    double dx  = 1.;
    double dxinv = 1.;
    double dy  = 1.;
    double dyinv = 1.;
    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dy * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_s , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_sw, theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_w , theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_sw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_s  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_sw + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_w  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_s , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_sw, theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_w , theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_sw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_s  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_sw + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_w  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_e], visc[p_se]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_s , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_se, theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_e , theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_se + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_s  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_se + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_e  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_se], visc[p_s]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dy * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_e , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_se, theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_s , theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);
    add_value(values, col_se + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_e  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_se + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_s  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dy * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_e , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_ne, theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);
    add_value(values, col_ne + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_e  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_ne + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_n  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_n], visc[p_ne], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_ne, theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_e , theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);
    add_value(values, col_ne + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_n  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_ne + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_e  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) );

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_n], visc[p_nw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_nw, theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_w , theta * scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);
    add_value(values, col_nw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_n  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) );
    add_value(values, col_nw + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) );
    add_value(values, col_w  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_3 face_7
    n_xi = -1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dyinv * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dyinv * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = 0.5 * dy * n_eta;
    add_value(values, col_0 , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_w , theta * scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_nw, theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.5 * 1. * dyinv) );
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);
    add_value(values, col_nw + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);

    add_value(values, col_0  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_w  + 2, theta * scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_nw + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );
    add_value(values, col_n  + 2, theta * scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.5 * 1. * dyinv) );

    rhs[row + 2] += - scvf_fac * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );
    }
//------------------------------------------------------------------------------

    return 0;
}
int viscosity_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& visc, size_t nx, size_t ny)                          // RHS vector [h, q, r]^{n}
{
    // Viscosity for post processing; WITHOUT integration over the control volumes, just the line-integral.
    double h;
    double q;
    double r;
    double nxi;
    double neta;
    double nxi_dl;
    double neta_dl;
    double dx = 1.;
    double dy = 1.;
    double dxinv = 1./dx;
    double dyinv = 1./dy;

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

            // scv_0 face_0
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            double visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
            double hn_0 = scvf_xi(hn[p_0], hn[p_w], hn[p_sw], hn[p_s]);
            double qn_0 = scvf_xi(qn[p_0], qn[p_w], qn[p_sw], qn[p_s]);
            double rn_0 = scvf_xi(rn[p_0], rn[p_w], rn[p_sw], rn[p_s]);

            double dhn_0 = dcdx_scvf_n(hn[p_0], hn[p_w], hn[p_s], hn[p_sw]);
            double dqn_0 = dcdx_scvf_n(qn[p_0], qn[p_w], qn[p_s], qn[p_sw]);
            double dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_t(hn[p_0], hn[p_s], hn[p_w], hn[p_sw]);
            double drn_0 = dcdy_scvf_t(rn[p_0], rn[p_s], rn[p_w], rn[p_sw]);
            double dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // scv_0 face_1
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
            hn_0 = scvf_eta(hn[p_0], hn[p_s], hn[p_sw], hn[p_w]);
            qn_0 = scvf_eta(qn[p_0], qn[p_s], qn[p_sw], qn[p_w]);
            rn_0 = scvf_eta(rn[p_0], rn[p_s], rn[p_sw], rn[p_w]);

            dhn_0 = dcdx_scvf_t(hn[p_0], hn[p_w], hn[p_s], hn[p_sw]);
            dqn_0 = dcdx_scvf_t(qn[p_0], qn[p_w], qn[p_s], qn[p_sw]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_n(hn[p_0], hn[p_s], hn[p_w], hn[p_sw]);
            drn_0 = dcdy_scvf_n(rn[p_0], rn[p_s], rn[p_w], rn[p_sw]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // sub control volume 1 ============================================
            // scv_1 face_2
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_e], visc[p_se]);
            hn_0 = scvf_eta(hn[p_0], hn[p_s], hn[p_e], hn[p_se]);
            qn_0 = scvf_eta(qn[p_0], qn[p_s], qn[p_e], qn[p_se]);
            rn_0 = scvf_eta(rn[p_0], rn[p_s], rn[p_e], rn[p_se]);

            dhn_0 = dcdx_scvf_t(hn[p_e], hn[p_0], hn[p_se], hn[p_s]);
            dqn_0 = dcdx_scvf_t(qn[p_e], qn[p_0], qn[p_se], qn[p_s]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_n(hn[p_0], hn[p_s], hn[p_e], hn[p_se]);
            drn_0 = dcdy_scvf_n(rn[p_0], rn[p_s], rn[p_e], rn[p_se]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // scv_1 face_3
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_se], visc[p_s]);
            hn_0 = scvf_xi(hn[p_0], hn[p_e], hn[p_se], hn[p_s]);
            qn_0 = scvf_xi(qn[p_0], qn[p_e], qn[p_se], qn[p_s]);
            rn_0 = scvf_xi(rn[p_0], rn[p_e], rn[p_se], rn[p_s]);

            dhn_0 = dcdx_scvf_n(hn[p_e], hn[p_0], hn[p_se], hn[p_s]);
            dqn_0 = dcdx_scvf_n(qn[p_e], qn[p_0], qn[p_se], qn[p_s]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_t(hn[p_0], hn[p_s], hn[p_e], hn[p_se]);
            drn_0 = dcdy_scvf_t(rn[p_0], rn[p_s], rn[p_e], rn[p_se]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // sub control volume 2 ============================================
            // scv_2 face_4
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
            hn_0 = scvf_xi(hn[p_0], hn[p_e], hn[p_ne], hn[p_n]);
            qn_0 = scvf_xi(qn[p_0], qn[p_e], qn[p_ne], qn[p_n]);
            rn_0 = scvf_xi(rn[p_0], rn[p_e], rn[p_ne], rn[p_n]);

            dhn_0 = dcdx_scvf_n(hn[p_e], hn[p_0], hn[p_ne], hn[p_n]);
            dqn_0 = dcdx_scvf_n(qn[p_e], qn[p_0], qn[p_ne], qn[p_n]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_t(hn[p_n], hn[p_0], hn[p_ne], hn[p_e]);
            drn_0 = dcdy_scvf_t(rn[p_n], rn[p_0], rn[p_ne], rn[p_e]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // scv_2 face_5
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            visc_0 = scvf_eta(visc[p_0], visc[p_n], visc[p_ne], visc[p_e]);
            hn_0 = scvf_eta(hn[p_0], hn[p_n], hn[p_ne], hn[p_e]);
            qn_0 = scvf_eta(qn[p_0], qn[p_n], qn[p_ne], qn[p_e]);
            rn_0 = scvf_eta(rn[p_0], rn[p_n], rn[p_ne], rn[p_e]);

            dhn_0 = dcdx_scvf_t(hn[p_e], hn[p_0], hn[p_ne], hn[p_n]);
            dqn_0 = dcdx_scvf_t(qn[p_e], qn[p_0], qn[p_ne], qn[p_n]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_n(hn[p_n], hn[p_0], hn[p_ne], hn[p_e]);
            drn_0 = dcdy_scvf_n(rn[p_n], rn[p_0], rn[p_ne], rn[p_e]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // sub control volume 3 ============================================
            // scv_3 face_6
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            visc_0 = scvf_eta(visc[p_0], visc[p_n], visc[p_nw], visc[p_w]);
            hn_0 = scvf_eta(hn[p_0], hn[p_n], hn[p_nw], hn[p_w]);
            qn_0 = scvf_eta(qn[p_0], qn[p_n], qn[p_nw], qn[p_w]);
                rn_0 = scvf_eta(rn[p_0], rn[p_n], rn[p_nw], rn[p_w]);

            dhn_0 = dcdx_scvf_t(hn[p_0], hn[p_w], hn[p_n], hn[p_nw]);
            dqn_0 = dcdx_scvf_t(qn[p_0], qn[p_w], qn[p_n], qn[p_nw]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;    
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_n(hn[p_n], hn[p_0], hn[p_nw], hn[p_w]);
            drn_0 = dcdy_scvf_n(rn[p_n], rn[p_0], rn[p_nw], rn[p_w]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);

            // scv_3 face_7
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
            hn_0 = scvf_xi(hn[p_0], hn[p_w], hn[p_nw], hn[p_n]);
            qn_0 = scvf_xi(qn[p_0], qn[p_w], qn[p_nw], qn[p_n]);
            rn_0 = scvf_xi(rn[p_0], rn[p_w], rn[p_nw], rn[p_n]);

            dhn_0 = dcdx_scvf_n(hn[p_0], hn[p_w], hn[p_n], hn[p_nw]);
            dqn_0 = dcdx_scvf_n(qn[p_0], qn[p_w], qn[p_n], qn[p_nw]);
            dd_1 = -2. * -visc_0 * qn_0/hn_0;
            rhs_q[p_0] += nxi_dl * (2. * -visc_0 * dxinv * dqn_0 - dd_1 * dxinv * dhn_0);

            dhn_0 = dcdy_scvf_t(hn[p_n], hn[p_0], hn[p_nw], hn[p_w]);
            drn_0 = dcdy_scvf_t(rn[p_n], rn[p_0], rn[p_nw], rn[p_w]);
            dd_4 = -2. * -visc_0 * rn_0/hn_0;
            rhs_r[p_0] += neta_dl * (2. * -visc_0 * dyinv * drn_0 - dd_4 * dyinv * dhn_0);
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


