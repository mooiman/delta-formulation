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

#include "convection.h"
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

int convection_matrix_and_rhs(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double theta, size_t nx, size_t ny)
{
    double h;
    double q;
    double r;
    double nxi = 0.0;  // xi-component of the outward normal vector
    double neta = 0.0;  // eta-component of the outward normal vector
    double scvf_fac;

    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    int p_0 = c_eq/(3*27);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (std::fmod(p_0, ny) == 0) { return 1; }  // south boundary
    if (std::fmod(p_0 + 1, ny) == 0) { return 2; }  // north boundary

    size_t p_sw = p_0 - ny - 1;
    size_t p_w  = p_0 - ny;
    size_t p_nw = p_0 - ny + 1;
    size_t p_s  = p_0 - 1; 
    size_t p_n  = p_0 + 1;
    size_t p_se = p_0 + ny - 1;
    size_t p_e  = p_0 + ny;
    size_t p_ne = p_0 + ny + 1;

    // at face 0
    double dy_dxi_f0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_dxi_f0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_deta_f0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);
    double dx_deta_f0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
    // at face 1
    double dy_dxi_f1  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_dxi_f1  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_deta_f1 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);
    double dx_deta_f1 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
    // at face 2
    double dy_dxi_f2  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_s]);
    double dx_dxi_f2  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_s]);
    double dy_deta_f2 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);
    double dx_deta_f2 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
    // at face 3
    double dy_dxi_f3  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
    double dx_dxi_f3  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
    double dy_deta_f3 = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);
    double dx_deta_f3 = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);
    // at face 4
    double dy_dxi_f4  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
    double dx_dxi_f4  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
    double dy_deta_f4 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);
    double dx_deta_f4 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
    // at face 5
    double dy_dxi_f5  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
    double dx_dxi_f5  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
    double dy_deta_f5 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);
    double dx_deta_f5 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
    // at face 6
    double dy_dxi_f6  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
    double dx_dxi_f6  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
    double dy_deta_f6 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);
    double dx_deta_f6 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
    // at face 7
    double dy_dxi_f7  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
    double dx_dxi_f7  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
    double dy_deta_f7 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);
    double dx_deta_f7 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // No contribution to the continuity equation
    //
    //--------------------------------------------------------------------------
    // q-momentum equation
    // 
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
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s , scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta);

    // scv_0 face_1
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w , scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta);
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta);

    // scv_1 face_3
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_se], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s , scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta);

    // sub control volume 2 ============================================
    // scv_2 face_4
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e , scvf_fac * 3.* convection_J_11(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n , scvf_fac * 1.* convection_J_11(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3.* convection_J_12(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1.* convection_J_12(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3.* convection_J_13(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1.* convection_J_13(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta);

    // scv_2 face_5
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_ne], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta);

    // sub control volume 3 ============================================
    // scv_3 face_6
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w , scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta);

    // scv_3 face_7
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w , scvf_fac * 3. * convection_J_11(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n , scvf_fac * 1. * convection_J_11(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3. * convection_J_12(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1. * convection_J_12(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3. * convection_J_13(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1. * convection_J_13(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta);

    //--------------------------------------------------------------------------
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
    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s , scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta);

    // scv_0 face_1
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w , scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta);
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta);

    // scv_1 face_3
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_se], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s , scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta);

    // sub control volume 2 ============================================
    // scv_2 face_4
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e , scvf_fac * 3.* convection_J_21(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n , scvf_fac * 1.* convection_J_21(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3.* convection_J_22(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1.* convection_J_22(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3.* convection_J_23(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1.* convection_J_23(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta);

    // scv_2 face_5
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_ne], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta);

    // sub control volume 3 ============================================
    // scv_3 face_6
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w , scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta));

    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta);

    // scv_3 face_7
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w , scvf_fac * 3. * convection_J_21(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n , scvf_fac * 1. * convection_J_21(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3. * convection_J_22(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1. * convection_J_22(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3. * convection_J_23(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1. * convection_J_23(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta));
    
    scvf_fac = 0.5;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta);
//------------------------------------------------------------------------------

    return 0;
}
int convection_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    size_t nx, size_t ny)                          // RHS vector [h, q, r]^{n}
{
    // Convection for post processing; WITHOUT integration over the control volumes, just the line-integral.
    double h;
    double q;
    double r;
    double nxi;
    double neta;
    double nxi_dl;
    double neta_dl;

    std::fill_n(rhs_q.data(), rhs_q.size(), 0.0);
    std::fill_n(rhs_r.data(), rhs_r.size(), 0.0);

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            size_t p_0  = convection_idx(i    , j    , ny); // central point of control volume
            size_t p_sw = convection_idx(i - 1, j - 1, ny);  
            size_t p_s  = convection_idx(i    , j - 1, ny);  
            size_t p_se = convection_idx(i + 1, j - 1, ny);  
            size_t p_w  = convection_idx(i - 1, j    , ny);  
            size_t p_e  = convection_idx(i + 1, j    , ny);  
            size_t p_nw = convection_idx(i - 1, j + 1, ny);  
            size_t p_n  = convection_idx(i    , j + 1, ny);  
            size_t p_ne = convection_idx(i + 1, j + 1, ny);  

            // at face 0
            double dy_dxi_f0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
            double dx_dxi_f0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
            double dy_deta_f0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);
            double dx_deta_f0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
            // at face 1
            double dy_dxi_f1  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
            double dx_dxi_f1  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
            double dy_deta_f1 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);
            double dx_deta_f1 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
            // at face 2
            double dy_dxi_f2  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_s]);
            double dx_dxi_f2  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_s]);
            double dy_deta_f2 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);
            double dx_deta_f2 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
            // at face 3
            double dy_dxi_f3  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
            double dx_dxi_f3  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
            double dy_deta_f3 = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);
            double dx_deta_f3 = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);
            // at face 4
            double dy_dxi_f4  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
            double dx_dxi_f4  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
            double dy_deta_f4 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);
            double dx_deta_f4 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
            // at face 5
            double dy_dxi_f5  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
            double dx_dxi_f5  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
            double dy_deta_f5 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);
            double dx_deta_f5 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
            // at face 6
            double dy_dxi_f6  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
            double dx_dxi_f6  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
            double dy_deta_f6 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);
            double dx_deta_f6 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
            // at face 7
            double dy_dxi_f7  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
            double dx_dxi_f7  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
            double dy_deta_f7 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);
            double dx_deta_f7 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);

            // scv_0 face_0
            h = scvf_xi(hn[p_0], hn[p_w], hn[p_sw], hn[p_s]);
            q = scvf_xi(qn[p_0], qn[p_w], qn[p_sw], qn[p_s]);
            r = scvf_xi(rn[p_0], rn[p_w], rn[p_sw], rn[p_s]);

            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f0, dx_deta_f0, dy_dxi_f0, dy_deta_f0, nxi, neta);

            // scv_0 face_1
            h = scvf_eta(hn[p_0], hn[p_s], hn[p_sw], hn[p_w]);
            q = scvf_eta(qn[p_0], qn[p_s], qn[p_sw], qn[p_w]);
            r = scvf_eta(rn[p_0], rn[p_s], rn[p_sw], rn[p_w]);

            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f1, dx_deta_f1, dy_dxi_f1, dy_deta_f1, nxi, neta);

            // sub control volume 1 ============================================
            // scv_1 face_2
            h = scvf_eta(hn[p_0], hn[p_s], hn[p_se], hn[p_e]);
            q = scvf_eta(qn[p_0], qn[p_s], qn[p_se], qn[p_e]);
            r = scvf_eta(rn[p_0], rn[p_s], rn[p_se], rn[p_e]);

            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f2, dx_deta_f2, dy_dxi_f2, dy_deta_f2, nxi, neta);

            // scv_1 face_3
            h = scvf_xi(hn[p_0], hn[p_e], hn[p_se], hn[p_s]);
            q = scvf_xi(qn[p_0], qn[p_e], qn[p_se], qn[p_s]);
            r = scvf_xi(rn[p_0], rn[p_e], rn[p_se], rn[p_s]);

            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f3, dx_deta_f3, dy_dxi_f3, dy_deta_f3, nxi, neta);

            // sub control volume 2 ============================================
            // scv_2 face_4
            h = scvf_xi(hn[p_0], hn[p_e], hn[p_ne], hn[p_n]);
            q = scvf_xi(qn[p_0], qn[p_e], qn[p_ne], qn[p_n]);
            r = scvf_xi(rn[p_0], rn[p_e], rn[p_ne], rn[p_n]);

            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f4, dx_deta_f4, dy_dxi_f4, dy_deta_f4, nxi, neta);

            // scv_2 face_5
            h = scvf_eta(hn[p_0], hn[p_n], hn[p_ne], hn[p_e]);
            q = scvf_eta(qn[p_0], qn[p_n], qn[p_ne], qn[p_e]);
            r = scvf_eta(rn[p_0], rn[p_n], rn[p_ne], rn[p_e]);

            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f5, dx_deta_f5, dy_dxi_f5, dy_deta_f5, nxi, neta);

            // sub control volume 3 ============================================
            // scv_3 face_6
            h = scvf_eta(hn[p_0], hn[p_n], hn[p_nw], hn[p_w]);
            q = scvf_eta(qn[p_0], qn[p_n], qn[p_nw], qn[p_w]);
            r = scvf_eta(rn[p_0], rn[p_n], rn[p_nw], rn[p_w]);

            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f6, dx_deta_f6, dy_dxi_f6, dy_deta_f6, nxi, neta);

            // scv_3 face_7
            h = scvf_xi(hn[p_0], hn[p_w], hn[p_nw], hn[p_n]);
            q = scvf_xi(qn[p_0], qn[p_w], qn[p_nw], qn[p_n]);
            r = scvf_xi(rn[p_0], rn[p_w], rn[p_nw], rn[p_n]);

            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl  * convection_J_10(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, dx_dxi_f7, dx_deta_f7, dy_dxi_f7, dy_deta_f7, nxi, neta);
        }
    }
    return 0;
}
inline size_t convection_idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}

inline void add_value(double * values, size_t col, double data){ 
    values[col] += data; 
}


