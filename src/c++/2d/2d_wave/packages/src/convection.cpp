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

int convection_matrix_and_rhs(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double theta, double dx, double dy, int nx, int ny)
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

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

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
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // scv_0 face_1
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // scv_1 face_3
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_se], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // sub control volume 2 ============================================
    // scv_2 face_4
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 3.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 1.* convection_J_11(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1.* convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1.* convection_J_13(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // scv_2 face_5
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_ne], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // sub control volume 3 ============================================
    // scv_3 face_6
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

    // scv_3 face_7
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 3. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 1. * convection_J_11(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1. * convection_J_12(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1. * convection_J_13(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += -scvf_fac * convection_J_10(h, q, r, nxi, neta);

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
    scvf_fac = theta * 0.5 * dy * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dy;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // scv_0 face_1
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_sw, scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_sw + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_sw + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dx;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    h = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    nxi =  0.0;
    neta = -1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dx;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // scv_1 face_3
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_se], qtheta[p_s]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_se, scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_s , scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_se + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_s  + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_se + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_s  + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dy;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // sub control volume 2 ============================================
    // scv_2 face_4
    h = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    nxi =  1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 3.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 1.* convection_J_21(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 3.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1.* convection_J_22(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 3.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1.* convection_J_23(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dy;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // scv_2 face_5
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_ne], qtheta[p_e]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_ne, scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_e , scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_ne + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_e  + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_ne + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_e  + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dx;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // sub control volume 3 ============================================
    // scv_3 face_6
    h = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    nxi =  0.0;
    neta =  1.0;
    scvf_fac = theta * 0.5 * dx * 0.125;

    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));

    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));

    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));

    scvf_fac = 0.5 * dx;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);

    // scv_3 face_7
    h = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    q = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    r = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    nxi = -1.0;
    neta =  0.0;
    scvf_fac = theta * 0.5 * dy * 0.125;
    
    add_value(values, col_0 , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_w , scvf_fac * 3. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_nw, scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    add_value(values, col_n , scvf_fac * 1. * convection_J_21(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_w  + 1, scvf_fac * 3. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_nw + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    add_value(values, col_n  + 1, scvf_fac * 1. * convection_J_22(h, q, r, nxi, neta));
    
    add_value(values, col_0  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_w  + 2, scvf_fac * 3. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_nw + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    add_value(values, col_n  + 2, scvf_fac * 1. * convection_J_23(h, q, r, nxi, neta));
    
    scvf_fac = 0.5 * dy;
    rhs[row + 2] += -scvf_fac * convection_J_20(h, q, r, nxi, neta);
//------------------------------------------------------------------------------

    return 0;
}
int convection_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double dx, double dy, int nx, int ny)                          // RHS vector [h, q, r]^{n}
{
    // Convection for post processing; WITHOUT integration over the control volumes, just the line-integral.
    double h;
    double q;
    double r;
    double nxi;
    double neta;
    double nxi_dl;
    double neta_dl;

    memset(rhs_q.data(), 0, rhs_q.size() * sizeof(double));
    memset(rhs_r.data(), 0, rhs_r.size() * sizeof(double));

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int p_0  = convection_idx(i    , j    , ny); // central point of control volume
            int p_sw = convection_idx(i - 1, j - 1, ny);  
            int p_s  = convection_idx(i    , j - 1, ny);  
            int p_se = convection_idx(i + 1, j - 1, ny);  
            int p_w  = convection_idx(i - 1, j    , ny);  
            int p_e  = convection_idx(i + 1, j    , ny);  
            int p_nw = convection_idx(i - 1, j + 1, ny);  
            int p_n  = convection_idx(i    , j + 1, ny);  
            int p_ne = convection_idx(i + 1, j + 1, ny);  

            // scv_0 face_0
            h = scvf_xi(hn[p_0], hn[p_w], hn[p_sw], hn[p_s]);
            q = scvf_xi(qn[p_0], qn[p_w], qn[p_sw], qn[p_s]);
            r = scvf_xi(rn[p_0], rn[p_w], rn[p_sw], rn[p_s]);

            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // scv_0 face_1
            h = scvf_eta(hn[p_0], hn[p_s], hn[p_sw], hn[p_w]);
            q = scvf_eta(qn[p_0], qn[p_s], qn[p_sw], qn[p_w]);
            r = scvf_eta(rn[p_0], rn[p_s], rn[p_sw], rn[p_w]);

            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // sub control volume 1 ============================================
            // scv_1 face_2
            h = scvf_eta(hn[p_0], hn[p_s], hn[p_se], hn[p_e]);
            q = scvf_eta(qn[p_0], qn[p_s], qn[p_se], qn[p_e]);
            r = scvf_eta(rn[p_0], rn[p_s], rn[p_se], rn[p_e]);

            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // scv_1 face_3
            h = scvf_xi(hn[p_0], hn[p_e], hn[p_se], hn[p_s]);
            q = scvf_xi(qn[p_0], qn[p_e], qn[p_se], qn[p_s]);
            r = scvf_xi(rn[p_0], rn[p_e], rn[p_se], rn[p_s]);

            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // sub control volume 2 ============================================
            // scv_2 face_4
            h = scvf_xi(hn[p_0], hn[p_e], hn[p_ne], hn[p_n]);
            q = scvf_xi(qn[p_0], qn[p_e], qn[p_ne], qn[p_n]);
            r = scvf_xi(rn[p_0], rn[p_e], rn[p_ne], rn[p_n]);

            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // scv_2 face_5
            h = scvf_eta(hn[p_0], hn[p_n], hn[p_ne], hn[p_e]);
            q = scvf_eta(qn[p_0], qn[p_n], qn[p_ne], qn[p_e]);
            r = scvf_eta(rn[p_0], rn[p_n], rn[p_ne], rn[p_e]);

            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // sub control volume 3 ============================================
            // scv_3 face_6
            h = scvf_eta(hn[p_0], hn[p_n], hn[p_nw], hn[p_w]);
            q = scvf_eta(qn[p_0], qn[p_n], qn[p_nw], qn[p_w]);
            r = scvf_eta(rn[p_0], rn[p_n], rn[p_nw], rn[p_w]);

            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);

            // scv_3 face_7
            h = scvf_xi(hn[p_0], hn[p_w], hn[p_nw], hn[p_n]);
            q = scvf_xi(qn[p_0], qn[p_w], qn[p_nw], qn[p_n]);
            r = scvf_xi(rn[p_0], rn[p_w], rn[p_nw], rn[p_n]);

            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;
            rhs_q[p_0] += nxi_dl * convection_J_10(h, q, r, nxi, neta);
            rhs_r[p_0] += neta_dl * convection_J_20(h, q, r, nxi, neta);
        }
    }
    return 0;
}
inline int convection_idx(int i, int j, int ny)
{
    return i * ny + j;
}

inline void add_value(double * values, int col, double data){ 
    values[col] += data; 
}


