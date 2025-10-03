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

int viscosity_matrix_and_rhs(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& visc, double theta, double dx, double dy, int nx, int ny)
{
    double n_xi = 0.0;  // xi-component of the outward normal vector
    double n_eta = 0.0;  // eta-component of the outward normal vector
    double scvf_fac;
    double dxinv = 1./dx;
    double dyinv = 1./dy;

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
    {
    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta =  0.0;

    double visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    double htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    double qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    double dhtheta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    double dqtheta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);

    double aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    double bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    double cc_1 =  2. * -visc_0;
    double dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_w , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_sw, scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_s , scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);

    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_w  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_sw + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_s  + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_w  + 2, 0.0);
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_sw + 2, 0.0);
    add_value(values, col_s  + 2, 0.0);

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_sw], qtheta[p_w]);
    dhtheta_0 = dcdx_scvf_t(htheta[p_0], htheta[p_w], htheta[p_s], htheta[p_sw]);
    dqtheta_0 = dcdx_scvf_t(qtheta[p_0], qtheta[p_w], qtheta[p_s], qtheta[p_sw]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_s , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_sw, scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_w , scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 3. * 0.25 * dxinv) * n_xi);
    add_value(values, col_s  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 3. * 0.25 * dxinv) * n_xi);
    add_value(values, col_sw + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 1. * 0.25 * dxinv) * n_xi);
    add_value(values, col_w  + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 1. * 0.25 * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_s  + 2, 0.0);
    add_value(values, col_sw + 2, 0.0);
    add_value(values, col_w  + 2, 0.0);
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    dhtheta_0 = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqtheta_0 = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);

    aa_1 = 2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_s , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_se, scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_e , scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_s  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_se + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_e  + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_s  + 2, 0.0);
    add_value(values, col_se + 2, 0.0);
    add_value(values, col_e  + 2, 0.0);
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_s], visc[p_se], visc[p_e]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_se], qtheta[p_s]);
    dhtheta_0 = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_se], htheta[p_s]);
    dqtheta_0 = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_se], qtheta[p_s]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_e , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_se, scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_s , scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);

    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_e  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_se + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_s  + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_e  + 2, 0.0);
    add_value(values, col_se + 2, 0.0);
    add_value(values, col_s  + 2, 0.0);

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    dhtheta_0 = dcdx_scvf_n(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqtheta_0 = dcdx_scvf_n(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_e , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_ne, scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_n , scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);

    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_e  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_ne + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_n  + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_e  + 2, 0.0);
    add_value(values, col_ne + 2, 0.0);
    add_value(values, col_n  + 2, 0.0);

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    dhtheta_0 = dcdx_scvf_t(htheta[p_e], htheta[p_0], htheta[p_ne], htheta[p_n]);
    dqtheta_0 = dcdx_scvf_t(qtheta[p_e], qtheta[p_0], qtheta[p_ne], qtheta[p_n]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_n , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_ne, scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_e , scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_n  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_ne + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_e  + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_n  + 2, 0.0);
    add_value(values, col_ne + 2, 0.0);
    add_value(values, col_e  + 2, 0.0);
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_n], visc[p_nw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    qtheta_0 = scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    dhtheta_0 = dcdx_scvf_t(htheta[p_w], htheta[p_0], htheta[p_n], htheta[p_nw]);
    dqtheta_0 = dcdx_scvf_t(qtheta[p_w], qtheta[p_0], qtheta[p_n], qtheta[p_nw]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_n , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_nw, scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_w , scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_n  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_nw + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_w  + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_n  + 2, 0.0);
    add_value(values, col_nw + 2, 0.0);
    add_value(values, col_w  + 2, 0.0);
    
    scvf_fac = 0.5 * dx;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );

    // scv_3 face_7
    n_xi =  -1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    qtheta_0 = scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_nw], qtheta[p_n]);
    dhtheta_0 = dcdx_scvf_n(htheta[p_0], htheta[p_w], htheta[p_n], htheta[p_nw]);
    dqtheta_0 = dcdx_scvf_n(qtheta[p_0], qtheta[p_w], qtheta[p_n], qtheta[p_nw]);

    aa_1 =  2. * -visc_0 * qtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_1 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_1 =  2. * -visc_0;
    dd_1 = -2. * -visc_0 * qtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_1  + dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_w , scvf_fac * (0.125 * 3.* aa_1  - dd_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_nw, scvf_fac * (0.125 * 1.* aa_1  - dd_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_n , scvf_fac * (0.125 * 1.* aa_1  + dd_1 * 0.25 * 1. * dxinv) * n_xi);

    add_value(values, col_0  + 1, scvf_fac * (0.125 * 3.* bb_1 + cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_w  + 1, scvf_fac * (0.125 * 3.* bb_1 - cc_1 * 0.25 * 3. * dxinv) * n_xi);
    add_value(values, col_nw + 1, scvf_fac * (0.125 * 1.* bb_1 - cc_1 * 0.25 * 1. * dxinv) * n_xi);
    add_value(values, col_n  + 1, scvf_fac * (0.125 * 1.* bb_1 + cc_1 * 0.25 * 1. * dxinv) * n_xi);
    
    add_value(values, col_0  + 2, 0.0);
    add_value(values, col_w  + 2, 0.0);
    add_value(values, col_nw + 2, 0.0);
    add_value(values, col_n  + 2, 0.0);

    scvf_fac = 0.5 * dy;
    rhs[row + 1] += - n_xi * ( 2. * -visc_0 * dxinv * dqtheta_0 - dd_1 * dxinv * dhtheta_0 );
    }
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

    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta =  0.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_w], visc[p_sw], visc[p_s]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_sw, scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w , scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_sw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_sw + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w  + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_sw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_sw], htheta[p_w]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_sw], rtheta[p_w]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_0], htheta[p_s], htheta[p_w], htheta[p_sw]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_w], rtheta[p_sw]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_sw, scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w , scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_sw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_sw + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w  + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    visc_0 = scvf_eta(visc[p_0], visc[p_s], visc[p_e], visc[p_se]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_se, scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_e , scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);
    add_value(values, col_se + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_s  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_se + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_e  + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_se], visc[p_s]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_se], htheta[p_s]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_se], rtheta[p_s]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_0], htheta[p_s], htheta[p_e], htheta[p_se]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_0], rtheta[p_s], rtheta[p_e], rtheta[p_se]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_e , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_se, scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_s , scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);
    add_value(values, col_se + 1, 0.0);
    add_value(values, col_s  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_e  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_se + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_s  + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    visc_0 = scvf_xi(visc[p_0], visc[p_e], visc[p_ne], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_e , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_ne, scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_n , scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);
    add_value(values, col_ne + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_e  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_ne + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_n  + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_n], visc[p_ne], visc[p_e]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_ne], htheta[p_e]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_ne], rtheta[p_e]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_ne], htheta[p_e]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_ne], rtheta[p_e]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_n , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_ne, scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_e , scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);
    add_value(values, col_ne + 1, 0.0);
    add_value(values, col_e  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_n  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_ne + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_e  + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    visc_0   = scvf_eta(visc[p_0], visc[p_n], visc[p_nw], visc[p_w]);
    htheta_0 = scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    rtheta_0 = scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    dhtheta_0 = dcdy_scvf_n(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    drtheta_0 = dcdy_scvf_n(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dx;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_n , scvf_fac * (0.125 * 3.* aa_4  + dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_nw, scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w , scvf_fac * (0.125 * 1.* aa_4  - dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);
    add_value(values, col_nw + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_n  + 2, scvf_fac * (0.125 * 3.* bb_4 + cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_nw + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_w  + 2, scvf_fac * (0.125 * 1.* bb_4 - cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );

    // scv_3 face_7
    n_xi = -1.0;
    n_eta =  0.0;

    visc_0   = scvf_xi(visc[p_0], visc[p_w], visc[p_nw], visc[p_n]);
    htheta_0 = scvf_xi(htheta[p_0], htheta[p_w], htheta[p_nw], htheta[p_n]);
    rtheta_0 = scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_nw], rtheta[p_n]);
    dhtheta_0 = dcdy_scvf_t(htheta[p_n], htheta[p_0], htheta[p_nw], htheta[p_w]);
    drtheta_0 = dcdy_scvf_t(rtheta[p_n], rtheta[p_0], rtheta[p_nw], rtheta[p_w]);

    aa_4 =  2. * -visc_0 * rtheta_0 /(htheta_0 * htheta_0) * dhtheta_0;
    bb_4 = -2. * -visc_0 / htheta_0 * dhtheta_0;
    cc_4 =  2. * -visc_0;
    dd_4 = -2. * -visc_0 * rtheta_0/htheta_0;

    scvf_fac = theta * 0.5 * dy;
    add_value(values, col_0 , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_w , scvf_fac * (0.125 * 3.* aa_4  - dd_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_nw, scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_n , scvf_fac * (0.125 * 1.* aa_4  + dd_4 * 0.25 * 1. * dyinv) * n_eta);
   
    add_value(values, col_0  + 1, 0.0);
    add_value(values, col_w  + 1, 0.0);
    add_value(values, col_nw + 1, 0.0);
    add_value(values, col_n  + 1, 0.0);

    add_value(values, col_0  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_w  + 2, scvf_fac * (0.125 * 3.* bb_4 - cc_4 * 0.25 * 3. * dyinv) * n_eta);
    add_value(values, col_nw + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);
    add_value(values, col_n  + 2, scvf_fac * (0.125 * 1.* bb_4 + cc_4 * 0.25 * 1. * dyinv) * n_eta);

    scvf_fac = 0.5 * dy;
    rhs[row + 2] += - n_eta * ( 2. * -visc_0 * dyinv * drtheta_0 - dd_4 * dyinv * dhtheta_0 );
    }
//------------------------------------------------------------------------------

    return 0;
}
int viscosity_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& visc, double dx, double dy, int nx, int ny)                          // RHS vector [h, q, r]^{n}
{
    // Viscosity for post processing; WITHOUT integration over the control volumes, just the line-integral.
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
            int p_0  = viscosity_idx(i    , j    , ny); // central point of control volume
            int p_sw = viscosity_idx(i - 1, j - 1, ny);  
            int p_s  = viscosity_idx(i    , j - 1, ny);  
            int p_se = viscosity_idx(i + 1, j - 1, ny);  
            int p_w  = viscosity_idx(i - 1, j    , ny);  
            int p_e  = viscosity_idx(i + 1, j    , ny);  
            int p_nw = viscosity_idx(i - 1, j + 1, ny);  
            int p_n  = viscosity_idx(i    , j + 1, ny);  
            int p_ne = viscosity_idx(i + 1, j + 1, ny);  

            // scv_0 face_0
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // scv_0 face_1
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // sub control volume 1 ============================================
            // scv_1 face_2
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // scv_1 face_3
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // sub control volume 2 ============================================
            // scv_2 face_4
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // scv_2 face_5
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // sub control volume 3 ============================================
            // scv_3 face_6
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = 0.5 * dx;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;

            // scv_3 face_7
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            rhs_q[p_0] += nxi_dl * 0.0;
            rhs_r[p_0] += neta_dl * 0.0;
        }
    }
    return 0;
}
inline int viscosity_idx(int i, int j, int ny)
{
    return i * ny + j;
}

inline void add_value(double * values, int col, double data){ 
    values[col] += data; 
}


