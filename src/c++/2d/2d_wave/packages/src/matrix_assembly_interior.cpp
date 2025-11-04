//
// programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formuation and Modified Newton iteration 
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

#include "interpolations.h"
#include "matrix_assembly_interior.h"

//------------------------------------------------------------------------------
//   nw - - - - - - - n - - - - - - - ne          2 - - - - - - - 5 - - - - - - - 8
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    w - - - - - - - c - - - - - - - e           1 - - - - - - - 4 - - - - - - - 7
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//   sw - - - - - - - s - - - - - - - se          0 - - - - - - - 3 - - - - - - - 6

//
//corner nodes
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

int interior(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, bool do_convection, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, std::vector<double>& mass)
{
    int p_0 = c_eq/(3*27);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (std::fmod(p_0, ny) == 0) { return 1; }  // south boundary
    if (std::fmod(p_0 + 1, ny) == 0) { return 2; }  // north boundary

    memset(&values[c_eq], 0, 3 * 27 * sizeof(double));  // set all coefficients for one row of Delta c-, Delta q- and Delta r-equation to zero

    // std::vector<int> p;
    // p.push_back(p_0 - ny - 1);
    // p.push_back(p_0 - ny);
    // p.push_back(p_0 - ny + 1);
    // p.push_back(p_0 - 1);
    // p.push_back(p_0);
    // p.push_back(p_0 + 1);
    // p.push_back(p_0 + ny - 1);
    // p.push_back(p_0 + ny);
    // p.push_back(p_0 + ny + 1);

    size_t p_sw = p_0 - ny - 1;
    size_t p_w  = p_0 - ny;
    size_t p_nw = p_0 - ny + 1;
    size_t p_s  = p_0 - 1; 
    size_t p_n  = p_0 + 1;
    size_t p_se = p_0 + ny - 1;
    size_t p_e  = p_0 + ny;
    size_t p_ne = p_0 + ny + 1;

    double htheta_0  = htheta[p_0 ];
    double htheta_sw = htheta[p_sw];
    double htheta_w  = htheta[p_w ];
    double htheta_nw = htheta[p_nw];
    double htheta_s  = htheta[p_s ];
    double htheta_n  = htheta[p_n ];
    double htheta_se = htheta[p_se];
    double htheta_e  = htheta[p_e ];
    double htheta_ne = htheta[p_ne];

    double qtheta_0  = qtheta[p_0 ];
    double qtheta_sw = qtheta[p_sw];
    double qtheta_w  = qtheta[p_w ];
    double qtheta_nw = qtheta[p_nw];
    double qtheta_s  = qtheta[p_s ];
    double qtheta_n  = qtheta[p_n ];
    double qtheta_se = qtheta[p_se];
    double qtheta_e  = qtheta[p_e ];
    double qtheta_ne = qtheta[p_ne];

    double rtheta_0  = rtheta[p_0 ];
    double rtheta_sw = rtheta[p_sw];
    double rtheta_w  = rtheta[p_w ];
    double rtheta_nw = rtheta[p_nw];
    double rtheta_s  = rtheta[p_s ];
    double rtheta_n  = rtheta[p_n ];
    double rtheta_se = rtheta[p_se];
    double rtheta_e  = rtheta[p_e ];
    double rtheta_ne = rtheta[p_ne];

    int col_sw = c_eq;
    int col_w  = c_eq + 3;
    int col_nw = c_eq + 6;
    int col_s  = c_eq + 9;
    int col_0  = c_eq + 12;
    int col_n  = c_eq + 15;
    int col_se = c_eq + 18;
    int col_e  = c_eq + 21;
    int col_ne = c_eq + 24;

    //==========================================================================
    // c-equation
    //==========================================================================
    // 
    rhs[row] = 0.0;
    // scv_0
    std::vector<double> x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
    std::vector<double> y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
    double scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_w , dtinv * scv_area * 3./16.);
    add_value(values, col_s , dtinv * scv_area * 3./16.);
    add_value(values, col_sw, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_w] - hn[p_w], hp[p_s] - hn[p_s], hp[p_sw] - hn[p_sw]);
    //
    // scv 1
    x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
    y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_s , dtinv * scv_area * 3./16.);
    add_value(values, col_e , dtinv * scv_area * 3./16.);
    add_value(values, col_se, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_s] - hn[p_s], hp[p_e] - hn[p_e], hp[p_se] - hn[p_se]);
    //
    // scv 2
    x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
    y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_e , dtinv * scv_area * 3./16.);
    add_value(values, col_n , dtinv * scv_area * 3./16.);
    add_value(values, col_ne, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_e] - hn[p_e], hp[p_n] - hn[p_n], hp[p_ne] - hn[p_ne]);
    //
    //scv 3
    x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
    y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_n , dtinv * scv_area * 3./16.);
    add_value(values, col_w , dtinv * scv_area * 3./16.);
    add_value(values, col_nw, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_n] - hn[p_n], hp[p_w] - hn[p_w], hp[p_nw] - hn[p_nw]);

    //
    // Mass flux
    // 
    double fluxx_0 = 0.0;
    double fluxy_1 = 0.0;
    double fluxy_2 = 0.0;
    double fluxx_3 = 0.0;
    double fluxx_4 = 0.0;
    double fluxy_5 = 0.0;
    double fluxy_6 = 0.0;
    double fluxx_7 = 0.0;

    // d(q)/dxi
    // 
    // sub control volume 0 ============================================
    // 
    // scv_0 face_0
    double dy_dxi  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_dxi  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_deta = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);
    double dx_deta = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);

    double n_xi = -1.0;
    double n_eta = 0.0;

    double dl_y = 0.5 * n_xi;  // length in y-direction
    double dl_x = 0.5 * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_w  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_sw + 1, theta * dl_y * dy_deta * 1./8.);
    add_value(values, col_s  + 1, theta * dl_y * dy_deta * 1./8.);

    add_value(values, col_0  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_w  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_sw + 2, -theta * dl_x * dx_deta * 1./8.);
    add_value(values, col_s  + 2, -theta * dl_x * dx_deta * 1./8.);

    // rhs
    fluxx_0 = dy_deta * (scvf_xi(qtheta_0, qtheta_w, qtheta_s, qtheta_sw)) * dl_y 
            - dx_deta * (scvf_xi(rtheta_0, rtheta_w, rtheta_s, rtheta_sw)) * dl_x;
    //
    // scv_0 face_1
    // No contribution to Delta r and rhs

    // sub control volume 1 ============================================
    // scv_1 face_2
    // No contribution to Delta r and rhs
    // scv_1 face_3
    dy_dxi  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
    dx_dxi  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
    dy_deta = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);
    dx_deta = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);

    n_xi = 1.0;
    n_eta = 0.0;

    dl_y = 0.5 * n_xi;  // length in y-direction
    dl_x = 0.5 * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_e  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_se + 1, theta * dl_y * dy_deta * 1./8.);
    add_value(values, col_s  + 1, theta * dl_y * dy_deta * 1./8.);

    add_value(values, col_0  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_s  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_se + 2, -theta * dl_x * dx_deta * 1./8.);
    add_value(values, col_e  + 2, -theta * dl_x * dx_deta * 1./8.);

    // rhs
    fluxx_3 = dy_deta * (scvf_xi(qtheta_e, qtheta_0, qtheta_se, qtheta_s)) * dl_y 
            - dx_deta * (scvf_xi(rtheta_e, rtheta_0, rtheta_se, rtheta_s)) * dl_x;
    
    // sub control volume 2 ============================================
    // scv_2 face_4
    dy_dxi  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_dxi  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_deta = dcdx_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);
    dx_deta = dcdx_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);

    n_xi =  1.0;
    n_eta = 0.0;

    dl_y = 0.5 * n_xi;  // length in y-direction
    dl_x = 0.5 * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_e  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_ne + 1, theta * dl_y * dy_deta * 1./8.);
    add_value(values, col_n  + 1, theta * dl_y * dy_deta * 1./8.);

    add_value(values, col_0  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_e  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_ne + 2, -theta * dl_x * dx_deta * 1./8.);
    add_value(values, col_n  + 2, -theta * dl_x * dx_deta * 1./8.);

    fluxx_4 = dy_deta * (scvf_xi(qtheta_e, qtheta_0, qtheta_ne, qtheta_n)) * dl_y 
            - dx_deta * (scvf_xi(rtheta_e, rtheta_0, rtheta_ne, rtheta_n)) * dl_x;
    // scv_2 face_5
    // No contribution to Delta r and rhs
                                
    // sub control volume 3 ============================================
    // scv_3 face_6
    // No contribution to Delta r and rhs
    // scv_3 face_7
    // 
    dy_dxi  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_dxi  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_deta = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);
    dx_deta = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);

    n_xi = -1.0;
    n_eta = 0.0;

    dl_y = 0.5 * n_xi;  // length in y-direction
    dl_x = 0.5 * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_w  + 1, theta * dl_y * dy_deta * 3./8.);
    add_value(values, col_n  + 1, theta * dl_y * dy_deta * 1./8.);
    add_value(values, col_nw + 1, theta * dl_y * dy_deta * 1./8.);

    add_value(values, col_0  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_w  + 2, -theta * dl_x * dx_deta * 3./8.);
    add_value(values, col_n  + 2, -theta * dl_x * dx_deta * 1./8.);
    add_value(values, col_nw + 2, -theta * dl_x * dx_deta * 1./8.);

    fluxx_7 = dy_deta * (scvf_xi(qtheta_0, qtheta_w, qtheta_n, qtheta_nw)) * dl_y 
            - dx_deta * (scvf_xi(rtheta_0, rtheta_w, rtheta_n, rtheta_nw)) * dl_x;
    //
    // mass flux
    // 
    // d(r)/deta
    //
    // sub control volume 0 ============================================
    //
    // scv_0 face_0
    // No contribution to Delta q and rhs
    // scv_0 face_1
    dy_dxi  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
    dx_dxi  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
    dy_deta = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);
    dx_deta = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);

    n_xi = 0.0;
    n_eta = -1.0;

    dl_y = 0.5 * n_eta;  // length in y-direction
    dl_x = 0.5 * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_s  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_w  + 1, -theta * dl_y * dy_dxi * 1./8.);
    add_value(values, col_sw + 1, -theta * dl_y * dy_dxi * 1./8.);

    add_value(values, col_0  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_s  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_w  + 2, theta * dl_x * dx_dxi * 1./8.);
    add_value(values, col_sw + 2, theta * dl_x * dx_dxi * 1./8.);

    // rhs
    fluxy_1 = - dy_dxi * (scvf_eta(qtheta_0, qtheta_s, qtheta_sw, qtheta_w)) * dl_y 
              + dx_dxi * (scvf_eta(rtheta_0, rtheta_s, rtheta_sw, rtheta_w)) * dl_x;
    //                                        
    // sub control volume 1 ============================================
    // scv_1 face_2
    dy_dxi  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_s]);
    dx_dxi  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_s]);
    dy_deta = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);
    dx_deta = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
    //
    n_xi = 0.0;
    n_eta = -1.0;

    dl_y = 0.5 * n_eta;  // length in y-direction
    dl_x = 0.5 * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_s  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_se + 1, -theta * dl_y * dy_dxi * 1./8.);
    add_value(values, col_e  + 1, -theta * dl_y * dy_dxi * 1./8.);

    add_value(values, col_0  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_s  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_se + 2, theta * dl_x * dx_dxi * 1./8.);
    add_value(values, col_e  + 2, theta * dl_x * dx_dxi * 1./8.);

    // rhs
    fluxy_2 = - dy_dxi * (scvf_eta(qtheta_0, qtheta_s, qtheta_se, qtheta_e)) * dl_y 
              + dx_dxi * (scvf_eta(rtheta_0, rtheta_s, rtheta_se, rtheta_e)) * dl_x;
    // scv_1 face_3
    // No contribution to r-momentum equation
                            
    // sub control volume 2 ============================================
    // scv_1 face_4
    // No contribution to Delta q and rhs
    // scv_1 face_5
    dy_dxi  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
    dx_dxi  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
    dy_deta = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);
    dx_deta = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);

    n_xi =  0.0;
    n_eta = 1.0;

    dl_y = 0.5 * n_eta;  // length in y-direction
    dl_x = 0.5 * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_n  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_e  + 1, -theta * dl_y * dy_dxi * 1./8.);
    add_value(values, col_ne + 1, -theta * dl_y * dy_dxi * 1./8.);

    add_value(values, col_0  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_n  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_e  + 2, theta * dl_x * dx_dxi * 1./8.);
    add_value(values, col_ne + 2, theta * dl_x * dx_dxi * 1./8.);

    // rhs
    fluxy_5 = - dy_dxi * (scvf_eta(qtheta_n, qtheta_0, qtheta_ne, qtheta_e)) * dl_y 
              + dx_dxi * (scvf_eta(rtheta_n, rtheta_0, rtheta_ne, rtheta_e)) * dl_x;
                           
    // sub control volume 3 ============================================
    dy_dxi  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
    dx_dxi  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
    dy_deta = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);
    dx_deta = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);
    // scv_1 face_6
    n_xi =  0.0;
    n_eta = 1.0;

    dl_y = 0.5 * n_eta;  // length in y-direction
    dl_x = 0.5 * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_n  + 1, -theta * dl_y * dy_dxi * 3./8.);
    add_value(values, col_nw + 1, -theta * dl_y * dy_dxi * 1./8.);
    add_value(values, col_w  + 1, -theta * dl_y * dy_dxi * 1./8.);

    add_value(values, col_0  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_n  + 2, theta * dl_x * dx_dxi * 3./8.);
    add_value(values, col_nw + 2, theta * dl_x * dx_dxi * 1./8.);
    add_value(values, col_w  + 2, theta * dl_x * dx_dxi * 1./8.);

    // rhs
    fluxy_6 = - dy_dxi * (scvf_eta(qtheta_n, qtheta_0, qtheta_nw, qtheta_w)) * dl_y 
              + dx_dxi * (scvf_eta(rtheta_n, rtheta_0, rtheta_nw, rtheta_w)) * dl_x;
    // scv_1 face_7
    // No contribution to Delta q and rhs

    //double fluxx_0 = -0.5 * dy * scvf_xi (qtheta_0, qtheta_w, qtheta_s, qtheta_sw);
    //double fluxy_1 = -0.5 * dx * scvf_eta(rtheta_0, rtheta_s, rtheta_w, rtheta_sw);
    //double fluxx_3 =  0.5 * dy * scvf_xi (qtheta_0, qtheta_e, qtheta_s, qtheta_se);
    //double fluxx_4 =  0.5 * dy * scvf_xi (qtheta_0, qtheta_e, qtheta_n, qtheta_ne);
    //double fluxy_5 =  0.5 * dx * scvf_eta(rtheta_0, rtheta_n, rtheta_e, rtheta_ne);
    //double fluxy_6 =  0.5 * dx * scvf_eta(rtheta_0, rtheta_n, rtheta_w, rtheta_nw);
    //double fluxx_7 = -0.5 * dy * scvf_xi (qtheta_0, qtheta_w, qtheta_n, qtheta_nw);

    double flux = (fluxx_0 + fluxy_1 + fluxy_2 + fluxx_3 + fluxx_4 + fluxy_5 + fluxy_6 + fluxx_7);
    rhs[row] += -flux;

    //
    //==========================================================================
    // q-equation
    //==========================================================================
    // 
    col_sw = q_eq;
    col_w  = q_eq + 3;
    col_nw = q_eq + 6;
    col_s  = q_eq + 9;
    col_0  = q_eq + 12;
    col_n  = q_eq + 15;
    col_se = q_eq + 18;
    col_e  = q_eq + 21;
    col_ne = q_eq + 24;
    // 
    // time-derivative: d(q)/dt
    //
    rhs[row + 1] = 0.0;
    // scv_0
    x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
    y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 1, dtinv * scv_area * 9./16.);
    add_value(values, col_w  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_s  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_sw + 1, dtinv * scv_area * 1./16.);
    rhs[row + 1] += -dtinv * scv_area * c_scv(qp[p_0] - qn[p_0], qp[p_w] - qn[p_w], qp[p_s] - qn[p_s], qp[p_sw] - qn[p_sw]);
    //
    // scv 1
    x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
    y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 1, dtinv * scv_area * 9./16.);
    add_value(values, col_s  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_e  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_se + 1, dtinv * scv_area * 1./16.);
    rhs[row + 1] += -dtinv * scv_area * c_scv(qp[p_0] - qn[p_0], qp[p_s] - qn[p_s], qp[p_e] - qn[p_e], qp[p_se] - qn[p_se]);
    //
    // scv 2
    x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
    y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 1, dtinv * scv_area * 9./16.);
    add_value(values, col_e  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_n  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_ne + 1, dtinv * scv_area * 1./16.);
    rhs[row + 1] += -dtinv * scv_area * c_scv(qp[p_0] - qn[p_0], qp[p_e] - qn[p_e], qp[p_n] - qn[p_n], qp[p_ne] - qn[p_ne]);
    //
    //scv 3
    x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
    y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 1, dtinv * scv_area * 9./16.);
    add_value(values, col_n  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_w  + 1, dtinv * scv_area * 3./16.);
    add_value(values, col_nw + 1, dtinv * scv_area * 1./16.);
    rhs[row + 1] += -dtinv * scv_area * c_scv(qp[p_0] - qn[p_0], qp[p_n] - qn[p_n], qp[p_w] - qn[p_w], qp[p_nw] - qn[p_nw]);
    //
    // pressure term: gh ( y_eta d(zeta)/dxi - y_xi d(zeta)/deta ) 
    //
    // sub control volume 0 ============================================

    double depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
    double depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
    double depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
    double depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

    double dy_deta_0 = 0.25 * (3. * (y[p_0] - y[p_s]) + 1. * (y[p_w ] - y[p_sw]));
    double dy_deta_1 = 0.25 * (3. * (y[p_0] - y[p_s]) + 1. * (y[p_e ] - y[p_se]));
    double dy_deta_2 = 0.25 * (3. * (y[p_n] - y[p_0]) + 1. * (y[p_ne] - y[p_e ]));
    double dy_deta_3 = 0.25 * (3. * (y[p_n] - y[p_0]) + 1. * (y[p_nw] - y[p_w ]));

    double dy_dxi_0 = 0.25 * (3. * (y[p_0] - y[p_w]) + 1. * (y[p_s ] - y[p_sw]));
    double dy_dxi_1 = 0.25 * (3. * (y[p_e] - y[p_0]) + 1. * (y[p_se] - y[p_s ]));
    double dy_dxi_2 = 0.25 * (3. * (y[p_e] - y[p_0]) + 1. * (y[p_ne] - y[p_n ]));
    double dy_dxi_3 = 0.25 * (3. * (y[p_0] - y[p_w]) + 1. * (y[p_n ] - y[p_nw]));
    
    double dzetadxi_0 = dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_s  + zb[p_s ], htheta_sw + zb[p_sw]);
    double dzetadxi_1 = dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_se + zb[p_se], htheta_s  + zb[p_s ]);
    double dzetadxi_2 = dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_n  + zb[p_n ]);
    double dzetadxi_3 = dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_n  + zb[p_n ], htheta_nw + zb[p_nw]);
    
    double dzetadeta_0 = dcdy_scv(htheta_0 + zb[p_0], htheta_w + zb[p_s], htheta_s  + zb[p_w ], htheta_sw + zb[p_sw]);
    double dzetadeta_1 = dcdy_scv(htheta_0 + zb[p_0], htheta_w + zb[p_s], htheta_s  + zb[p_e ], htheta_sw + zb[p_se]);
    double dzetadeta_2 = dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_e  + zb[p_e ]);
    double dzetadeta_3 = dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_nw + zb[p_nw], htheta_w  + zb[p_w ]);
    //
    // theta * g ( dy_deta dzeta/dxi - dx_dxi dzeta/deta) Delta h
    //
    // Contribution to Delta h
    double fac = theta * 0.25 * g;
    // scv_0
    add_value(values, col_0 , fac * (dy_deta_0 * dzetadxi_0 * 9./16. + -dy_dxi_0 * dzetadeta_0 * 9./16.));   
    add_value(values, col_w , fac * (dy_deta_0 * dzetadxi_0 * 3./16. + -dy_dxi_0 * dzetadeta_0 * 3./16.));
    add_value(values, col_sw, fac * (dy_deta_0 * dzetadxi_0 * 1./16. + -dy_dxi_0 * dzetadeta_0 * 1./16.));
    add_value(values, col_s , fac * (dy_deta_0 * dzetadxi_0 * 3./16. + -dy_dxi_0 * dzetadeta_0 * 3./16.)); 
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (dy_deta_1 * dzetadxi_1 * 9./16. + -dy_dxi_1 * dzetadeta_1 * 9./16.));   
    add_value(values, col_s , fac * (dy_deta_1 * dzetadxi_1 * 3./16. + -dy_dxi_1 * dzetadeta_1 * 3./16.)); 
    add_value(values, col_se, fac * (dy_deta_1 * dzetadxi_1 * 1./16. + -dy_dxi_1 * dzetadeta_1 * 1./16.));
    add_value(values, col_e , fac * (dy_deta_1 * dzetadxi_1 * 3./16. + -dy_dxi_1 * dzetadeta_1 * 3./16.));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (dy_deta_2 * dzetadxi_2 * 9./16. + -dy_dxi_2 * dzetadeta_2 * 9./16.));   
    add_value(values, col_e , fac * (dy_deta_2 * dzetadxi_2 * 3./16. + -dy_dxi_2 * dzetadeta_2 * 3./16.));
    add_value(values, col_ne, fac * (dy_deta_2 * dzetadxi_2 * 1./16. + -dy_dxi_2 * dzetadeta_2 * 1./16.));
    add_value(values, col_n , fac * (dy_deta_2 * dzetadxi_2 * 3./16. + -dy_dxi_2 * dzetadeta_2 * 3./16.));   
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (dy_deta_3 * dzetadxi_3 * 9./16. + -dy_dxi_3 * dzetadeta_3 * 9./16.));   
    add_value(values, col_n , fac * (dy_deta_3 * dzetadxi_3 * 3./16. + -dy_dxi_3 * dzetadeta_3 * 3./16.));   
    add_value(values, col_nw, fac * (dy_deta_3 * dzetadxi_3 * 1./16. + -dy_dxi_3 * dzetadeta_3 * 1./16.));
    add_value(values, col_w , fac * (dy_deta_3 * dzetadxi_3 * 3./16. + -dy_dxi_3 * dzetadeta_3 * 3./16.)); 
    //
    // theta * g * h * ( dy_deta d(Delta zeta)/dxi  - dy_dxi d(Delta zeta)/deta)
    //
    // sub control volume 0 ============================================
    // scv_0
    fac = theta * 0.25 * g * 0.25;
    add_value(values, col_0 , fac * (+3. * dy_deta_0 * depth_0 +  3. * -dy_dxi_0 * depth_0));
    add_value(values, col_w , fac * (-3. * dy_deta_0 * depth_0 +  1. * -dy_dxi_0 * depth_0));
    add_value(values, col_sw, fac * (-1. * dy_deta_0 * depth_0 + -1. * -dy_dxi_0 * depth_0));
    add_value(values, col_s , fac * (+1. * dy_deta_0 * depth_0 + -3. * -dy_dxi_0 * depth_0));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (-3. * dy_deta_1 * depth_1 +  3. * -dy_dxi_1 * depth_1));
    add_value(values, col_s , fac * (-1. * dy_deta_1 * depth_1 + -3. * -dy_dxi_1 * depth_1));
    add_value(values, col_se, fac * (+1. * dy_deta_1 * depth_1 + -1. * -dy_dxi_1 * depth_1));
    add_value(values, col_e , fac * (+3. * dy_deta_1 * depth_1 +  1. * -dy_dxi_1 * depth_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (-3. * dy_deta_2 * depth_2 + -3. * -dy_dxi_2 * depth_2));
    add_value(values, col_e , fac * (+3. * dy_deta_2 * depth_2 + -1. * -dy_dxi_2 * depth_2));
    add_value(values, col_ne, fac * (+1. * dy_deta_2 * depth_2 +  1. * -dy_dxi_2 * depth_2));
    add_value(values, col_n , fac * (-1. * dy_deta_2 * depth_2 +  3. * -dy_dxi_2 * depth_2));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (+3. * dy_deta_3 * depth_3 + -3. * -dy_dxi_3 * depth_3));
    add_value(values, col_n , fac * (+1. * dy_deta_3 * depth_3 +  3. * -dy_dxi_3 * depth_3));
    add_value(values, col_nw, fac * (-1. * dy_deta_3 * depth_3 +  1. * -dy_dxi_3 * depth_3));
    add_value(values, col_w , fac * (-3. * dy_deta_3 * depth_3 + -1. * -dy_dxi_3 * depth_3));
    //
    // RHS q-momentum equation
    //
    rhs[row + 1] += - 0.25 * g * (
        depth_0 * (dy_deta_0 * dzetadxi_0 - dy_dxi_0 * dzetadeta_0) + 
        depth_1 * (dy_deta_1 * dzetadxi_1 - dy_dxi_1 * dzetadeta_1) +
        depth_2 * (dy_deta_2 * dzetadxi_2 - dy_dxi_2 * dzetadeta_2) +
        depth_3 * (dy_deta_3 * dzetadxi_3 - dy_dxi_3 * dzetadeta_3)
        );
    //
    //==========================================================================
    // r-equation
    //==========================================================================
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
    // time-derivative: d(r)/dt
    //
    rhs[row + 2] = 0.0;
    // scv_0
    x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
    y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 2, dtinv * scv_area * 9./16.);
    add_value(values, col_w  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_s  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_sw + 2, dtinv * scv_area * 1./16.);
    rhs[row + 2] += -dtinv * scv_area * c_scv(rp[p_0] - rn[p_0], rp[p_w] - rn[p_w], rp[p_s] - rn[p_s], rp[p_sw] - rn[p_sw]);
    //
    // scv 1
    x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
    y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 2, dtinv * scv_area * 9./16.);
    add_value(values, col_s  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_e  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_se + 2, dtinv * scv_area * 1./16.);
    rhs[row + 2] += -dtinv * scv_area * c_scv(rp[p_0] - rn[p_0], rp[p_s] - rn[p_s], rp[p_e] - rn[p_e], rp[p_se] - rn[p_se]);
    //
    // scv 2
    x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
    y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 2, dtinv * scv_area * 9./16.);
    add_value(values, col_e  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_n  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_ne + 2, dtinv * scv_area * 1./16.);
    rhs[row + 2] += -dtinv * scv_area * c_scv(rp[p_0] - rn[p_0], rp[p_e] - rn[p_e], rp[p_n] - rn[p_n], rp[p_ne] - rn[p_ne]);
    //
    //scv 3
    x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
    y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
    scv_area = polygon_area(x_pol, y_pol);
    add_value(values, col_0  + 2, dtinv * scv_area * 9./16.);
    add_value(values, col_n  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_w  + 2, dtinv * scv_area * 3./16.);
    add_value(values, col_nw + 2, dtinv * scv_area * 1./16.);
    rhs[row + 2] += -dtinv * scv_area * c_scv(rp[p_0] - rn[p_0], rp[p_n] - rn[p_n], rp[p_w] - rn[p_w], rp[p_nw] - rn[p_nw]);
    //
    // pressure term: gh ( -x_eta d(zeta)/dxi + x_xi d(zeta)/deta)
    //
    depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
    depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
    depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
    depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

    double dx_dxi_0 = 0.25 * (3. * (x[p_0] - x[p_w]) + 1. * (x[p_s ] - x[p_sw]));
    double dx_dxi_1 = 0.25 * (3. * (x[p_e] - x[p_0]) + 1. * (x[p_se] - x[p_s ]));
    double dx_dxi_2 = 0.25 * (3. * (x[p_e] - x[p_0]) + 1. * (x[p_ne] - x[p_n ]));
    double dx_dxi_3 = 0.25 * (3. * (x[p_0] - x[p_w]) + 1. * (x[p_n ] - x[p_nw]));

    double dx_deta_0 = 0.25 * (3. * (x[p_0] - x[p_s]) + 1. * (x[p_w ] - x[p_sw]));
    double dx_deta_1 = 0.25 * (3. * (x[p_0] - x[p_s]) + 1. * (x[p_e ] - x[p_se]));
    double dx_deta_2 = 0.25 * (3. * (x[p_n] - x[p_0]) + 1. * (x[p_ne] - x[p_e]));
    double dx_deta_3 = 0.25 * (3. * (x[p_n] - x[p_0]) + 1. * (x[p_nw] - x[p_w]));

    dzetadxi_0 = dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_s  + zb[p_s ], htheta_sw + zb[p_sw]);
    dzetadxi_1 = dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_se + zb[p_se], htheta_s  + zb[p_s ]);
    dzetadxi_2 = dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_n  + zb[p_n ]);
    dzetadxi_3 = dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_n  + zb[p_n ], htheta_nw + zb[p_nw]);

    dzetadeta_0 = dcdy_scv(htheta_0 + zb[p_0], htheta_s + zb[p_s], htheta_w  + zb[p_w ], htheta_sw + zb[p_sw]);
    dzetadeta_1 = dcdy_scv(htheta_0 + zb[p_0], htheta_s + zb[p_s], htheta_e  + zb[p_e ], htheta_se + zb[p_se]);
    dzetadeta_2 = dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_e  + zb[p_e ]);
    dzetadeta_3 = dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_nw + zb[p_nw], htheta_w  + zb[p_w ]);
    //
    // theta * g * (-dx_deta dzeta/dxi + dx_dxi dzeta/deta) * Delta h
    //
    // sub control volume 0 ============================================
    fac = theta * 0.25 * g;
    // Contribution to Delta h
    // scv_0
    add_value(values, col_0 , fac * (-dx_deta_0 * dzetadxi_0 * 9./16. + dx_dxi_0 * dzetadeta_0 * 9./16.));
    add_value(values, col_w , fac * (-dx_deta_0 * dzetadxi_0 * 3./16. + dx_dxi_0 * dzetadeta_0 * 3./16.));
    add_value(values, col_sw, fac * (-dx_deta_0 * dzetadxi_0 * 1./16. + dx_dxi_0 * dzetadeta_0 * 1./16.));
    add_value(values, col_s , fac * (-dx_deta_0 * dzetadxi_0 * 3./16. + dx_dxi_0 * dzetadeta_0 * 3./16.));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (-dx_deta_1 * dzetadxi_1 * 9./16. + dx_dxi_1 * dzetadeta_1) * 9./16.);
    add_value(values, col_s , fac * (-dx_deta_1 * dzetadxi_1 * 3./16. + dx_dxi_1 * dzetadeta_1) * 3./16.);
    add_value(values, col_se, fac * (-dx_deta_1 * dzetadxi_1 * 1./16. + dx_dxi_1 * dzetadeta_1) * 1./16.);
    add_value(values, col_e , fac * (-dx_deta_1 * dzetadxi_1 * 3./16. + dx_dxi_1 * dzetadeta_1) * 3./16.);
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (-dx_deta_2 * dzetadxi_2 * 9./16. + dx_dxi_2 * dzetadeta_2 * 9./16.));
    add_value(values, col_e , fac * (-dx_deta_2 * dzetadxi_2 * 3./16. + dx_dxi_2 * dzetadeta_2 * 3./16.));
    add_value(values, col_ne, fac * (-dx_deta_2 * dzetadxi_2 * 1./16. + dx_dxi_2 * dzetadeta_2 * 1./16.));
    add_value(values, col_n , fac * (-dx_deta_2 * dzetadxi_2 * 3./16. + dx_dxi_2 * dzetadeta_2 * 3./16.));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (-dx_deta_3 * dzetadxi_3 * 9./16. + dx_dxi_3 * dzetadeta_3 * 9./16.));
    add_value(values, col_n , fac * (-dx_deta_3 * dzetadxi_3 * 3./16. + dx_dxi_3 * dzetadeta_3 * 3./16.));
    add_value(values, col_nw, fac * (-dx_deta_3 * dzetadxi_3 * 1./16. + dx_dxi_3 * dzetadeta_3 * 1./16.));
    add_value(values, col_w , fac * (-dx_deta_3 * dzetadxi_3 * 3./16. + dx_dxi_3 * dzetadeta_3 * 3./16.));
    //
    // theta * g * h * d(Delta zeta)/dy
    //
    fac = theta * 0.25 * g * 0.25;
    // sub control volume 0 ============================================
    // scv_0
    add_value(values, col_0 , fac * ( 3. * -dx_deta_0 * depth_0 + 3. * dx_dxi_0 * depth_0));
    add_value(values, col_w , fac * (-3. * -dx_deta_0 * depth_0 + 1. * dx_dxi_0 * depth_0));
    add_value(values, col_sw, fac * (-1. * -dx_deta_0 * depth_0 - 1. * dx_dxi_0 * depth_0));
    add_value(values, col_s , fac * ( 1. * -dx_deta_0 * depth_0 - 3. * dx_dxi_0 * depth_0));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (-3. * -dx_deta_1 * depth_0 + 3. * dx_dxi_1 * depth_1));
    add_value(values, col_s , fac * (-1. * -dx_deta_1 * depth_0 - 3. * dx_dxi_1 * depth_1));
    add_value(values, col_se, fac * ( 1. * -dx_deta_1 * depth_0 - 1. * dx_dxi_1 * depth_1));
    add_value(values, col_e , fac * ( 3. * -dx_deta_1 * depth_0 + 1. * dx_dxi_1 * depth_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (-3. * -dx_deta_2 * depth_0 - 3. * dx_dxi_2 * depth_2));
    add_value(values, col_e , fac * ( 3. * -dx_deta_2 * depth_0 - 1. * dx_dxi_2 * depth_2));
    add_value(values, col_ne, fac * ( 1. * -dx_deta_2 * depth_0 + 1. * dx_dxi_2 * depth_2));
    add_value(values, col_n , fac * (-1. * -dx_deta_2 * depth_0 + 3. * dx_dxi_2 * depth_2));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * ( 3. * -dx_deta_3 * depth_0 - 3. * dx_dxi_3 * depth_3));
    add_value(values, col_n , fac * ( 1. * -dx_deta_3 * depth_0 + 3. * dx_dxi_3 * depth_3));
    add_value(values, col_nw, fac * (-1. * -dx_deta_3 * depth_0 + 1. * dx_dxi_3 * depth_3));
    add_value(values, col_w , fac * (-3. * -dx_deta_3 * depth_0 - 1. * dx_dxi_3 * depth_3));
    // 
    // RHS r-momentum equation
    //
    rhs[row + 2] += - 0.25 * g * (
        depth_0 * (-dx_deta_0 * dzetadxi_0 + dx_dxi_0 * dzetadeta_0) + 
        depth_1 * (-dx_deta_1 * dzetadxi_1 + dx_dxi_1 * dzetadeta_1) + 
        depth_2 * (-dx_deta_2 * dzetadxi_2 + dx_dxi_2 * dzetadeta_2) + 
        depth_3 * (-dx_deta_3 * dzetadxi_3 + dx_dxi_3 * dzetadeta_3)
        );
    return 0;
}
        
inline void add_value(double * values, int col, double data){ 
    values[col] += data; 
}
