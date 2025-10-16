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
    double & dtinv, double & dxinv, double & theta, double & g, bool do_convection, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double dx, double dy, double dxdy, std::vector<double>& mass)
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

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

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

    //------------------------------------------------------------------------
    // c-equation
    // 
    rhs[row] = 0.0;
    // scv_0
    double scv_area = 0.25 * dxdy;
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_w , dtinv * scv_area * 3./16.);
    add_value(values, col_s , dtinv * scv_area * 3./16.);
    add_value(values, col_sw, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_w] - hn[p_w], hp[p_s] - hn[p_s], hp[p_sw] - hn[p_sw]);
    //
    // scv 1
    scv_area = 0.25 * dxdy; 
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_s , dtinv * scv_area * 3./16.);
    add_value(values, col_e , dtinv * scv_area * 3./16.);
    add_value(values, col_se, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_s] - hn[p_s], hp[p_e] - hn[p_e], hp[p_se] - hn[p_se]);
    //
    // scv 2
    scv_area = 0.25 * dxdy;
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_e , dtinv * scv_area * 3./16.);
    add_value(values, col_n , dtinv * scv_area * 3./16.);
    add_value(values, col_ne, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_e] - hn[p_e], hp[p_n] - hn[p_n], hp[p_ne] - hn[p_ne]);
    //
    //scv 3
    scv_area = 0.25 * dxdy;  
    add_value(values, col_0 , dtinv * scv_area * 9./16.);
    add_value(values, col_n , dtinv * scv_area * 3./16.);
    add_value(values, col_w , dtinv * scv_area * 3./16.);
    add_value(values, col_nw, dtinv * scv_area * 1./16.);
    rhs[row] += -dtinv * scv_area * c_scv(hp[p_0] - hn[p_0], hp[p_n] - hn[p_n], hp[p_w] - hn[p_w], hp[p_nw] - hn[p_nw]);
    //
    // Mass flux
    // 
    // d(q)/dxi
    // 
    // sub control volume 0 ============================================
    double dy_dxi = 0.5 * (y[p_0] - y[p_w]) + 0.5 * (y[p_s] - y[p_sw]);
    double dx_dxi = 0.5 * (x[p_0] - x[p_w]) + 0.5 * (x[p_s] - x[p_sw]);
    double dy_deta = 0.5 * (y[p_0] - y[p_s]) + 0.5 * (y[p_w] - y[p_sw]);
    double dx_deta = 0.5 * (x[p_0] - x[p_s]) + 0.5 * (x[p_w] - x[p_sw]);
    // 
    // scv_0 face_0
    double n_xi = -1.0;
    double n_eta = 0.0;

    double dl_y = 0.5 * dy_deta * n_xi;  // length in y-direction
    double dl_x = 0.5 * dx_deta * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_w  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_sw + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_s  + 1, theta * dl_y * 0.125 * 1.);

    add_value(values, col_0  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_w  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_sw + 2, -theta * dl_x * 0.125 * 1.);
    add_value(values, col_s  + 2, -theta * dl_x * 0.125 * 1.);

    // rhs
    double fluxx_0 = (scvf_xi(qtheta_0, qtheta_w, qtheta_s, qtheta_sw)) * dl_y 
                   - (scvf_xi(rtheta_0, rtheta_w, rtheta_s, rtheta_sw)) * dl_x;
    //
    // scv_0 face_1
    // No contribution to Delta r and rhs

    // sub control volume 1 ============================================
    dy_dxi = 0.5 * (y[p_e] - y[p_0]) + 0.5 * (y[p_se] - y[p_s]);
    dx_dxi = 0.5 * (x[p_e] - x[p_0]) + 0.5 * (x[p_se] - x[p_s]);
    dy_deta = 0.5 * (y[p_0] - y[p_s]) + 0.5 * (y[p_e] - y[p_se]);
    dx_deta = 0.5 * (x[p_0] - x[p_s]) + 0.5 * (x[p_e] - x[p_se]);
    //
    // scv_1 face_2
    // No contribution to Delta r and rhs
    // scv_1 face_3
    n_xi = 1.0;
    n_eta = 0.0;

    dl_y = 0.5 * dy_deta * n_xi;  // length in y-direction
    dl_x = 0.5 * dx_deta * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_s  + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_se + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_e  + 1, theta * dl_y * 0.125 * 3.);

    add_value(values, col_0  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_s  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_se + 2, -theta * dl_x * 0.125 * 1.);
    add_value(values, col_e  + 2, -theta * dl_x * 0.125 * 1.);

    // rhs
    double fluxx_3 = (scvf_xi(qtheta_e, qtheta_0, qtheta_se, qtheta_s)) * dl_y 
                   - (scvf_xi(rtheta_e, rtheta_0, rtheta_se, rtheta_s)) * dl_x;

                                
    // sub control volume 2 ============================================
    dy_dxi = 0.5 * (y[p_e] - y[p_0]) + 0.5 * (y[p_ne] - y[p_n]);
    dx_dxi = 0.5 * (x[p_e] - x[p_0]) + 0.5 * (x[p_ne] - x[p_n]);
    dy_deta = 0.5 * (y[p_n] - y[p_0]) + 0.5 * (y[p_ne] - y[p_e]);
    dx_deta = 0.5 * (x[p_n] - x[p_0]) + 0.5 * (x[p_ne] - x[p_e]);
    // scv_2 face_4
    n_xi =  1.0;
    n_eta = 0.0;

    dl_y = 0.5 * dy_deta * n_xi;  // length in y-direction
    dl_x = 0.5 * dx_deta * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_e  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_ne + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_n  + 1, theta * dl_y * 0.125 * 1.);

    add_value(values, col_0  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_s  + 2, -theta * dl_x * 0.125 * 3.);
    add_value(values, col_se + 2, -theta * dl_x * 0.125 * 1.);
    add_value(values, col_e  + 2, -theta * dl_x * 0.125 * 1.);

    double fluxx_4 = (scvf_xi(qtheta_e, qtheta_0, qtheta_ne, qtheta_n)) * dl_y 
                   - (scvf_xi(rtheta_e, rtheta_0, rtheta_ne, rtheta_n)) * dl_x;
    // scv_2 face_5
    // No contribution to Delta r and rhs
                                
    // sub control volume 3 ============================================
    dy_dxi = 0.5 * (y[p_0] - y[p_w]) + 0.5 * (y[p_n] - y[p_nw]);
    dx_dxi = 0.5 * (x[p_0] - x[p_w]) + 0.5 * (x[p_n] - x[p_nw]);
    dy_deta = 0.5 * (y[p_n] - y[p_0]) + 0.5 * (y[p_nw] - y[p_w]);
    dx_deta = 0.5 * (x[p_n] - x[p_0]) + 0.5 * (x[p_nw] - x[p_w]);
    // scv_3 face_6
    // No contribution to Delta r and rhs
    // scv_3 face_7
    // 
    n_xi = -1.0;
    n_eta = 0.0;

    dl_y = 0.5 * dy_deta * n_xi;  // length in y-direction
    dl_x = 0.5 * dx_deta * n_xi;  // length in x-direction

    add_value(values, col_0  + 1, theta * dl_y * 0.125 * 3.);
    add_value(values, col_n  + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_nw + 1, theta * dl_y * 0.125 * 1.);
    add_value(values, col_w  + 1, theta * dl_y * 0.125 * 3.);

    add_value(values, col_0  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_n  + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_nw + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_w  + 2, theta * dl_x * 0.125 * 3.);

    double fluxx_7 = (scvf_xi(qtheta_0, qtheta_w, qtheta_nw, qtheta_n)) * dl_y 
                   - (scvf_xi(rtheta_0, rtheta_w, rtheta_nw, rtheta_n)) * dl_x;
    //
    // mass flux
    // 
    // d(r)/deta
    //
    //
    // sub control volume 0 ============================================
    dy_dxi = 0.5 * (y[p_0] - y[p_w]) + 0.5 * (y[p_s] - y[p_sw]);
    dx_dxi = 0.5 * (x[p_0] - x[p_w]) + 0.5 * (x[p_s] - x[p_sw]);
    dy_deta = 0.5 * (y[p_0] - y[p_s]) + 0.5 * (y[p_w] - y[p_sw]);
    dx_deta = 0.5 * (x[p_0] - x[p_s]) + 0.5 * (x[p_w] - x[p_sw]);
    //
    // scv_0 face_0
    // No contribution to Delta q and rhs
    // scv_0 face_1
    n_xi = 0.0;
    n_eta = -1.0;

    dl_y = 0.5 * dy_dxi * n_eta;  // length in y-direction
    dl_x = 0.5 * dx_dxi * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_w  + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_sw + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_s  + 1, -theta * dl_y * 0.125 * 3.);

    add_value(values, col_0  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_w  + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_sw + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_s  + 2, theta * dl_x * 0.125 * 3.);
    // rhs
    double fluxy_1 = - (scvf_eta(qtheta_0, qtheta_s, qtheta_sw, qtheta_w)) * dl_y 
                     + (scvf_eta(rtheta_0, rtheta_s, rtheta_sw, rtheta_w)) * dl_x;
    //                        
    //                        
    // sub control volume 1 ============================================
    // scv_1 face_2
    
    // sub control volume 1 ============================================
    dy_dxi = 0.5 * (y[p_e] - y[p_0]) + 0.5 * (y[p_se] - y[p_s]);
    dx_dxi = 0.5 * (x[p_e] - x[p_0]) + 0.5 * (x[p_se] - x[p_s]);
    dy_deta = 0.5 * (y[p_0] - y[p_s]) + 0.5 * (y[p_e] - y[p_se]);
    dx_deta = 0.5 * (x[p_0] - x[p_s]) + 0.5 * (x[p_e] - x[p_se]);
    //
    // scv_1 face_2
    n_xi = 0.0;
    n_eta = -1.0;

    dl_y = 0.5 * dy_dxi * n_eta;  // length in y-direction
    dl_x = 0.5 * dx_dxi * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_s  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_se + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_e  + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_0  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_s  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_se + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_e  + 2, theta * dl_x * 0.125 * 1.);
    // rhs
    double fluxy_2 = - (scvf_eta(qtheta_0, qtheta_s, qtheta_se, qtheta_e)) * dl_y 
                        + (scvf_eta(rtheta_0, rtheta_s, rtheta_se, rtheta_e)) * dl_x;
    // scv_1 face_3
    // No contribution to r-momentum equation
                            
    // sub control volume 2 ============================================
    // scv_1 face_4
    // No contribution to Delta q and rhs
    // scv_1 face_5
    dy_dxi = 0.5 * (y[p_e] - y[p_0]) + 0.5 * (y[p_ne] - y[p_n]);
    dx_dxi = 0.5 * (x[p_e] - x[p_0]) + 0.5 * (x[p_ne] - x[p_n]);
    dy_deta = 0.5 * (y[p_n] - y[p_0]) + 0.5 * (y[p_ne] - y[p_e]);
    dx_deta = 0.5 * (x[p_n] - x[p_0]) + 0.5 * (x[p_ne] - x[p_e]);

    n_xi =  0.0;
    n_eta = 1.0;

    dl_y = 0.5 * dy_dxi * n_eta;  // length in y-direction
    dl_x = 0.5 * dx_dxi * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_e  + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_ne + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_n  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_0  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_e  + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_ne + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_n  + 2, theta * dl_x * 0.125 * 3.);
    // rhs
    double fluxy_5 = - (scvf_eta(qtheta_n, qtheta_0, qtheta_ne, qtheta_e)) * dl_y 
                     + (scvf_eta(rtheta_n, rtheta_0, rtheta_ne, rtheta_e)) * dl_x;
                           
    // sub control volume 3 ============================================
    dy_dxi = 0.5 * (y[p_0] - y[p_w]) + 0.5 * (y[p_n] - y[p_nw]);
    dx_dxi = 0.5 * (x[p_0] - x[p_w]) + 0.5 * (x[p_n] - x[p_nw]);
    dy_deta = 0.5 * (y[p_n] - y[p_0]) + 0.5 * (y[p_nw] - y[p_w]);
    dx_deta = 0.5 * (x[p_n] - x[p_0]) + 0.5 * (x[p_nw] - x[p_w]);
    // scv_1 face_6
    n_xi =  0.0;
    n_eta = 1.0;

    dl_y = 0.5 * dy_dxi * n_eta;  // length in y-direction
    dl_x = 0.5 * dx_dxi * n_eta;  // length in x-direction

    add_value(values, col_0  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_n  + 1, -theta * dl_y * 0.125 * 3.);
    add_value(values, col_nw + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_w  + 1, -theta * dl_y * 0.125 * 1.);
    add_value(values, col_0  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_n  + 2, theta * dl_x * 0.125 * 3.);
    add_value(values, col_nw + 2, theta * dl_x * 0.125 * 1.);
    add_value(values, col_w  + 2, theta * dl_x * 0.125 * 1.);
    // rhs
    double fluxy_6 = - (scvf_eta(qtheta_n, qtheta_0, qtheta_nw, qtheta_w)) * dl_y 
                     + (scvf_eta(rtheta_n, rtheta_0, rtheta_nw, rtheta_w)) * dl_x;
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


    //--------------------------------------------------------------------------
    // q-equation
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
    add_value(values, col_sw + 1, dtinv * dxdy * mass[0] * mass[0]);  // Delta h
    add_value(values, col_w  + 1, dtinv * dxdy * mass[1] * mass[0]);  // Delta q
    add_value(values, col_nw + 1, dtinv * dxdy * mass[2] * mass[0]);  // Delta r

    add_value(values, col_s  + 1 ,dtinv * dxdy * mass[0] * mass[1]);   // Delta h
    add_value(values, col_0  + 1 ,dtinv * dxdy * mass[1] * mass[1]);   // Delta q
    add_value(values, col_n  + 1 ,dtinv * dxdy * mass[2] * mass[1]);   // Delta r

    add_value(values, col_se + 1, dtinv * dxdy * mass[0] * mass[2]);   // Delta h
    add_value(values, col_e  + 1, dtinv * dxdy * mass[1] * mass[2]);   // Delta q
    add_value(values, col_ne + 1, dtinv * dxdy * mass[2] * mass[2]);   // Delta r

    rhs[row + 1] = -(
        dtinv * dxdy * mass[0] * mass[0] * (qp[p_sw] - qn[p_sw]) +
        dtinv * dxdy * mass[1] * mass[0] * (qp[p_s ] - qn[p_s ]) +
        dtinv * dxdy * mass[2] * mass[0] * (qp[p_se] - qn[p_se]) +
        //
        dtinv * dxdy * mass[0] * mass[1] * (qp[p_w] - qn[p_w]) +
        dtinv * dxdy * mass[1] * mass[1] * (qp[p_0] - qn[p_0]) +
        dtinv * dxdy * mass[2] * mass[1] * (qp[p_e] - qn[p_e]) +
        //
        dtinv * dxdy * mass[0] * mass[2] * (qp[p_nw] - qn[p_nw]) +
        dtinv * dxdy * mass[1] * mass[2] * (qp[p_n ] - qn[p_n ]) +
        dtinv * dxdy * mass[2] * mass[2] * (qp[p_ne] - qn[p_ne])
    );
    //
    // pressure term: gh d(zeta)/dx
    //
    double depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
    double depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
    double depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
    double depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

    scv_area = 0.25 * dxdy;
    double dzetadx_0 = 1.0 / dx * dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_s  + zb[p_s ], htheta_sw + zb[p_sw]);
    double dzetadx_1 = 1.0 / dx * dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_se + zb[p_se], htheta_s  + zb[p_s ]);
    double dzetadx_2 = 1.0 / dx * dcdx_scv(htheta_e + zb[p_e], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_n  + zb[p_n ]);
    double dzetadx_3 = 1.0 / dx * dcdx_scv(htheta_0 + zb[p_0], htheta_w + zb[p_w], htheta_n  + zb[p_n ], htheta_nw + zb[p_nw]);
    //
    // theta * g * dzeta/dx * Delta h
    // Contribution to Delta h
    double fac = theta * scv_area * g * 0.0625;
    // sub control volume 0 ============================================
    // scv_0
    add_value(values, col_0 , fac * (9. * dzetadx_0));   
    add_value(values, col_w , fac * (3. * dzetadx_0));
    add_value(values, col_sw, fac * (1. * dzetadx_0));
    add_value(values, col_s , fac * (3. * dzetadx_0)); 
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (9. * dzetadx_1));   
    add_value(values, col_s , fac * (3. * dzetadx_1)); 
    add_value(values, col_se, fac * (1. * dzetadx_1));
    add_value(values, col_e , fac * (3. * dzetadx_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (9. * dzetadx_2));   
    add_value(values, col_e , fac * (3. * dzetadx_2));
    add_value(values, col_ne, fac * (1. * dzetadx_2));
    add_value(values, col_n , fac * (3. * dzetadx_2));   
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (9. * dzetadx_3));   
    add_value(values, col_n , fac * (3. * dzetadx_3));   
    add_value(values, col_nw, fac * (1. * dzetadx_3));
    add_value(values, col_w , fac * (3. * dzetadx_3)); 
    //
    // theta * g * h * d(Delta zeta)/dx 
    //
    fac = theta * scv_area * g * 0.25 / dx;
    // sub control volume 0 ============================================
    // scv_0
    add_value(values, col_0 , fac * (+3. * depth_0));
    add_value(values, col_w , fac * (-3. * depth_0));
    add_value(values, col_sw, fac * (-1. * depth_0));
    add_value(values, col_s , fac * (+1. * depth_0));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (-3. * depth_1));
    add_value(values, col_s , fac * (-1. * depth_1));
    add_value(values, col_se, fac * (+1. * depth_1));
    add_value(values, col_e , fac * (+3. * depth_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (-3. * depth_2));
    add_value(values, col_e , fac * (+3. * depth_2));
    add_value(values, col_ne, fac * (+1. * depth_2));
    add_value(values, col_n , fac * (-1. * depth_2));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (+3. * depth_3));
    add_value(values, col_n , fac * (+1. * depth_3));
    add_value(values, col_nw, fac * (-1. * depth_3));
    add_value(values, col_w , fac * (-3. * depth_3));
    //
    // RHS q-momentum equation
    //
    rhs[row + 1] += -scv_area * g * (depth_0 * dzetadx_0 + depth_1 * dzetadx_1 + depth_2 * dzetadx_2 + depth_3 * dzetadx_3);;
    //
    //--------------------------------------------------------------------------
    // r-equation
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
    add_value(values, col_sw + 2, dtinv * dxdy * mass[0] * mass[0]);  // Delta h
    add_value(values, col_w  + 2, dtinv * dxdy * mass[1] * mass[0]);  // Delta q
    add_value(values, col_nw + 2, dtinv * dxdy * mass[2] * mass[0]);  // Delta r

    add_value(values, col_s  + 2, dtinv * dxdy * mass[0] * mass[1]);   // Delta h
    add_value(values, col_0  + 2, dtinv * dxdy * mass[1] * mass[1]);   // Delta q
    add_value(values, col_n  + 2, dtinv * dxdy * mass[2] * mass[1]);   // Delta r

    add_value(values, col_se + 2, dtinv * dxdy * mass[0] * mass[2]);   // Delta h
    add_value(values, col_e  + 2, dtinv * dxdy * mass[1] * mass[2]);   // Delta q
    add_value(values, col_ne + 2, dtinv * dxdy * mass[2] * mass[2]);   // Delta r

    rhs[row + 2] = -(
        dtinv * dxdy * mass[0] * mass[0] * (rp[p_sw] - rn[p_sw]) +
        dtinv * dxdy * mass[1] * mass[0] * (rp[p_s ] - rn[p_s ]) +
        dtinv * dxdy * mass[2] * mass[0] * (rp[p_se] - rn[p_se]) +
        //
        dtinv * dxdy * mass[0] * mass[1] * (rp[p_w] - rn[p_w]) +
        dtinv * dxdy * mass[1] * mass[1] * (rp[p_0] - rn[p_0]) +
        dtinv * dxdy * mass[2] * mass[1] * (rp[p_e] - rn[p_e]) +
        //
        dtinv * dxdy * mass[0] * mass[2] * (rp[p_nw] - rn[p_nw]) +
        dtinv * dxdy * mass[1] * mass[2] * (rp[p_n ] - rn[p_n ]) +
        dtinv * dxdy * mass[2] * mass[2] * (rp[p_ne] - rn[p_ne])
    );
    //
    // pressure term: gh d(zeta)/dy
    //
    depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
    depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
    depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
    depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

    scv_area = 0.25 * dxdy;
    double dzetady_0 = 1.0 / dy * dcdy_scv(htheta_0 + zb[p_0], htheta_s + zb[p_s], htheta_w  + zb[p_w ], htheta_sw + zb[p_sw]);
    double dzetady_1 = 1.0 / dy * dcdy_scv(htheta_0 + zb[p_0], htheta_s + zb[p_s], htheta_e  + zb[p_e ], htheta_se + zb[p_se]);
    double dzetady_2 = 1.0 / dy * dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_ne + zb[p_ne], htheta_e  + zb[p_e ]);
    double dzetady_3 = 1.0 / dy * dcdy_scv(htheta_n + zb[p_n], htheta_0 + zb[p_0], htheta_nw + zb[p_nw], htheta_w  + zb[p_w ]);

    // theta * dzeta/dy * Delta h
    // sub control volume 0 ============================================
    fac = theta * scv_area * g * 0.0625;
    // Contribution to Delta h
    // scv_0
    add_value(values, col_0 , fac * (9. * dzetady_0));
    add_value(values, col_w , fac * (3. * dzetady_0));
    add_value(values, col_sw, fac * (1. * dzetady_0));
    add_value(values, col_s , fac * (3. * dzetady_0));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (9. * dzetady_1));
    add_value(values, col_s , fac * (3. * dzetady_1));
    add_value(values, col_se, fac * (1. * dzetady_1));
    add_value(values, col_e , fac * (3. * dzetady_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (9. * dzetady_2));
    add_value(values, col_e , fac * (3. * dzetady_2));
    add_value(values, col_ne, fac * (1. * dzetady_2));
    add_value(values, col_n , fac * (3. * dzetady_2));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (9. * dzetady_3));
    add_value(values, col_n , fac * (3. * dzetady_3));
    add_value(values, col_nw, fac * (1. * dzetady_3));
    add_value(values, col_w , fac * (3. * dzetady_3));
    //
    // theta * h * d(Delta zeta)/dy
    fac = theta * scv_area * g * 0.25 / dy;
    // sub control volume 0 ============================================
    // scv_0
    add_value(values, col_0 , fac * (+ 3. * depth_0));
    add_value(values, col_w , fac * (+ 1. * depth_0));
    add_value(values, col_sw, fac * (- 1. * depth_0));
    add_value(values, col_s , fac * (- 3. * depth_0));
    // sub control volume 1 ============================================
    // scv_1
    add_value(values, col_0 , fac * (+ 3. * depth_1));
    add_value(values, col_s , fac * (- 3. * depth_1));
    add_value(values, col_se, fac * (- 1. * depth_1));
    add_value(values, col_e , fac * (+ 1. * depth_1));
    // sub control volume 2 ============================================
    // scv_2
    add_value(values, col_0 , fac * (- 3. * depth_2));
    add_value(values, col_e , fac * (- 1. * depth_2));
    add_value(values, col_ne, fac * (+ 1. * depth_2));
    add_value(values, col_n , fac * (+ 3. * depth_2));
    // sub control volume 3 ============================================
    // scv_3
    add_value(values, col_0 , fac * (- 3. * depth_3));
    add_value(values, col_n , fac * (+ 3. * depth_3));
    add_value(values, col_nw, fac * (+ 1. * depth_3));
    add_value(values, col_w , fac * (- 1. * depth_3));
    // 
    // RHS r-momentum equation
    //
    rhs[row + 2] += -scv_area * g * (depth_0 * dzetady_0 + depth_1 * dzetady_1 + depth_2 * dzetady_2 + depth_3 * dzetady_3);

    return 0;
}
        
inline void add_value(double * values, int col, double data){ 
    values[col] += data; 
}
