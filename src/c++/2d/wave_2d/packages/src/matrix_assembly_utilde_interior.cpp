//
// programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formuation and Modified Newton iteration 
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
//
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

#include "interpolations.h"
#include "matrix_assembly_utilde_interior.h"

//==============================================================================
// 
//  \int_\Omega utilde - \int_\Omega \nabla \dotp \Psi utilde 
//
//==============================================================================

int utilde_interior_matrix(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    double psi_11, double psi_22, struct _grid_metric & metric)
{
    size_t ny = metric.ny;
    std::vector<double> x = metric.x;
    std::vector<double> y = metric.y;

    size_t p_0 = c_eq/(9);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (p_0 % ny == 0) { return 1; }  // south boundary
    if ((p_0 + 1) % ny == 0) { return 2; }  // north boundary

    std::fill_n(&values[c_eq], 9, 0.0);  // set all coefficients for one row of Delta T to zero

    size_t p_sw = p_0 - ny - 1;
    size_t p_w  = p_0 - ny;
    size_t p_nw = p_0 - ny + 1;
    size_t p_s  = p_0 - 1; 
    size_t p_n  = p_0 + 1;
    size_t p_se = p_0 + ny - 1;
    size_t p_e  = p_0 + ny;
    size_t p_ne = p_0 + ny + 1;

    size_t col_sw = c_eq;
    size_t col_w  = c_eq + 1;
    size_t col_nw = c_eq + 2;
    size_t col_s  = c_eq + 3;
    size_t col_0  = c_eq + 4;
    size_t col_n  = c_eq + 5;
    size_t col_se = c_eq + 6;
    size_t col_e  = c_eq + 7;
    size_t col_ne = c_eq + 8;

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
    //==============================================================================
    // 
    //  \int_\Omega utilde 
    //
    //==============================================================================
    // scv_0
    add_value(values, col_0 , scv_area_0 * 9./16.);
    add_value(values, col_w , scv_area_0 * 3./16.);
    add_value(values, col_s , scv_area_0 * 3./16.);
    add_value(values, col_sw, scv_area_0 * 1./16.);
    //
    // scv 1
    add_value(values, col_0 , scv_area_1 * 9./16.);
    add_value(values, col_s , scv_area_1 * 3./16.);
    add_value(values, col_e , scv_area_1 * 3./16.);
    add_value(values, col_se, scv_area_1 * 1./16.);
    //
    // scv 2
    add_value(values, col_0 , scv_area_2 * 9./16.);
    add_value(values, col_e , scv_area_2 * 3./16.);
    add_value(values, col_n , scv_area_2 * 3./16.);
    add_value(values, col_ne, scv_area_2 * 1./16.);
    //
    //scv 3
    add_value(values, col_0 , scv_area_3 * 9./16.);
    add_value(values, col_n , scv_area_3 * 3./16.);
    add_value(values, col_w , scv_area_3 * 3./16.);
    add_value(values, col_nw, scv_area_3 * 1./16.);

    //==============================================================================
    // 
    //  \int_\Omega \nabla \dotp \nalbla \Psi utilde =
    //
    //==============================================================================

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
    // diffusion part
    // 

    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    double n_xi = -1.0;
    double n_eta =  0.0;

    double scvf_fac = 0.5 * dy_deta_f0 * n_xi ;
    add_value(values, col_0 , scvf_fac * ( + psi_11 *  dy_deta_f0 * 3./4. + psi_11 * -dy_dxi_f0 * 1./2.) );
    add_value(values, col_w , scvf_fac * ( - psi_11 *  dy_deta_f0 * 3./4. + psi_11 * -dy_dxi_f0 * 1./2.) );
    add_value(values, col_s , scvf_fac * ( + psi_11 *  dy_deta_f0 * 1./4. - psi_11 * -dy_dxi_f0 * 1./2.) );
    add_value(values, col_sw, scvf_fac * ( - psi_11 *  dy_deta_f0 * 1./4. - psi_11 * -dy_dxi_f0 * 1./2.) );

    scvf_fac = 0.5 * -dx_deta_f0 * n_xi ;
    add_value(values, col_0 , scvf_fac * ( + psi_22 * -dx_deta_f0 * 3./4. + psi_22 *  dx_dxi_f0 * 1./2.) );
    add_value(values, col_w , scvf_fac * ( - psi_22 * -dx_deta_f0 * 3./4. + psi_22 *  dx_dxi_f0 * 1./2.) );
    add_value(values, col_s , scvf_fac * ( + psi_22 * -dx_deta_f0 * 1./4. - psi_22 *  dx_dxi_f0 * 1./2.) );
    add_value(values, col_sw, scvf_fac * ( - psi_22 * -dx_deta_f0 * 1./4. - psi_22 *  dx_dxi_f0 * 1./2.) );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    scvf_fac = 0.5 * -dy_dxi_f1 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( + psi_11 *  dy_deta_f1 * 1./2. + psi_11 * -dy_dxi_f1 * 3./4. ) );
    add_value(values, col_w , scvf_fac * ( - psi_11 *  dy_deta_f1 * 1./2. + psi_11 * -dy_dxi_f1 * 1./4. ) );
    add_value(values, col_s , scvf_fac * ( + psi_11 *  dy_deta_f1 * 1./2. - psi_11 * -dy_dxi_f1 * 3./4. ) );
    add_value(values, col_sw, scvf_fac * ( - psi_11 *  dy_deta_f1 * 1./2. - psi_11 * -dy_dxi_f1 * 1./4. ) );

    scvf_fac = 0.5 * dx_dxi_f1 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( + psi_22 * -dx_deta_f1 * 1./2. + psi_22 *  dx_dxi_f1 * 3./4. ) );
    add_value(values, col_w , scvf_fac * ( + psi_22 * -dx_deta_f1 * 1./2. + psi_22 *  dx_dxi_f1 * 1./4. ) );
    add_value(values, col_s , scvf_fac * ( - psi_22 * -dx_deta_f1 * 1./2. - psi_22 *  dx_dxi_f1 * 3./4. ) );
    add_value(values, col_sw, scvf_fac * ( - psi_22 * -dx_deta_f1 * 1./2. - psi_22 *  dx_dxi_f1 * 1./4. ) );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    scvf_fac = 0.5 * -dy_dxi_f2 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_11 * dy_deta_f2 * 1./2. + psi_11 * -dy_dxi_f2 * 3./4. ) );
    add_value(values, col_s , scvf_fac * ( - psi_11 * dy_deta_f2 * 1./2. - psi_11 * -dy_dxi_f2 * 3./4. ) );
    add_value(values, col_e , scvf_fac * ( + psi_11 * dy_deta_f2 * 1./2. + psi_11 * -dy_dxi_f2 * 1./4. ) );
    add_value(values, col_se, scvf_fac * ( + psi_11 * dy_deta_f2 * 1./2. - psi_11 * -dy_dxi_f2 * 1./4. ) );

    scvf_fac = 0.5 * dx_dxi_f2 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_22 * -dx_deta_f2 * 1./2. + psi_22 * dx_dxi_f2 * 3./4. ) );
    add_value(values, col_s , scvf_fac * ( - psi_22 * -dx_deta_f2 * 1./2. - psi_22 * dx_dxi_f2 * 3./4. ) );
    add_value(values, col_e , scvf_fac * ( + psi_22 * -dx_deta_f2 * 1./2. + psi_22 * dx_dxi_f2 * 1./4. ) );
    add_value(values, col_se, scvf_fac * ( + psi_22 * -dx_deta_f2 * 1./2. - psi_22 * dx_dxi_f2 * 1./4. ) );
 
    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    scvf_fac = 0.5 * dy_deta_f3 * n_xi ;
    add_value(values, col_e , scvf_fac * ( + psi_11 * dy_deta_f3 * 3./4. + psi_11 * -dy_dxi_f3 * 1./2.) );
    add_value(values, col_0 , scvf_fac * ( - psi_11 * dy_deta_f3 * 3./4. + psi_11 * -dy_dxi_f3 * 1./2.) );
    add_value(values, col_se, scvf_fac * ( + psi_11 * dy_deta_f3 * 1./4. - psi_11 * -dy_dxi_f3 * 1./2.) );
    add_value(values, col_s , scvf_fac * ( - psi_11 * dy_deta_f3 * 1./4. - psi_11 * -dy_dxi_f3 * 1./2.) );

    scvf_fac = 0.5 * -dx_deta_f3 * n_xi ;
    add_value(values, col_e , scvf_fac * ( + psi_22 * -dx_deta_f3 * 3./4. + psi_22 * dx_dxi_f3 * 1./2.) );
    add_value(values, col_0 , scvf_fac * ( - psi_22 * -dx_deta_f3 * 3./4. + psi_22 * dx_dxi_f3 * 1./2.) );
    add_value(values, col_se, scvf_fac * ( + psi_22 * -dx_deta_f3 * 1./4. - psi_22 * dx_dxi_f3 * 1./2.) );
    add_value(values, col_s , scvf_fac * ( - psi_22 * -dx_deta_f3 * 1./4. - psi_22 * dx_dxi_f3 * 1./2.) );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    scvf_fac = 0.5 * dy_deta_f4 * n_xi ;
    add_value(values, col_e , scvf_fac * ( + psi_11 * 3./4. * dy_deta_f4 + psi_11 * -dy_dxi_f4 * 1./2.) );
    add_value(values, col_0 , scvf_fac * ( - psi_11 * 3./4. * dy_deta_f4 + psi_11 * -dy_dxi_f4 * 1./2.) );
    add_value(values, col_ne, scvf_fac * ( + psi_11 * 1./4. * dy_deta_f4 - psi_11 * -dy_dxi_f4 * 1./2.) );
    add_value(values, col_n , scvf_fac * ( - psi_11 * 1./4. * dy_deta_f4 - psi_11 * -dy_dxi_f4 * 1./2.) );

    scvf_fac = 0.5 * -dx_deta_f4 * n_xi ;
    add_value(values, col_e , scvf_fac * ( + psi_22 * 3./4. * -dx_deta_f4 + psi_22 * dx_dxi_f4 * 1./2) );
    add_value(values, col_0 , scvf_fac * ( - psi_22 * 3./4. * -dx_deta_f4 + psi_22 * dx_dxi_f4 * 1./2) );
    add_value(values, col_ne, scvf_fac * ( + psi_22 * 1./4. * -dx_deta_f4 - psi_22 * dx_dxi_f4 * 1./2) );
    add_value(values, col_n , scvf_fac * ( - psi_22 * 1./4. * -dx_deta_f4 - psi_22 * dx_dxi_f4 * 1./2) );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    scvf_fac = 0.5 * -dy_dxi_f5 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_11 * dy_deta_f5 * 1./2. - psi_11 * -dy_dxi_f5 * 3./4. ) );
    add_value(values, col_n , scvf_fac * ( - psi_11 * dy_deta_f5 * 1./2. + psi_11 * -dy_dxi_f5 * 3./4. ) );
    add_value(values, col_e , scvf_fac * ( + psi_11 * dy_deta_f5 * 1./2. - psi_11 * -dy_dxi_f5 * 1./4. ) );
    add_value(values, col_ne, scvf_fac * ( + psi_11 * dy_deta_f5 * 1./2. + psi_11 * -dy_dxi_f5 * 1./4. ) );

    scvf_fac = 0.5 * dx_dxi_f5 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_22 * -dx_deta_f5 * 1./2. - psi_22 * dx_dxi_f5 * 3./4. ) );
    add_value(values, col_n , scvf_fac * ( - psi_22 * -dx_deta_f5 * 1./2. + psi_22 * dx_dxi_f5 * 3./4. ) );
    add_value(values, col_e , scvf_fac * ( + psi_22 * -dx_deta_f5 * 1./2. - psi_22 * dx_dxi_f5 * 1./4. ) );
    add_value(values, col_ne, scvf_fac * ( + psi_22 * -dx_deta_f5 * 1./2. + psi_22 * dx_dxi_f5 * 1./4. ) );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    scvf_fac = 0.5 * -dy_dxi_f6 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_11 * dy_deta_f6 * 1./2. - psi_11 * -dy_dxi_f6 * 3./4. ) );
    add_value(values, col_n , scvf_fac * ( + psi_11 * dy_deta_f6 * 1./2. + psi_11 * -dy_dxi_f6 * 3./4. ) );
    add_value(values, col_w , scvf_fac * ( + psi_11 * dy_deta_f6 * 1./2. - psi_11 * -dy_dxi_f6 * 1./4. ) );
    add_value(values, col_nw, scvf_fac * ( - psi_11 * dy_deta_f6 * 1./2. + psi_11 * -dy_dxi_f6 * 1./4. ) );

    scvf_fac = 0.5 * dx_dxi_f6 * n_eta ;
    add_value(values, col_0 , scvf_fac * ( - psi_22 * -dx_deta_f6 * 1./2. - psi_22 * dx_dxi_f6 * 3./4. ) );
    add_value(values, col_n , scvf_fac * ( + psi_22 * -dx_deta_f6 * 1./2. + psi_22 * dx_dxi_f6 * 3./4. ) );
    add_value(values, col_w , scvf_fac * ( - psi_22 * -dx_deta_f6 * 1./2. - psi_22 * dx_dxi_f6 * 1./4. ) );
    add_value(values, col_nw, scvf_fac * ( + psi_22 * -dx_deta_f6 * 1./2. + psi_22 * dx_dxi_f6 * 1./4. ) );

    // scv_3 face_7
    n_xi =  -1.0;
    n_eta =  0.0;

    scvf_fac = 0.5 * dy_deta_f7 * n_xi ;
    add_value(values, col_0 , scvf_fac * ( + psi_11 * dy_deta_f7 * 3./4. - psi_11 * -dy_dxi_f7 * 1./2.) );
    add_value(values, col_w , scvf_fac * ( - psi_11 * dy_deta_f7 * 3./4. - psi_11 * -dy_dxi_f7 * 1./2.) );
    add_value(values, col_n , scvf_fac * ( + psi_11 * dy_deta_f7 * 1./4. + psi_11 * -dy_dxi_f7 * 1./2.) );
    add_value(values, col_nw, scvf_fac * ( - psi_11 * dy_deta_f7 * 1./4. + psi_11 * -dy_dxi_f7 * 1./2.) );

    scvf_fac = 0.5 * -dx_deta_f7 * n_xi ;
    add_value(values, col_0 , scvf_fac * ( + psi_22 * -dx_deta_f7 * 3./4 - psi_22 * dx_dxi_f7 * 1./2.) );
    add_value(values, col_w , scvf_fac * ( - psi_22 * -dx_deta_f7 * 3./4 - psi_22 * dx_dxi_f7 * 1./2.) );
    add_value(values, col_n , scvf_fac * ( + psi_22 * -dx_deta_f7 * 1./4 + psi_22 * dx_dxi_f7 * 1./2.) );
    add_value(values, col_nw, scvf_fac * ( - psi_22 * -dx_deta_f7 * 1./4 + psi_22 * dx_dxi_f7 * 1./2.) );

    rhs[row] = 0.0;

    return 0;
}
int utilde_interior_rhs(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& u_giv, struct _grid_metric & metric)
{
    size_t ny = metric.ny;
    std::vector<double> x = metric.x;
    std::vector<double> y = metric.y;

    size_t p_0 = c_eq/(9);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (p_0 % ny == 0) { return 1; }  // south boundary
    if ((p_0 + 1) % ny == 0) { return 2; }  // north boundary

    int refine = 2;
    // loop over finite volumes 
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

    rhs[row]  = scv_area_0 * c_scv(u_giv[p_0], u_giv[p_w], u_giv[p_s], u_giv[p_sw]);
    rhs[row] += scv_area_1 * c_scv(u_giv[p_0], u_giv[p_s], u_giv[p_e], u_giv[p_se]);
    rhs[row] += scv_area_2 * c_scv(u_giv[p_0], u_giv[p_e], u_giv[p_n], u_giv[p_ne]);
    rhs[row] += scv_area_3 * c_scv(u_giv[p_0], u_giv[p_n], u_giv[p_w], u_giv[p_nw]);

    return 0;
}
inline void add_value(double * values, size_t col, double data){
    values[col] += data; 
}
