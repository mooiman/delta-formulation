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
#include "matrix_assembly_psi_interior.h"

int reg_interior_matrix_psi(double* values, size_t row, size_t c_eq,
    double c_psi, struct _grid_metric metric)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    std:: vector<double>& x = metric.x;
    std:: vector<double>& y = metric.y;

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

    double scv_area_0 = 0.25;  // computational space
    double scv_area_1 = 0.25;  // computational space
    double scv_area_2 = 0.25;  // computational space
    double scv_area_3 = 0.25;  // computational space
    //------------------------------------------------------------------------
    // 
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

    
    // Diffusion part
    double n_xi = 0.0;
    double n_eta = 0.0;

    double c_psi_org = c_psi;
    c_psi = -c_psi;
    double fac = c_psi * 0.5;
    //
    // scv_0 face_0
    n_xi = -1.0;
    n_eta = 0.0;
    fac = c_psi * 0.5;
    add_value(values, col_0 , fac * n_xi *  3./4.);
    add_value(values, col_w , fac * n_xi * -3./4.);
    add_value(values, col_s , fac * n_xi *  1./4.);
    add_value(values, col_sw, fac * n_xi * -1./4.);
    //
    // scv_0 face_1
    n_xi = 0.0;
    n_eta = -1.0;
    fac = c_psi * 0.5;
    add_value(values, col_0 , fac * n_eta *  3./4.);
    add_value(values, col_s , fac * n_eta * -3./4.);
    add_value(values, col_w , fac * n_eta *  1./4.);
    add_value(values, col_sw, fac * n_eta * -1./4.);
    //
    // scv 1 face_2
    n_xi = 0.0;
    n_eta = -1.0;
    fac = c_psi * 0.5;
    add_value(values, col_0 , fac * n_eta *  3./4.);
    add_value(values, col_s , fac * n_eta * -3./4.);
    add_value(values, col_e , fac * n_eta *  1./4.);
    add_value(values, col_se, fac * n_eta * -1./4.);
    //
    // scv 1 face_3
    n_xi = 1.0;
    n_eta = 0.0;
    fac = c_psi * 0.5;
    add_value(values, col_e , fac * n_xi *  3./4.);
    add_value(values, col_0 , fac * n_xi * -3./4.);
    add_value(values, col_se, fac * n_xi *  1./4.);
    add_value(values, col_s , fac * n_xi * -1./4.);
    //
    // scv 1 face_4
    n_xi = 1.0;
    n_eta = 0.0;
    fac = c_psi * 0.5;
    add_value(values, col_e , fac * n_xi *  3./4.);
    add_value(values, col_0 , fac * n_xi * -3./4.);
    add_value(values, col_ne, fac * n_xi *  1./4.);
    add_value(values, col_n , fac * n_xi * -1./4.);
    //
    // scv 2 face_5
    n_xi = 0.0;
    n_eta = 1.0;
    fac = c_psi * 0.5;
    add_value(values, col_n , fac * n_eta *  3./4.);
    add_value(values, col_0 , fac * n_eta * -3./4.);
    add_value(values, col_ne, fac * n_eta *  1./4.);
    add_value(values, col_e , fac * n_eta * -1./4.);
    //
    // scv 3 face_6
    n_xi = 0.0;
    n_eta = 1.0;
    fac = c_psi * 0.5;
    add_value(values, col_n , fac * n_eta *  3./4.);
    add_value(values, col_0 , fac * n_eta * -3./4.);
    add_value(values, col_nw, fac * n_eta *  1./4.);
    add_value(values, col_w , fac * n_eta * -1./4.);
    //
    //scv 3 face_7
    n_xi = -1.0;
    n_eta = 0.0;
    fac = c_psi * 0.5;
    add_value(values, col_0 , fac * n_xi *  3./4.);
    add_value(values, col_w , fac * n_xi * -3./4.);
    add_value(values, col_n , fac * n_xi *  1./4.);
    add_value(values, col_nw, fac * n_xi * -1./4.);

    c_psi_org = c_psi;
    return 0;
}
//------------------------------------------------------------------------------
int reg_interior_rhs_psi( size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& h, std::vector<double>& q, std::vector<double>& r,
    double c_psi, double g, struct _grid_metric &  metric)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;
    std::vector<double> Err_psi(nxny, 0.0);

    std::vector<double> d2h_dxi2(nxny, 0.0);
    std::vector<double> d2h_dxideta(nxny, 0.0);
    std::vector<double> d2h_deta2(nxny, 0.0);
    std::vector<double> d2q_dxi2(nxny, 0.0);
    std::vector<double> d2q_dxideta(nxny, 0.0);
    std::vector<double> d2q_deta2(nxny, 0.0);
    std::vector<double> d2r_dxi2(nxny, 0.0);
    std::vector<double> d2r_dxideta(nxny, 0.0);
    std::vector<double> d2r_deta2(nxny, 0.0);
    std::vector<double> d2s_dxi2(nxny, 0.0);
    std::vector<double> d2s_dxideta(nxny, 0.0);
    std::vector<double> d2s_deta2(nxny, 0.0);
    std::vector<double> s(nxny, 0.0);
    std::vector<size_t> p(9);

    if ( row != c_eq/(9) ) { std::cerr << "Jan Mooiman" << std::endl; }
    p[0] = row - ny - 1;
    p[1] = row - ny ;
    p[2] = row - ny + 1;
    p[3] = row - 1;
    p[4] = row;
    p[5] = row + 1;
    p[6] = row + ny - 1;
    p[7] = row + ny;
    p[8] = row + ny + 1;

    //
    // Error based on potential energy
    //
    //   \sqrt{g \widehat{h}}
    // 

    double f1 = F1(h, p, metric);
    double f2 = F2(h, p, metric);
    double f3 = F3(h, p, metric);
    f2 = 0.0;
    rhs[row] = c_psi * std::sqrt( g/h[row] ) * (
        1.0/16.0 * f1 + 1.0/8.0 * f2 + 1.0/16.0 * f3
        );
    return 0;
}

inline double F1(std::vector<double> & u, std::vector<size_t>& p, 
    struct _grid_metric & metric)
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dx_deta = 1.0;
    double dxi_dx = 1.0;
    double deta_dx = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid
    double d2xi_dx2 = 0.0;
    double d2eta_dx2 = 0.0;

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dx_dxi * dx_dxi * (
          dxi_dx * dxi_dx * d2u_dxi2 
//        + 2.0 * dxi_dx * deta_dx * d2u_dxideta 
//        + deta_dx * deta_dx * d2u_deta2 
//        + d2xi_dx2 * du_dxi 
//        + d2eta_dx2 * du_deta
        );

    return retval;
}
inline double F2(std::vector<double> & u, std::vector<size_t>& p,
    struct _grid_metric & metric)
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dy_deta = 1.0;
    double dxi_dx = 1.0;
    double dxi_dy = 1.0;
    double deta_dx = 1.0;
    double deta_dy = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dx_dxi * dy_deta * (
        dxi_dx * dxi_dy * d2u_dxi2 
        + (dxi_dx * deta_dy + deta_dx * dxi_dy) * d2u_dxideta 
        + deta_dx * deta_dy * d2u_deta2
        + d2xi_dxdy * du_dxi
        + d2eta_dxdy * du_deta
        );

    return retval;
}
inline double F3(std::vector<double> & u, std::vector<size_t>& p, 
    struct _grid_metric & metric)
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dy_deta = 1.0;
    double dxi_dx = 1.0;
    double dxi_dy = 1.0;
    double deta_dx = 1.0;
    double deta_dy = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dy_deta * dy_deta * (
//         dxi_dy * dxi_dy * d2u_dxi2 
//       + 2.0 * dxi_dy * deta_dy * d2u_dxideta 
        + deta_dy * deta_dy * d2u_deta2 
//        + d2xi_dy2 * du_dxi 
//        + d2eta_dy2 * du_deta
        );

    return retval;
}
inline double d2udxi2(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    double dxi = 1.0;
    double deta = 1.0;
    double retval = 0.0;

    retval =  1./8. * u[p[0]] +
              6./8. * u[p[1]] +
              1./8. * u[p[2]] +
             -6./8. * u[p[3]] +
            -12./8. * u[p[4]] +
             -6./8. * u[p[5]] +
              1./8. * u[p[6]] +
              6./8. * u[p[7]] +
              1./8. * u[p[8]];
    return retval;
}
inline double d2udxideta(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    double dxi = 1.0;
    double deta = 1.0;
    double retval = 0.0;

    retval = 1.0 * u[p[0]] +
             0.0 * u[p[1]] +
            -1.0 * u[p[2]] +
             0.0 * u[p[3]] +
             0.0 * u[p[4]] +
             0.0 * u[p[5]] +
            -1.0 * u[p[6]] +
             0.0 * u[p[7]] +
             1.0 * u[p[8]];

    return 0.0;
    //return retval/(4.0 * dxi * deta);
}
inline double d2udeta2(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    double retval;

    retval =  1./8. * u[p[0]] +
             -6./8. * u[p[1]] +
              1./8. * u[p[2]] +
              6./8. * u[p[3]] +
            -12./8. * u[p[4]] +
              6./8. * u[p[5]] +
              1./8. * u[p[6]] +
             -6./8. * u[p[7]] +
              1./8. * u[p[8]];
    return retval;
    //return d2udxi2(u, p);
}

//------------------------------------------------------------------------------        
inline void add_value(double * values, size_t col, double data){ 
    values[col] += data; 
}
