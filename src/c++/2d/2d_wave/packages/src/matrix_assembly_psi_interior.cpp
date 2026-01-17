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
    double c_psi, std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny)
{
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

    return 0;
}
        
inline void add_value(double * values, size_t col, double data){ 
    values[col] += data; 
}
