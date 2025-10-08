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

int interior_time(double* values, int row, int c_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double theta, 
    int nx, int ny,
    std::vector<double>& Tn,
    std::vector<double>& Tp, 
    double dx, double dy, double dxdy, std::vector<double>& mass)
{
    int p_0 = c_eq/(9);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (std::fmod(p_0, ny) == 0) { return 1; }  // south boundary
    if (std::fmod(p_0 + 1, ny) == 0) { return 2; }  // north boundary

    memset(&values[c_eq], 0, 9 * sizeof(double));  // set all coefficients for one row of Delta T to zero

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

    int col_sw = c_eq;
    int col_w  = c_eq + 1;
    int col_nw = c_eq + 2;
    int col_s  = c_eq + 3;
    int col_0  = c_eq + 4;
    int col_n  = c_eq + 5;
    int col_se = c_eq + 6;
    int col_e  = c_eq + 7;
    int col_ne = c_eq + 8;

    //------------------------------------------------------------------------
    // heat-equation
    // 
    set_value(values, col_sw, dtinv * dxdy * mass[0] * mass[0]);
    set_value(values, col_w , dtinv * dxdy * mass[1] * mass[0]);
    set_value(values, col_nw, dtinv * dxdy * mass[2] * mass[0]);

    set_value(values, col_s  ,dtinv * dxdy * mass[0] * mass[1]);
    set_value(values, col_0  ,dtinv * dxdy * mass[1] * mass[1]);
    set_value(values, col_n  ,dtinv * dxdy * mass[2] * mass[1]);

    set_value(values, col_se, dtinv * dxdy * mass[0] * mass[2]);
    set_value(values, col_e , dtinv * dxdy * mass[1] * mass[2]);
    set_value(values, col_ne, dtinv * dxdy * mass[2] * mass[2]);

    rhs[row] = -(
        dtinv * dxdy * mass[0] * mass[0] * (Tp[p_sw] - Tn[p_sw]) +
        dtinv * dxdy * mass[1] * mass[0] * (Tp[p_s ] - Tn[p_s ]) +
        dtinv * dxdy * mass[2] * mass[0] * (Tp[p_se] - Tn[p_se]) +
        //
        dtinv * dxdy * mass[0] * mass[1] * (Tp[p_w ] - Tn[p_w]) +
        dtinv * dxdy * mass[1] * mass[1] * (Tp[p_0 ] - Tn[p_0]) +
        dtinv * dxdy * mass[2] * mass[1] * (Tp[p_e]  - Tn[p_e]) +
        //
        dtinv * dxdy * mass[0] * mass[2] * (Tp[p_nw] - Tn[p_nw]) +
        dtinv * dxdy * mass[1] * mass[2] * (Tp[p_n ] - Tn[p_n ]) +
        dtinv * dxdy * mass[2] * mass[2] * (Tp[p_ne] - Tn[p_ne])
        );

    return 0;
}
        

inline void set_value(double * values, int col, double data){ 
    values[col] += data; 
}
