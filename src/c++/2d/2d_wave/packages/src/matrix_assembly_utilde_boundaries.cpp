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
#include "matrix_assembly_utilde_boundaries.h"

int reg_boundary_north_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& u_giv, double psi_11, double psi_22, 
    double theta, struct _grid_metric & metric)
{
    size_t ny = metric.ny;

    size_t p_5 = c_eq/(9);  // node number of boundary point, ie north point of molecule
    size_t p_4 = p_5 - 1;
    size_t p_3 = p_5 - 2;
    size_t p_0 = p_5 - ny - 2;
    size_t p_1 = p_5 - ny - 1;
    size_t p_2 = p_5 - ny;
    size_t p_6 = p_5 + ny - 2;
    size_t p_7 = p_5 + ny - 1;
    size_t p_8 = p_5 + ny;

    size_t col_b  = c_eq + 5 * 1;
    size_t col_s  = c_eq + 4 * 1;
    size_t col_ss = c_eq + 3 * 1;

    // Contribution Delta h
    values[c_eq     ] = 0.0;
    values[c_eq +  1] = 0.0;
    values[c_eq +  2] = 0.0;
    values[c_eq +  3] = 0.0;
    values[c_eq +  4] = -1.0;
    values[c_eq +  5] = 1.0;
    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;
    rhs[row    ] = 0.0;

    return 0;
}
int reg_boundary_east_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& u_giv, double psi_11, double psi_22, 
    double theta, struct _grid_metric & metric)
{
    size_t p_0 = c_eq/(9);  // p_e

    values[c_eq     ] = 0.0;
    values[c_eq +  1] = 0.0;
    values[c_eq +  2] = 0.0;
    values[c_eq +  3] = 0.0;
    values[c_eq +  4] = -1.0;
    values[c_eq +  5] = 0.0;
    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 1.0;
    values[c_eq +  8] = 0.0;
    rhs[row    ] = 0.0;

    return 0;
}
int reg_boundary_south_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& u_giv, double psi_11, double psi_22, 
    double theta, struct _grid_metric & metric)
{
    size_t p_0 = c_eq/(9);  // p_s

    values[c_eq     ] = 0.0;
    values[c_eq +  1] = 0.0;
    values[c_eq +  2] = 0.0;
    values[c_eq +  3] = 1.0;
    values[c_eq +  4] = -1.0;
    values[c_eq +  5] = 0.0;
    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;
    rhs[row    ] = 0.0;

    return 0;
}
int reg_boundary_west_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& u_giv, double psi_11, double psi_22, 
    double theta, struct _grid_metric & metric)
{
    size_t p_0 = c_eq/(9);  // p_w

    values[c_eq     ] = 0.0;
    values[c_eq +  1] = -1.0;
    values[c_eq +  2] = 0.0;
    values[c_eq +  3] = 0.0;
    values[c_eq +  4] = 1.0;
    values[c_eq +  5] = 0.0;
    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;
    rhs[row    ] = 0.0;

    return 0;
}
