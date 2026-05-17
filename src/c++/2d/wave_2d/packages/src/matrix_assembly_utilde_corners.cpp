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

#include "matrix_assembly_utilde_corners.h"

//
//corner nodes
//
int reg_corner_north_east_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, struct _grid_metric & metric)
{
    size_t ny = metric.ny;

    size_t p0 = c_eq/(9);
    size_t p1 = p0 - ny;
    size_t p2 = p0 - 1;
    rhs[row    ] = utilde[p1] - 2.0 * utilde[p0] + utilde[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    values[c_eq     ] = 0.0;
    values[c_eq +  1] = 0.0;
    values[c_eq +  2] = 0.0;

    values[c_eq +  3] = 0.0;
    values[c_eq +  4] = 0.0;
    values[c_eq +  5] = -theta;

    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = -theta;
    values[c_eq +  8] = 2.0 * theta;

    return 0;
}
int reg_corner_south_east_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, struct _grid_metric & metric)
{
    size_t ny = metric.ny;

    size_t p0 = c_eq/(9);
    size_t p1 = p0 - ny;
    size_t p2 = p0 + 1;
    rhs[row    ] = utilde[p1] - 2.0 * utilde[p0] + utilde[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    values[c_eq     ] = 0.0;
    values[c_eq +  1] = 0.0;
    values[c_eq +  2] = 0.0;

    values[c_eq +  3] = -theta;
    values[c_eq +  4] = 0.0;
    values[c_eq +  5] = 0.0;

    values[c_eq +  6] = 2.0 * theta;
    values[c_eq +  7] = -theta;
    values[c_eq +  8] = 0.0;

    return 0;
}
int reg_corner_south_west_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, struct _grid_metric & metric)
{
    size_t ny = metric.ny;

    size_t p0 = c_eq/(9);
    size_t p1 = p0 + 1;
    size_t p2 = p0 + ny;
    rhs[row    ] = utilde[p1] - 2.0 * utilde[p0] + utilde[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    values[c_eq     ] = 2.0 * theta;
    values[c_eq +  1] = -theta;
    values[c_eq +  2] = 0.0;

    values[c_eq +  3] = -theta;
    values[c_eq +  4] = 0.0;
    values[c_eq +  5] = 0.0;

    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;

    return 0;
}
int reg_corner_north_west_utilde(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, struct _grid_metric & metric)
{
    size_t ny = metric.ny;

    size_t p0 = c_eq/(9);
    size_t p1 = p0 - 1;
    size_t p2 = p0 + ny;
    rhs[row    ] = utilde[p1] - 2.0 * utilde[p0] + utilde[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    values[c_eq     ] = 0.0;
    values[c_eq +  1] = -theta;
    values[c_eq +  2] = 2.0 * theta;

    // w
    values[c_eq +  3] = 0.0;   
    values[c_eq +  4] = 0.0;   
    values[c_eq +  5] = -theta;

    // nw
    values[c_eq +  6] = 0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;

    return 0;
}
