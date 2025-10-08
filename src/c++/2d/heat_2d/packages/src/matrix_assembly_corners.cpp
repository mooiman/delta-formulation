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

#include "matrix_assembly_corners.h"

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
int corner_north_east(double* values, size_t row, int c_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta)
{
    size_t p0 = c_eq/(9);
    size_t p1 = p0 - ny;
    int p2 = p0 - 1;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];

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
int corner_south_east(double* values, size_t row, int c_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta)
{
    size_t p0 = c_eq/(9);
    size_t p1 = p0 - ny;
    size_t p2 = p0 + 1;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];

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
int corner_south_west(double* values, size_t row, int c_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta)
{
    size_t p0 = c_eq/(9);
    size_t p1 = p0 + 1;
    size_t p2 = p0 + ny;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // sw
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
int corner_north_west(double* values, size_t row, int c_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta)
{
    size_t p0 = c_eq/(9);
    size_t p1 = p0 - 1;
    size_t p2 = p0 + ny;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    values[c_eq     ] = 0.0;
    values[c_eq +  1] = -theta;
    values[c_eq +  2] = 2.0 * theta;

    values[c_eq +  3] = 0.0;
    values[c_eq +  4] = 0.0;   
    values[c_eq +  5] = -theta;   

    values[c_eq +  6] = 0.0;
    values[c_eq +  7] = 0.0;
    values[c_eq +  8] = 0.0;

    return 0;
}
        
inline size_t ma_index(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}

