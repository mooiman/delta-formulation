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
int corner_north_east(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, int nx, int ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta)
{
    int p0 = c_eq/(3*27);
    int p1 = p0 - ny;
    int p2 = p0 - 1;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];
    rhs[row + 1] = qtheta[p1] - 2.0 * qtheta[p0] + qtheta[p2];
    rhs[row + 2] = rtheta[p1] - 2.0 * rtheta[p0] + rtheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // sw
    values[c_eq     ] = 0.0;   // Delta h
    values[c_eq +  1] = 0.0;   // Delta q
    values[c_eq +  2] = 0.0;   // Delta r

    // w
    values[c_eq +  3] = 0.0;   // Delta h
    values[c_eq +  4] = 0.0;   // Delta q
    values[c_eq +  5] = 0.0;   // Delta r

    // nw
    values[c_eq +  6] = 0.0;   // Delta h
    values[c_eq +  7] = 0.0;   // Delta q
    values[c_eq +  8] = 0.0;   // Delta r

    // s
    values[c_eq +  9] = 0.0;   // Delta h
    values[c_eq + 10] = 0.0;   // Delta q
    values[c_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[c_eq + 12] = 0.0;   // Delta h
    values[c_eq + 13] = 0.0;   // Delta q
    values[c_eq + 14] = 0.0;   // Delta r

    // n
    values[c_eq + 15] = -theta;   // Delta h
    values[c_eq + 16] = 0.0;   // Delta q
    values[c_eq + 17] = 0.0;   // Delta r

    // se
    values[c_eq + 18] = 0.0;   // Delta h
    values[c_eq + 19] = 0.0;   // Delta q
    values[c_eq + 20] = 0.0;   // Delta r

    // e
    values[c_eq + 21] = -theta;   // Delta h
    values[c_eq + 22] = 0.0;   // Delta q
    values[c_eq + 23] = 0.0;   // Delta r

    // ne
    values[c_eq + 24] = 2.0 * theta;   // Delta h
    values[c_eq + 25] = 0.0;   // Delta q
    values[c_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // q-equation
    // 
    // sw
    values[q_eq     ] = 0.0;   // Delta h
    values[q_eq +  1] = 0.0;   // Delta q
    values[q_eq +  2] = 0.0;   // Delta r

    // w   
    values[q_eq +  3] = 0.0;   // Delta h
    values[q_eq +  4] = 0.0;   // Delta q
    values[q_eq +  5] = 0.0;   // Delta r

    // nw
    values[q_eq +  6] = 0.0;   // Delta h
    values[q_eq +  7] = 0.0;   // Delta q
    values[q_eq +  8] = 0.0;   // Delta r

    // s
    values[q_eq +  9] = 0.0;   // Delta h
    values[q_eq + 10] = 0.0;   // Delta q
    values[q_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[q_eq + 12] = 0.0;   // Delta h
    values[q_eq + 13] = 0.0;   // Delta q
    values[q_eq + 14] = 0.0;   // Delta r

    // n
    values[q_eq + 15] = 0.0;   // Delta h
    values[q_eq + 16] = -theta;   // Delta q
    values[q_eq + 17] = 0.0;   // Delta r

    // se
    values[q_eq + 18] = 0.0;   // Delta h
    values[q_eq + 19] = 0.0;   // Delta q
    values[q_eq + 20] = 0.0;   // Delta r

    // e
    values[q_eq + 21] = 0.0;   // Delta h
    values[q_eq + 22] = -theta;   // Delta q
    values[q_eq + 23] = 0.0;   // Delta r

    // ne
    values[q_eq + 24] = 0.0;   // Delta h
    values[q_eq + 25] = 2.0 * theta;   // Delta q
    values[q_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // r-equation
    // 
    // sw
    values[r_eq     ] = 0.0;   // Delta h
    values[r_eq +  1] = 0.0;   // Delta q
    values[r_eq +  2] = 0.0;   // Delta r

    // w   
    values[r_eq +  3] = 0.0;   // Delta h
    values[r_eq +  4] = 0.0;   // Delta q
    values[r_eq +  5] = 0.0;   // Delta r

    // nw
    values[r_eq +  6] = 0.0;   // Delta h
    values[r_eq +  7] = 0.0;   // Delta q
    values[r_eq +  8] = 0.0;   // Delta r

    // s
    values[r_eq +  9] = 0.0;   // Delta h
    values[r_eq + 10] = 0.0;   // Delta q
    values[r_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[r_eq + 12] = 0.0;   // Delta h
    values[r_eq + 13] = 0.0;   // Delta q
    values[r_eq + 14] = 0.0;   // Delta r

    // n
    values[r_eq + 15] = 0.0;   // Delta h
    values[r_eq + 16] = 0.0;   // Delta q
    values[r_eq + 17] = -theta;   // Delta r

    // se
    values[r_eq + 18] = 0.0;   // Delta h
    values[r_eq + 19] = 0.0;   // Delta q
    values[r_eq + 20] = 0.0;   // Delta r

    // e
    values[r_eq + 21] = 0.0;   // Delta h
    values[r_eq + 22] = 0.0;   // Delta q
    values[r_eq + 23] = -theta;   // Delta r

    // ne
    values[r_eq + 24] = 0.0;   // Delta h
    values[r_eq + 25] = 0.0;   // Delta q
    values[r_eq + 26] = 2.0 * theta;   // Delta r

    return 0;
}
int corner_south_east(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, int nx, int ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta)
{
    int p0 = c_eq/(3*27);
    int p1 = p0 - ny;
    int p2 = p0 + 1;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];
    rhs[row + 1] = qtheta[p1] - 2.0 * qtheta[p0] + qtheta[p2];
    rhs[row + 2] = rtheta[p1] - 2.0 * rtheta[p0] + rtheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // sw
    values[c_eq     ] = 0.0;   // Delta h
    values[c_eq +  1] = 0.0;   // Delta q
    values[c_eq +  2] = 0.0;   // Delta r

    // w
    values[c_eq +  3] = 0.0;   // Delta h
    values[c_eq +  4] = 0.0;   // Delta q
    values[c_eq +  5] = 0.0;   // Delta r

    // nw
    values[c_eq +  6] = 0.0;   // Delta h
    values[c_eq +  7] = 0.0;   // Delta q
    values[c_eq +  8] = 0.0;   // Delta r

    // s
    values[c_eq +  9] = -theta;   // Delta h
    values[c_eq + 10] = 0.0;   // Delta q
    values[c_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[c_eq + 12] = 0.0;   // Delta h
    values[c_eq + 13] = 0.0;   // Delta q
    values[c_eq + 14] = 0.0;   // Delta r

    // n
    values[c_eq + 15] = 0.0;   // Delta h
    values[c_eq + 16] = 0.0;   // Delta q
    values[c_eq + 17] = 0.0;   // Delta r

    // se
    values[c_eq + 18] = 2.0 * theta;   // Delta h
    values[c_eq + 19] = 0.0;   // Delta q
    values[c_eq + 20] = 0.0;   // Delta r

    // e
    values[c_eq + 21] = -theta;   // Delta h
    values[c_eq + 22] = 0.0;   // Delta q
    values[c_eq + 23] = 0.0;   // Delta r

    // ne
    values[c_eq + 24] = 0.0;   // Delta h
    values[c_eq + 25] = 0.0;   // Delta q
    values[c_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // q-equation
    // 
    // sw
    values[q_eq     ] = 0.0;   // Delta h
    values[q_eq +  1] = 0.0;   // Delta q
    values[q_eq +  2] = 0.0;   // Delta r

    // w   
    values[q_eq +  3] = 0.0;   // Delta h
    values[q_eq +  4] = 0.0;   // Delta q
    values[q_eq +  5] = 0.0;   // Delta r

    // nw
    values[q_eq +  6] = 0.0;   // Delta h
    values[q_eq +  7] = 0.0;   // Delta q
    values[q_eq +  8] = 0.0;   // Delta r

    // s
    values[q_eq +  9] = 0.0;   // Delta h
    values[q_eq + 10] = -theta;   // Delta q
    values[q_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[q_eq + 12] = 0.0;   // Delta h
    values[q_eq + 13] = 0.0;   // Delta q
    values[q_eq + 14] = 0.0;   // Delta r

    // n
    values[q_eq + 15] = 0.0;   // Delta h
    values[q_eq + 16] = 0.0;   // Delta q
    values[q_eq + 17] = 0.0;   // Delta r

    // se
    values[q_eq + 18] = 0.0;   // Delta h
    values[q_eq + 19] = 2.0 * theta;   // Delta q
    values[q_eq + 20] = 0.0;   // Delta r

    // e
    values[q_eq + 21] = 0.0;   // Delta h
    values[q_eq + 22] = -theta;   // Delta q
    values[q_eq + 23] = 0.0;   // Delta r

    // ne
    values[q_eq + 24] = 0.0;   // Delta h
    values[q_eq + 25] = 0.0;   // Delta q
    values[q_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // r-equation
    // 
    // sw
    values[r_eq     ] = 0.0;   // Delta h
    values[r_eq +  1] = 0.0;   // Delta q
    values[r_eq +  2] = 0.0;   // Delta r

    // w   
    values[r_eq +  3] = 0.0;   // Delta h
    values[r_eq +  4] = 0.0;   // Delta q
    values[r_eq +  5] = 0.0;   // Delta r

    // nw
    values[r_eq +  6] = 0.0;   // Delta h
    values[r_eq +  7] = 0.0;   // Delta q
    values[r_eq +  8] = 0.0;   // Delta r

    // s
    values[r_eq +  9] = 0.0;   // Delta h
    values[r_eq + 10] = 0.0;   // Delta q
    values[r_eq + 11] = -theta;   // Delta r

    // c centre point of discretisation molecule
    values[r_eq + 12] = 0.0;   // Delta h
    values[r_eq + 13] = 0.0;   // Delta q
    values[r_eq + 14] = 0.0;   // Delta r

    // n
    values[r_eq + 15] = 0.0;   // Delta h
    values[r_eq + 16] = 0.0;   // Delta q
    values[r_eq + 17] = 0.0;   // Delta r

    // se
    values[r_eq + 18] = 0.0;   // Delta h
    values[r_eq + 19] = 0.0;   // Delta q
    values[r_eq + 20] = 2.0 * theta;   // Delta r

    // e
    values[r_eq + 21] = 0.0;   // Delta h
    values[r_eq + 22] = 0.0;   // Delta q
    values[r_eq + 23] = -theta;   // Delta r

    // ne
    values[r_eq + 24] = 0.0;   // Delta h
    values[r_eq + 25] = 0.0;   // Delta q
    values[r_eq + 26] = 0.0;   // Delta r

    return 0;
}
int corner_south_west(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, int nx, int ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta)
{
    int p0 = c_eq/(3*27);
    int p1 = p0 + 1;
    int p2 = p0 + ny;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];
    rhs[row + 1] = qtheta[p1] - 2.0 * qtheta[p0] + qtheta[p2];
    rhs[row + 2] = rtheta[p1] - 2.0 * rtheta[p0] + rtheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // sw
    values[c_eq     ] = 2.0 * theta;   // Delta h
    values[c_eq +  1] = 0.0;   // Delta q
    values[c_eq +  2] = 0.0;   // Delta r

    // w
    values[c_eq +  3] = -theta;   // Delta h
    values[c_eq +  4] = 0.0;   // Delta q
    values[c_eq +  5] = 0.0;   // Delta r

    // nw
    values[c_eq +  6] = 0.0;   // Delta h
    values[c_eq +  7] = 0.0;   // Delta q
    values[c_eq +  8] = 0.0;   // Delta r

    // s
    values[c_eq +  9] = -theta;   // Delta h
    values[c_eq + 10] = 0.0;   // Delta q
    values[c_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[c_eq + 12] = 0.0;   // Delta h
    values[c_eq + 13] = 0.0;   // Delta q
    values[c_eq + 14] = 0.0;   // Delta r

    // n
    values[c_eq + 15] = 0.0;   // Delta h
    values[c_eq + 16] = 0.0;   // Delta q
    values[c_eq + 17] = 0.0;   // Delta r

    // se
    values[c_eq + 18] = 0.0;   // Delta h
    values[c_eq + 19] = 0.0;   // Delta q
    values[c_eq + 20] = 0.0;   // Delta r

    // s
    values[c_eq + 21] = 0.0;   // Delta h
    values[c_eq + 22] = 0.0;   // Delta q
    values[c_eq + 23] = 0.0;   // Delta r

    // ne
    values[c_eq + 24] = 0.0;   // Delta h
    values[c_eq + 25] = 0.0;   // Delta q
    values[c_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // q-equation
    // 
    // sw
    values[q_eq     ] = 0.0;   // Delta h
    values[q_eq +  1] = 2.0 * theta;   // Delta q
    values[q_eq +  2] = 0.0;   // Delta r

    // w   
    values[q_eq +  3] = 0.0;   // Delta h
    values[q_eq +  4] = -theta;   // Delta q
    values[q_eq +  5] = 0.0;   // Delta r

    // nw
    values[q_eq +  6] = 0.0;   // Delta h
    values[q_eq +  7] = 0.0;   // Delta q
    values[q_eq +  8] = 0.0;   // Delta r

    // s
    values[q_eq +  9] = 0.0;   // Delta h
    values[q_eq + 10] = -theta;   // Delta q
    values[q_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[q_eq + 12] = 0.0;   // Delta h
    values[q_eq + 13] = 0.0;   // Delta q
    values[q_eq + 14] = 0.0;   // Delta r

    // n
    values[q_eq + 15] = 0.0;   // Delta h
    values[q_eq + 16] = 0.0;   // Delta q
    values[q_eq + 17] = 0.0;   // Delta r

    // se
    values[q_eq + 18] = 0.0;   // Delta h
    values[q_eq + 19] = 0.0;   // Delta q
    values[q_eq + 20] = 0.0;   // Delta r

    // s
    values[q_eq + 21] = 0.0;   // Delta h
    values[q_eq + 22] = 0.0;   // Delta q
    values[q_eq + 23] = 0.0;   // Delta r

    // ne
    values[q_eq + 24] = 0.0;   // Delta h
    values[q_eq + 25] = 0.0;   // Delta q
    values[q_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // r-equation
    // 
    // sw
    values[r_eq     ] = 0.0;   // Delta h
    values[r_eq +  1] = 0.0;   // Delta q
    values[r_eq +  2] = 2.0 * theta;   // Delta r

    // w   
    values[r_eq +  3] = 0.0;   // Delta h
    values[r_eq +  4] = 0.0;   // Delta q
    values[r_eq +  5] = -theta;   // Delta r

    // nw
    values[r_eq +  6] = 0.0;   // Delta h
    values[r_eq +  7] = 0.0;   // Delta q
    values[r_eq +  8] = 0.0;   // Delta r

    // s
    values[r_eq +  9] = 0.0;   // Delta h
    values[r_eq + 10] = 0.0;   // Delta q
    values[r_eq + 11] = -theta;   // Delta r

    // c centre point of discretisation molecule
    values[r_eq + 12] = 0.0;   // Delta h
    values[r_eq + 13] = 0.0;   // Delta q
    values[r_eq + 14] = 0.0;   // Delta r

    // n
    values[r_eq + 15] = 0.0;   // Delta h
    values[r_eq + 16] = 0.0;   // Delta q
    values[r_eq + 17] = 0.0;   // Delta r

    // se
    values[r_eq + 18] = 0.0;   // Delta h
    values[r_eq + 19] = 0.0;   // Delta q
    values[r_eq + 20] = 0.0;   // Delta r

    // s
    values[r_eq + 21] = 0.0;   // Delta h
    values[r_eq + 22] = 0.0;   // Delta q
    values[r_eq + 23] = 0.0;   // Delta r

    // ne
    values[r_eq + 24] = 0.0;   // Delta h
    values[r_eq + 25] = 0.0;   // Delta q
    values[r_eq + 26] = 0.0;   // Delta r
    return 0;
}
int corner_north_west(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, int nx, int ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta)
{
    int p0 = c_eq/(3*27);
    int p1 = p0 - 1;
    int p2 = p0 + ny;
    rhs[row    ] = htheta[p1] - 2.0 * htheta[p0] + htheta[p2];
    rhs[row + 1] = qtheta[p1] - 2.0 * qtheta[p0] + qtheta[p2];
    rhs[row + 2] = rtheta[p1] - 2.0 * rtheta[p0] + rtheta[p2];

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // sw
    values[c_eq     ] = 0.0;   // Delta h
    values[c_eq +  1] = 0.0;   // Delta q
    values[c_eq +  2] = 0.0;   // Delta r

    // w
    values[c_eq +  3] = -theta;   // Delta h
    values[c_eq +  4] = 0.0;   // Delta q
    values[c_eq +  5] = 0.0;   // Delta r

    // nw
    values[c_eq +  6] = 2.0 * theta;   // Delta h
    values[c_eq +  7] = 0.0;   // Delta q
    values[c_eq +  8] = 0.0;   // Delta r

    // s
    values[c_eq +  9] = 0.0;   // Delta h
    values[c_eq + 10] = 0.0;   // Delta q
    values[c_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[c_eq + 12] = 0.0;   // Delta h
    values[c_eq + 13] = 0.0;   // Delta q
    values[c_eq + 14] = 0.0;   // Delta r

    // n
    values[c_eq + 15] = -theta;   // Delta h
    values[c_eq + 16] = 0.0;   // Delta q
    values[c_eq + 17] = 0.0;   // Delta r

    // se
    values[c_eq + 18] = 0.0;   // Delta h
    values[c_eq + 19] = 0.0;   // Delta q
    values[c_eq + 20] = 0.0;   // Delta r

    // e
    values[c_eq + 21] = 0.0;   // Delta h
    values[c_eq + 22] = 0.0;   // Delta q
    values[c_eq + 23] = 0.0;   // Delta r

    // ne
    values[c_eq + 24] = 0.0;   // Delta h
    values[c_eq + 25] = 0.0;   // Delta q
    values[c_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // q-equation
    // 
    // sw
    values[q_eq     ] = 0.0;   // Delta h
    values[q_eq +  1] = 0,0;   // Delta q
    values[q_eq +  2] = 0.0;   // Delta r

    // w   
    values[q_eq +  3] = 0.0;   // Delta h
    values[q_eq +  4] = -theta;   // Delta q
    values[q_eq +  5] = 0.0;   // Delta r

    // nw
    values[q_eq +  6] = 0.0;     // Delta h
    values[q_eq +  7] = 2.0 * theta;   // Delta q
    values[q_eq +  8] = 0.0;     // Delta r

    // s
    values[q_eq +  9] = 0.0;   // Delta h
    values[q_eq + 10] = 0.0;   // Delta q
    values[q_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[q_eq + 12] = 0.0;   // Delta h
    values[q_eq + 13] = 0.0;   // Delta q
    values[q_eq + 14] = 0.0;   // Delta r

    // n
    values[q_eq + 15] = 0.0;   // Delta h
    values[q_eq + 16] = -theta;   // Delta q
    values[q_eq + 17] = 0.0;   // Delta r

    // se
    values[q_eq + 18] = 0.0;   // Delta h
    values[q_eq + 19] = 0.0;   // Delta q
    values[q_eq + 20] = 0.0;   // Delta r

    // e
    values[q_eq + 21] = 0.0;   // Delta h
    values[q_eq + 22] = 0.0;   // Delta q
    values[q_eq + 23] = 0.0;   // Delta r

    // ne
    values[q_eq + 24] = 0.0;   // Delta h
    values[q_eq + 25] = 0.0;   // Delta q
    values[q_eq + 26] = 0.0;   // Delta r

    //--------------------------------------------------------------------------
    // r-equation
    // 
    // sw
    values[r_eq     ] = 0.0;   // Delta h
    values[r_eq +  1] = 0.0;   // Delta q
    values[r_eq +  2] = 0.0;   // Delta r

    // w   
    values[r_eq +  3] = 0.0;   // Delta h
    values[r_eq +  4] = 0.0;   // Delta q
    values[r_eq +  5] = -theta;   // Delta r

    // nw
    values[r_eq +  6] = 0.0;      // Delta h
    values[r_eq +  7] = 0.0;     // Delta q
    values[r_eq +  8] = 2.0 * theta;   // Delta r

    // s
    values[r_eq +  9] = 0.0;   // Delta h
    values[r_eq + 10] = 0.0;   // Delta q
    values[r_eq + 11] = 0.0;   // Delta r

    // c centre point of discretisation molecule
    values[r_eq + 12] = 0.0;   // Delta h
    values[r_eq + 13] = 0.0;   // Delta q
    values[r_eq + 14] = 0.0;   // Delta r

    // n
    values[r_eq + 15] = 0.0;   // Delta h
    values[r_eq + 16] = 0.0;   // Delta q
    values[r_eq + 17] = -theta;   // Delta r

    // se
    values[r_eq + 18] = 0.0;   // Delta h
    values[r_eq + 19] = 0.0;   // Delta q
    values[r_eq + 20] = 0.0;   // Delta r
                                
    // e                        
    values[r_eq + 21] = 0.0;   // Delta h
    values[r_eq + 22] = 0.0;   // Delta q
    values[r_eq + 23] = 0.0;   // Delta r
                                
    // ne                       
    values[r_eq + 24] = 0.0;   // Delta h
    values[r_eq + 25] = 0.0;   // Delta q
    values[r_eq + 26] = 0.0;   // Delta r
    return 0;
}
        
inline int ma_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}

