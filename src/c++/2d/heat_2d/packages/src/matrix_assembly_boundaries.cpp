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

#include "interpolations.h"
#include "matrix_assembly_boundaries.h"

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
// boundary nodes
//
//==============================================================================
int boundary_north(double* values, int row, int c_eq, Eigen::VectorXd& rhs, 
    double & dtinv,double & theta, double eps_bc_corr, 
    bool stationary, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& Tn, 
    std::vector<double>& Tp, 
    std::vector<double>& Ttheta,
    std::vector<std::string> bc_type,  std::vector<std::string> bc_vars, int BC_NORTH, std::vector<double> bc)
{
    memset(&values[c_eq], 0, 9 * sizeof(double));  // set all coefficients for one row of r-equation to zero

    int p_5 = c_eq/(9);  // node number of boundary point, ie north point of molecule
    int p_4 = p_5 - 1;
    int p_3 = p_5 - 2;
    int p_0 = p_5 - ny - 2;
    int p_1 = p_5 - ny - 1;
    int p_2 = p_5 - ny;
    int p_6 = p_5 + ny - 2;
    int p_7 = p_5 + ny - 1;
    int p_8 = p_5 + ny;

    double T_given = bc[BC_NORTH] ;

    // nnorth
    int col_b  = c_eq + 5 * 1;
    int col_s  = c_eq + 4 * 1;
    int col_ss = c_eq + 3 * 1;
    if (bc_type[BC_NORTH] == "neumann")
    {
        // Contribution Delta h
        values[col_b ] =  1.0 * theta;
        values[col_s ] = -1.0 * theta;
        values[col_ss] = 0.0;

        rhs[row] = -(Ttheta[p_5] - Ttheta[p_4]);
    }
    if (bc_type[BC_NORTH] == "dirichlet")
    {
        // Contribution Delta h
        values[col_b ] = 1.0 * theta;
        values[col_s ] = 0.0;
        values[col_ss] = 0.0;

        rhs[row] = -(Ttheta[p_5] - bc[BC_NORTH]);
    }
    return 0;
}
//==============================================================================
int boundary_east(double* values, int row, int c_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double eps_bc_corr, 
    bool stationary, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& Tn, 
    std::vector<double>& Tp, 
    std::vector<double>& Ttheta,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_EAST, std::vector<double> bc)
{
    memset(&values[c_eq], 0, 9 * sizeof(double));  // set all coefficients for one row r-equation to zero

    int p_7 = c_eq/(9);  // node number of boundary point, ie east point of molecule
    int p_8 = p_7 + 1;
    int p_6 = p_7 - 1;
    int p_5 = p_7 - ny + 1;
    int p_4 = p_7 - ny;
    int p_3 = p_7 - ny - 1;
    int p_2 = p_7 - 2 * ny + 1;
    int p_1 = p_7 - 2 * ny;
    int p_0 = p_7 - 2 * ny - 1;

    double T_given = bc[BC_EAST] ;

    // eeast
    int col_b  = c_eq + 7 * 1;
    int col_w  = c_eq + 4 * 1;
    int col_ww = c_eq + 1 * 1;
    if (bc_type[BC_EAST] == "neumann")
    {
        // Contribution Delta h
        values[col_b ] =  1.0 * theta;
        values[col_w ] = -1.0 * theta;
        values[col_ww] = 0.0;

        rhs[row] = -(Ttheta[p_7] - Ttheta[p_4]);
    }
    if (bc_type[BC_EAST] == "dirichlet")
    {
        // Contribution Delta h
        values[col_b ] = 1.0 * theta;
        values[col_w ] = 0.0;
        values[col_ww] = 0.0;

        rhs[row] = -(Ttheta[p_7] - bc[BC_EAST]);
    }
    return 0;
}
//==============================================================================
int boundary_south(double* values, int row, int c_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double eps_bc_corr, 
    bool stationary, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& Tn, 
    std::vector<double>& Tp, 
    std::vector<double>& Ttheta,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_SOUTH, std::vector<double> bc)
{
    memset(&values[c_eq], 0, 9 * sizeof(double));  // set all coefficients for one row of r-equation to zero

    int p_3 = c_eq/(9);  // node number of boundary point, ie south point of molecule
    int p_4 = p_3 + 1;
    int p_5 = p_3 + 2;
    int p_0 = p_3 - ny;
    int p_1 = p_3 - ny + 1;
    int p_2 = p_3 - ny + 2;
    int p_6 = p_3 + ny;
    int p_7 = p_3 + ny + 1;
    int p_8 = p_3 + ny + 2;

    double T_given = bc[BC_SOUTH];
    // ssouth
    int col_b  = c_eq + 3 * 1;
    int col_n  = c_eq + 4 * 1;
    int col_nn = c_eq + 5 * 1;
    if (bc_type[BC_SOUTH] == "neumann")
    {
        // Contribution Delta h
        values[col_b ] = -1.0 * theta;
        values[col_n ] =  1.0 * theta;
        values[col_nn] = 0.0;

        rhs[row] = -(Ttheta[p_4]- Ttheta[p_3]);
    }
    if (bc_type[BC_SOUTH] == "dirichlet")
    {
        // Contribution Delta h
        values[col_b ] = 1.0 * theta;
        values[col_n ] = 0.0;
        values[col_nn] = 0.0;

        rhs[row] = -(Ttheta[p_3] - bc[BC_SOUTH]);
    }
    return 0;
}
//==============================================================================
int boundary_west(double* values, int row, int c_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double eps_bc_corr, 
    bool stationary, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& Tn, 
    std::vector<double>& Tp, 
    std::vector<double>& Ttheta,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_WEST, std::vector<double> bc)
{
    memset(&values[c_eq], 0, 9 * sizeof(double));  // set all coefficients for one row of r-equation to zero

    int p_1 = c_eq/(9);  // node number of boundary point, ie west point of molucule
    int p_0 = p_1 - 1;
    int p_2 = p_1 + 1;
    int p_3 = p_1 + ny - 1;
    int p_4 = p_1 + ny;
    int p_5 = p_1 + ny + 1;
    int p_6 = p_1 + 2 * ny - 1;
    int p_7 = p_1 + 2 * ny;
    int p_8 = p_1 + 2 * ny + 1;

    double T_given = bc[BC_WEST];
    // wwest
    int col_b  = c_eq + 1 * 1; // point of boundary, ie west point of molecule
    int col_e  = c_eq + 4 * 1;
    int col_ee = c_eq + 7 * 1;
    if (bc_type[BC_WEST] == "neumann")
    {
        // Contribution Delta h
        values[col_b ] = -1.0 * theta;
        values[col_e ] =  1.0 * theta;
        values[col_ee] = 0.0;

        rhs[row] = -(Ttheta[p_4] - Ttheta[p_1]);
    }
    if (bc_type[BC_WEST] == "dirichlet")
    {
        values[col_b ] = 1.0 * theta;
        values[col_e ] = 0.0;
        values[col_ee] = 0.0;

        rhs[row] = -(Ttheta[p_1] - bc[BC_WEST]);
    }
    return 0;
}
//==============================================================================
void  molecule(std::vector<int>& p, int p_sw, int ny)
{
        p[0] = p_sw;
        p[1] = p_sw + 1;
        p[2] = p_sw + 2;
        p[3] = p_sw + ny;
        p[4] = p_sw + ny + 1;
        p[5] = p_sw + ny + 2;
        p[6] = p_sw + 2 * ny;
        p[7] = p_sw + 2 * ny + 1;
        p[8] = p_sw + 2 * ny + 2;
}
inline void set_value(double * values, int col, double data)
{ 
    values[col] += data; 
}
inline int ma_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}

