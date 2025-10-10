//
// Programmer: Jan Mooiman
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

#ifndef __MATRIX_ASSEMBLY_CORNERS_H__
#define __MATRIX_ASSEMBLY_CORNERS_H__

#include <vector>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

int corner_north_east(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta
    );
int corner_south_east(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta
    );
int corner_south_west(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta
    );
int corner_north_west(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double theta, size_t nx, size_t ny,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta
    );
inline void set_value(double * values, size_t col, double data)
{ 
    values[col] += data; 
}

#endif  // __MATRIX_ASSEMBLY_CORNERS_H__
