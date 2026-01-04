//
// Programmer: Jan Mooiman
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

#ifndef __MATRIX_ASSEMBLY_PSI_CORNERS_H__
#define __MATRIX_ASSEMBLY_PSI_CORNERS_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

int reg_corner_north_east_psi(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, size_t nx, size_t ny);

int reg_corner_south_east_psi(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, size_t nx, size_t ny);

int reg_corner_south_west_psi(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, size_t nx, size_t ny);

int reg_corner_north_west_psi(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& utilde, double theta, size_t nx, size_t ny);

#endif
