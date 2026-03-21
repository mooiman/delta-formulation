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

#ifndef __MATRIX_ASSEMBLY_PSI_INTERIOR_H__
#define __MATRIX_ASSEMBLY_PSI_INTERIOR_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "grid_metric.h"

int reg_interior_matrix_psi(double* values, size_t row, size_t c_eq,
    double c_psi, struct _grid_metric metric);
int reg_interior_rhs_psi( size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& h, std::vector<double>& q, std::vector<double>& r,
    double c_psi, double g, struct _grid_metric & metric);

inline double F1(std::vector<double> & u, std::vector<size_t>& p, struct _grid_metric & metric );
inline double F2(std::vector<double> & u, std::vector<size_t>& p, struct _grid_metric & metric );
inline double F3(std::vector<double> & u, std::vector<size_t>& p, struct _grid_metric & metric );

inline double d2udxi2(std::vector<double> & u, std::vector<size_t>& p);
inline double d2udxideta(std::vector<double> & u, std::vector<size_t>& p);
inline double d2udeta2(std::vector<double> & u, std::vector<size_t>& p);

inline void add_value(double * values, size_t col, double data);

#endif
