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

#ifndef __VISCOSITY_H__
#define __VISCOSITY_H__

#include <vector>
#include <Eigen/Sparse>

int viscosity_matrix_and_rhs(double* values, int row, int c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& Ttheta, 
    std::vector<double>& visc, double theta, double dx, double dy, int nx, int ny);
int viscosity_post_rhs(std::vector<double>& rhs_q, std::vector<double>& Tn,
    std::vector<double>& visc, double dx, double dy, int nx, int ny);

inline int viscosity_idx(int i, int j, int nx);
inline void add_value(double * values, int col, double data);

#endif  // __VISCOSITY_H__
