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

#ifndef __DIFFUSION2D_H__
#define __DIFFUSION2D_H__

#include <vector>
#include <Eigen/Sparse>

#include "perf_timer.h"

int diffusion_matrix_and_rhs(double* values, size_t row, size_t c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& u_giv, double psi_11, double psi_22, 
    double theta, size_t nx, size_t ny);
int diffusion_post_rhs(std::vector<double>& rhs_q, 
    std::vector<double>& Tn, std::vector<double>& psi_11, std::vector<double>& psi_22, 
    size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y);

inline size_t diffusion_idx(size_t i, size_t j, size_t nx);
inline void add_value(double * values, size_t col, double data);

#endif