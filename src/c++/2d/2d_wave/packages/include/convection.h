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

#ifndef __CONVECTION_H__
#define __CONVECTION_H__

#include <vector>
#include <Eigen/Sparse>

int convection_matrix_and_rhs(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double theta, size_t nx, size_t ny);
int convection_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    size_t nx, size_t ny );                // RHS vector [h, q, r]^{n}

inline size_t convection_idx(size_t i, size_t j, size_t nx);
inline void add_value(double * values, size_t col, double data);

#endif  // __CONVECTION_H__
