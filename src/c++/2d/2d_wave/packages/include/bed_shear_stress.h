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

#ifndef __BEDSHEARSTRESS_H__
#define __BEDSHEARSTRESS_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

int bed_shear_stress_matrix_and_rhs(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double cf, double theta, double dx, double dy, int nx, int ny);
void bed_shear_stress_post_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double cf, int nx, int ny);

    inline int bed_shear_stress_idx(int i, int j, int ny);
    inline void add_value(double * values, int col, double data);

#endif  // __BEDSHEARSTRESS_H__
