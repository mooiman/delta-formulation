//
//    Solving the 1D advection/diffusion equation, fully implicit with delta-formuation and Modified Newton iteration 
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
#ifndef __ADV_DIFF_REGULARIZATION__
#define __ADV_DIFF_REGULARIZATION__
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

void adv_diff_regularization(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& u_giv, double dx, double c_psi);

#endif __ADV_DIFF_REGULARIZATION__