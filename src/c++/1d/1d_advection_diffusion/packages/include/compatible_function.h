//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D advection/diffusion equation, fully implicit with delta-formulation and Modified Newton iteration 
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
#ifndef __COMPATIBLE_FUNCTION_H__
#define __COMPATIBLE_FUNCTION_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>

// for bicgstab  solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

void compatible_function(std::vector<double>& x, std::vector<double>& u_out, std::vector<double>& u_in, std::vector<double>& mass);

#endif __COMPATIBLE_FUNCTION_H__
