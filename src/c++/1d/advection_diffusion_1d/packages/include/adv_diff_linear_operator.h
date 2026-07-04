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
#ifndef __ADV_DIFF_LINEAR_OPERATOR_H__
#define __ADV_DIFF_LINEAR_OPERATOR_H__

#include <cstdlib>
#include <vector>

void adv_diff_linear_operator(std::vector<double>& la, std::vector<double>& lb, std::vector<double>& lc, const std::vector<double> u, const std::vector<double> eps, const double dx, const int);

#endif __ADV_DIFF_LINEAR_OPERATOR_H__
