//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D wave equation, fully implicit with delta-formulation and Modified Newton iteration 
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

#pragma once

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>

struct _grid_metric {
    size_t nx;
    size_t ny;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> dx_dxi ; 
    std::vector<double> dy_dxi ; 
    std::vector<double> dx_deta; 
    std::vector<double> dy_deta;
    std::vector<double> ddx_dxi2 ; 
    std::vector<double> ddy_dxi2;
    std::vector<double> ddx_deta2; 
    std::vector<double> ddy_deta2; 
};
int grid_metric(struct _grid_metric & metric); 
inline size_t idx(size_t i, size_t j, size_t ny_in);
