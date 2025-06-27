//
// Programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
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

#ifndef __INITIAL_CONDITIONS_H__
#define __INITIAL_CONDITIONS_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>

void initial_conditions(std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny, 
    std::vector<double>& s, std::vector<double>& u, std::vector<double>& v,
    std::vector<std::string>& ini_vars, double gauss_amp, double gauss_mu_x, double gauss_mu_y, double gauss_sigma_x, double gauss_sigma_y);

inline int p_index(int i, int j, int ny);

#endif __INITIAL_CONDITIONS_H__
