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
#ifndef __ADV_DIFF_INIT_CONCENTRATION_H__
#define __ADV_DIFF_INIT_CONCENTRATION_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>

enum class SHAPE_CONC
{
    NONE = 0,
    Constant,
    Envelope,
    NR_SHAPES
};

void adv_diff_init_concentration(std::vector<double>& mass, std::vector<double>& x, double Lx, SHAPE_CONC shape, std::vector<double>& d);
void control_volumes(std::vector<double>& u_ana, std::vector<double>& cv, double dx, size_t refine);
void compatible_function(std::vector<double>& mass, std::vector<double>& cv, std::vector<double>& u_ana, std::vector<double>& u_out, 
    double dx, size_t refine);


#endif __ADV_DIFF_INIT_CONCENTRATION_H__
