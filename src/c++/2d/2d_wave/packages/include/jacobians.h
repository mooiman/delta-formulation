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

#ifndef __JACOBIANS_H__
#define __JACOBIANS_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define UNUSED(x) (void)x;

double bed_shear_stress_J_10(double& h, double& q, double& r, double& cf);  // right hand side q-equation
double bed_shear_stress_J_11(double& h, double& q, double& r, double& cf);
double bed_shear_stress_J_12(double& h, double& q, double& r, double& cf);
double bed_shear_stress_J_13(double& h, double& q, double& r, double& cf);
double bed_shear_stress_J_20(double& h, double& q, double& r, double& cf);  // right hand side r-equation
double bed_shear_stress_J_21(double& h, double& q, double& r, double& cf);
double bed_shear_stress_J_22(double& h, double& q, double& r, double& cf);
double bed_shear_stress_J_23(double& h, double& q, double& r, double& cf);
double abs_vecq(double& qp, double& rp, double a);

double convection_J_10(double& h, double& q, double& r, double nxi, double neta);
double convection_J_20(double& h, double& q, double& r, double nxi, double neta);
double convection_J_11(double& h, double& q, double& r, double nxi, double neta);
double convection_J_21(double& h, double& q, double& r, double nxi, double neta);
double convection_J_12(double& h, double& q, double& r, double nxi, double neta);
double convection_J_22(double& h, double& q, double& r, double nxi, double neta);
double convection_J_13(double& h, double& q, double& r, double nxi, double neta);
double convection_J_23(double& h, double& q, double& r, double nxi, double neta);

#endif  // __JACOBIANS_H__
