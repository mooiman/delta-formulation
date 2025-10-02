//
// programmer: Jan Mooiman
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

#include "interpolations.h"

double bed_shear_stress_scv(double c0, double c1, double c2, double c3)
{
    // value at the quadrature point of a sub control volume
// scv 0 c_{i-1/4, j-1/4}
// scv 1 c_{i+1/4, j-1/4}
// scv 2 c_{i+1/4, j+1/4}
// scv 3 c_{i-1/4, j+1/4}

//  2 - - - - - - - 3 3 - - - - - - - 2
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  | - - - - - - - | | - - - - - - - | 
//  |       |       | |       |       | 
//  |       |   x   | |   x   |       | 
//  |       |       | |       |       | 
//  1 - - - - - - - 0 0 - - - - - - - 1 
//  1 - - - - - - - 0 0 - - - - - - - 1
//  |       |       | |       |       | 
//  |       |   x   | |   x   |       | 
//  |       |       | |       |       | 
//  | - - - - - - - | | - - - - - - - | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  2 - - - - - - - 3 3 - - - - - - - 2 

    return 1./16. * (9.0 * c0 + 3.0 * c1 + 1.0 * c2 + 3.0 * c3);
}

double scvf_xi(double c0, double c1, double c2, double c3)
{
// face 0 c_{i-1/2, j-1/4}
// face 3 c_{i+1/2, j-1/4}
// face 4 c_{i+1/2, j+1/4}
// face 7 c_{i-1/2, j+1/4}

//  2 - - - - - - - 3 3 - - - - - - - 2
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  | - - - - - - - | | - - - - - - - | 
//  |       |       | |       |       | 
//  |       x       | |       x       | 
//  |       |       | |       |       | 
//  1 - - - - - - - 0 0 - - - - - - - 1 
//  1 - - - - - - - 0 0 - - - - - - - 1
//  |       |       | |       |       | 
//  |       x       | |       x       | 
//  |       |       | |       |       | 
//  | - - - - - - - | | - - - - - - - | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  |       |       | |       |       | 
//  2 - - - - - - - 3 3 - - - - - - - 2 

    return 0.125 * (3. * c0 + 3. * c1 + 1. * c2 + 1. * c3);
}
double scvf_eta(double c0, double c1, double c2, double c3)
{
// face 1 c_{i-1/4, j-1/2}
// face 2 c_{i+1/4, j-1/2}
// face 5 c_{i+1/4, j+1/2}
// face 6 c_{i-1/4, j+1/2}

//   2 - - - - - - - 1 1 - - - - - - - 2
//   |       |       | |       |       |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   | - - - - - x - | | - x - - - - - |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   3 - - - - - - - 0 0 - - - - - - - 3
//   3 - - - - - - - 0 0 - - - - - - - 3
//   |       |       | |       |       |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   | - - - - - x - | | - x - - - - - |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   |       |       | |       |       |
//   2 - - - - - - - 1 1 - - - - - - - 2

    return 0.125 * (3. * c0 + 3. * c1 + 1. * c2 + 1. * c3);
}
double c_scv(double c0, double c1, double c2, double c3)
{
    // value at subcontrol volume
    return 0.0625 * (9. * c0 + 3. * c1 +  3. * c2 + c3);
}
double dcdx_scv(double c0, double c1, double c2, double c3)
{
    // value quadrature point (i+1/4, j+1/4) at subcontrol volume
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
double dcdy_scv(double c0, double c1, double c2, double c3)
{
    // value quadrature point (i+1/4, j+1/4) at subcontrol volume
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
double dcdx_scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdx normal at subcontrol volume edge
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
double dcdx_scvf_t(double c0, double c1, double c2, double c3)
{
    // dcdx tangential at subcontrol volume edge
    return 0.5 * (c0 - c1 + c2 - c3);
}
double dcdy_scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdy normal at subcontrol volume edge
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
double dcdy_scvf_t(double c0, double c1, double c2, double c3)
{
    // dcdy tangential at subcontrol volume edge
    return 0.5 * (c0 - c1 + c2 - c3);
}

