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
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Righthand side 0f  the Advection-Diffusion equation
//
//   L_h C = d/dx (u.C - eps.d(C)/dx) =q
//

#include "adv_diff_source.h"

void adv_diff_source(std::vector<double> & q, int select)
{
    switch (select)
    {
    case 1:
        //
        // Homogeneous righthand side
        //
        for (size_t i = 0; i < q.size(); ++i)
        {
            q[i] = 0.0;
        }
        break;
    }
}
