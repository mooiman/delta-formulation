//---------------------------------------------------------------
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
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
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Initial concentration for the Advection-Diffusion equation
//
//   L_h C = d/dx (u.C - eps.d(C)/dx)
//

#include "adv_diff_init_concentration.h"
#include "adv_diff_init_velocity.h"

void adv_diff_init_velocity(std::vector<double>& u, const double g,  const std::vector<double>& zb, const std::vector<double>& x, SHAPE_CONC shape_conc)
{
    switch (shape_conc)
    {
    case SHAPE_CONC::Constant:
    case SHAPE_CONC::Envelope:
    case SHAPE_CONC::EnvelopePhi:
    {
        for (size_t i = 0; i < x.size(); ++i)
        {
            u[i] = sqrt(g*std::abs(zb[i]));
        }
        break;
    }
    case SHAPE_CONC::NONE:
    {
        double umax = 9.01;
        for (size_t i = 0; i < x.size(); ++i)
        {
            if (x[i] < 0.4) {
                u[i] = x[i] / 0.4 * umax;
            }
            else if (x[i] >= 0.4 && x[i] < 0.6) {
                u[i] = umax;
            }
            else {
                u[i] = (1.0 - x[i]) / 0.4 * umax;
            }
            u[i] = 0.50;
        }
        break;
    }
    default:
    {
        break;
    }
    }
}
