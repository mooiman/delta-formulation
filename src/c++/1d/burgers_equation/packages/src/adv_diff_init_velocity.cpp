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
//   Initial concentration for the Advection-Diffusion equation
//
//   L_h C = d/dx (u.C - eps.d(C)/dx)
//

#include "adv_diff_init_velocity.h"

int adv_diff_init_velocity(std::vector<double>& u, const double u_initial, const std::vector<double>& x, std::string ini_var)
{
    int status = 1;
    if (ini_var == "constant")
    {
        for (size_t i = 0; i < x.size(); ++i)
        {
            u[i] = u_initial;
        }
        status = 0;
    }
    else if (ini_var == "colombo")
    {
        // u(x,0) = \exp{x} + 0.3 \exp( -200 (x + 0.5)^2)
        for (size_t i = 0; i < x.size(); ++i)
        {
            double mu = -0.5;
            double sigma = 0.05;
            u[i] = std::exp(x[i]) + 0.3 * std::exp( - (x[i] - mu) * (x[i] - mu) / (2. * sigma * sigma) );
        }
        status = 0;
    }
    else
    {
        status = 1;
    }
    return status;
}
