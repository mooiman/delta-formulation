//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formulation and Modified Newton iteration 
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
//   DESCRIPTION
//
//   Analytic solution of the 2D heat equation with a Dirac-delta function in centre of the domain
//   2D: theta = \frac{1}{4D\pi t} \exp\left(\frac{-r^2}{4Dt}\right), \qquad r^2 = x^2 + y^2
//
#include "analytic_solution.h"

void analytic_solution(double time, double conductivity, std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny, std::vector<double>& h_ana)
{
    size_t k;
    size_t k_centre;
    k_centre = nx/2 * ny + ny/2;
    if (time == 0.0)
    {
        h_ana[k_centre] = 1.0;
        return;
    }
    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            k = i * ny + j;
            double dx = x[k] - x[k_centre];
            double dy = y[k] - y[k_centre];
            double r2 = dx*dx + dy*dy;
            
            h_ana[k] = 1. / (4. * conductivity * M_PI * time) * std::exp(-r2/(4.0 * conductivity * time));
        }
    }
    return;
}

