//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2025-02-19
//   copyright © 2025 Mooiman
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
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Boundary conditions for the 1D wave equation
//

#include "initial_conditions.h"

void initial_conditions(std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp, std::vector<double>& zb_giv,
    std::vector<std::string> & ini_vars, double amp, double gauss_mu, double gauss_sigma)
{
    // 
    //initialize h, q and r
    //
    double s_giv = 0.0;
    double u_giv = 0.0;
    double v_giv = 0.0;
    size_t k;
    for (size_t j = 0; j < ny; j++)
    {
        for (size_t i = 0; i < nx; i++)
        {
            // 
            // initialize via water level, u-velocity and v-velocity
            //
            k = j * nx + i;
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_constant")
            {
                s_giv = amp;  // initial water level
            }
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_linear_y")
            {
                size_t k0 = 0 * nx + i;
                size_t k1 = (ny - 1) * nx + i;
                s_giv = amp * (y[k] - y[0]) / (y[k1] - y[k0]); // initial water level
            }
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_linear_x")
            {
                size_t k0 = j * nx + 0;
                size_t k1 = j * nx + nx - 1;
                s_giv = amp * (x[k] - x[0]) / (x[k1] - x[k0]); // initial water level
            }
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_GaussHump")
            {
                s_giv = amp * std::exp(-((x[k] - gauss_mu) * (x[k] - gauss_mu) + (y[k] - gauss_mu) * (y[k] - gauss_mu)) / (2. * gauss_sigma * gauss_sigma));  // initial water level
            }
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_GaussHump_x")
            {
                s_giv = amp * std::exp(-((x[k] - gauss_mu) * (x[k] - gauss_mu)) / (2. * gauss_sigma * gauss_sigma));  // initial water level
            }
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_GaussHump_y")
            {
                s_giv = amp * std::exp(-((y[k] - gauss_mu) * (y[k] - gauss_mu)) / (2. * gauss_sigma * gauss_sigma));  // initial water level
            }
            u_giv = 0.0;
            v_giv = 0.0;
            hn[k] = s_giv - zb_giv[k];  // Initial water depth
            qn[k] = hn[k] * u_giv;  // Initial q=hu -velocity
            rn[k] = hn[k] * v_giv;  // Initial r=hv -velocity
            // 
            hp[k] = hn[k]; 
            qp[k] = qn[k]; 
            rp[k] = rn[k]; 
        }
    }
}
