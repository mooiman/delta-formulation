//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
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
//   Initial conditions for the 2D wave equation
//

#include "initial_conditions.h"

void initial_conditions(std::vector<double>& x, size_t nx,
    std::vector<double>& s, std::vector<double>& u,
    std::vector<std::string> & ini_vars, double amp, double gauss_mu_x, double gauss_sigma_x)
{
    // 
    //initialize h, q and r
    //
    double s_giv = 0.0;
    double u_giv = 0.0;
    size_t k;
    for (size_t var_i = 0; var_i < ini_vars.size(); var_i ++)
    {
        for (size_t i = 0; i < nx; i++)
        {
            // 
            // initialization via water level, u-velocity and v-velocity
            //
            k = i;
            if (ini_vars[var_i] == "zeta_constant")
            {
                s_giv = amp;  // initial water level
            }
            else if (ini_vars[var_i] == "zeta_linear_x")
            {
                size_t k0 = 0;
                size_t k1 = nx - 1;
                s_giv = amp * (x[k] - x[0]) / (x[k1] - x[k0]); // initial water level
            }
            else if (ini_vars[var_i] == "zeta_gauss_hump")
            {
                s_giv = amp * std::exp( -(x[k] - gauss_mu_x) * (x[k] - gauss_mu_x) / (2. * gauss_sigma_x * gauss_sigma_x) );  // initial water level
            }
            else if (ini_vars[var_i] == "zeta_gauss_hump_x")
            {
                s_giv = amp * std::exp( -(x[k] - gauss_mu_x) * (x[k] - gauss_mu_x) / (2. * gauss_sigma_x * gauss_sigma_x) );  // initial water level
            }
            else if (ini_vars[var_i] == "u_constant")
            {
                u_giv = amp;  // initial water level
            }
            else
            {
                std::cout << "----------------------------" << std::endl;
                std::cout << "Initialization variables." << std::endl;
                std::cout << "Option: \"" << ini_vars[0] << "\" is not supported." << std::endl;
                //std::cout << "Press Enter to finish";
                //std::cin.ignore();
                std::chrono::duration<int, std::milli> timespan(3000);
                std::this_thread::sleep_for(timespan);
                exit(1);
            }
            s[k] = s_giv;
            u[k] = u_giv;  // Initial q=hu -velocity
        }
    }
}
