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
//   Initial conditions for the 2D heat equation
//

#include "initial_conditions.h"

void initial_conditions(std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny,
    std::vector<double>& T_ini,
    std::vector<std::string> & ini_vars, double amp, double gauss_mu_x, double gauss_mu_y, double gauss_sigma_x, double gauss_sigma_y)
{
    // 
    // initialize T (temperature)
    //
    double T_giv = 0.0;
    size_t k;
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            // 
            // initialization via water level, u-velocity and v-velocity
            //
            k = p_idx(i, j, ny);
            if (ini_vars[0] == "zeta")
            {
                if (ini_vars[1] == "zeta_constant")
                {
                    T_giv = amp;  // initial water level
                }
                else if (ini_vars[1] == "zeta_linear_y")
                {
                    size_t k0 = p_idx(i, 0, ny);
                    size_t k1 = p_idx(j, ny - 1, ny);
                    T_giv = amp * (y[k] - y[0]) / (y[k1] - y[k0]); // initial water level
                }
                else if (ini_vars[1] == "zeta_linear_x")
                {
                    size_t k0 = p_idx(0     , 0, ny);
                    size_t k1 = p_idx(nx - 1, 0, ny);
                    T_giv = amp * (x[k] - x[0]) / (x[k1] - x[k0]); // initial water level
                }
                else if (ini_vars[1] == "zeta_gauss_hump")
                {
                    T_giv = amp * std::exp(
                        -(  (x[k] - gauss_mu_x) * (x[k] - gauss_mu_x) / (2. * gauss_sigma_x * gauss_sigma_x) +
                            (y[k] - gauss_mu_y) * (y[k] - gauss_mu_y) / (2. * gauss_sigma_y * gauss_sigma_y)
                            ) );  // initial water level
                }
                else if (ini_vars[1] == "zeta_gauss_hump_x")
                {
                    T_giv = amp * std::exp( -(x[k] - gauss_mu_x) * (x[k] - gauss_mu_x) / (2. * gauss_sigma_x * gauss_sigma_x) );  // initial water level
                }
                else if (ini_vars[1] == "zeta_gauss_hump_y")
                {
                    T_giv = amp * std::exp( -(y[k] - gauss_mu_y) * (y[k] - gauss_mu_y) / (2. * gauss_sigma_y * gauss_sigma_y) );  // initial water level
                }
                else
                {
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "Initialization variables." << std::endl;
                    std::cout << "Option: " << ini_vars[1] << "\' is not supported to intialize the water level ." << std::endl;
                    //std::cout << "Press Enter to finish";
                    //std::cin.ignore();
                    std::chrono::duration<int, std::milli> timespan(3000);
                    std::this_thread::sleep_for(timespan);
                    exit(1);
                }
            }
            else
            {
                std::cout << "----------------------------" << std::endl;
                std::cout << "Initialization variables." << std::endl;
                std::cout << "Option: " << ini_vars[0] << "\' is not supported." << std::endl;
                //std::cout << "Press Enter to finish";
                //std::cin.ignore();
                std::chrono::duration<int, std::milli> timespan(3000);
                std::this_thread::sleep_for(timespan);
                exit(1);
            }
            T_ini[k] = T_giv;
        }
    }
}
inline size_t p_idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}
