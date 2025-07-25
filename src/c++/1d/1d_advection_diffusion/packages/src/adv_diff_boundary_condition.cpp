//---------------------------------------------------------------
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
//   Boundary conditions for the Advection-Diffusion equation
//

#include "adv_diff_init_concentration.h"
#include "adv_diff_boundary_condition.h"

void adv_diff_boundary_condition(double& bc0_out, double& bc1_out, double& bc0_in, double& bc1_in, double& time, double& treg, std::string bc_signal)
{
    double reg_a = 0.0;
    double reg_b = 1.0;
    double reg_factor = 1.0;
    double reg_interp = 0.0;
    if (time < treg)
    {
        reg_factor = 0.5 * (std::cos(M_PI * (treg - time) / treg) + 1.0);
    }
    reg_interp = reg_a + (reg_b - reg_a) * reg_factor;  // 0 <= reg_factor <= 1

    if (bc_signal == "constant")
    {
        //
        // Given value at both sides
        //
        bc0_out = reg_interp * bc0_in;
        bc1_out = reg_interp * bc1_in;
    }
    else if (bc_signal == "sine")
    {
        //
        // given sine function at left boundary
        //
        if (time < treg)
        {
            bc0_out = reg_interp;
        }
        else
        {
            bc0_out = -std::cos(M_PI * (time) / treg);
        }
        bc1_out = reg_interp * bc1_in;
    }
    else
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Boundary signal not supported" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
        exit(1);
    }
}
