//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       $date$
//   version:    $version$
//   copyright © 2025 Mooiman
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Boundary conditions for the Advection-Diffusion equation
//

#include "adv_diff_boundary_condition.h"

void adv_diff_boundary_condition(double& bc0, double& bc1, double& time, double& treg, int select)
{
    bc0 = 1.0;
    bc1 = 0.0;
    //
    double reg_a = 0.0;
    double reg_b = 1.0;
    double reg_factor = 1.0;
    double reg_interp = 0.0;
    if (time < treg)
    {
        reg_factor = 0.5 * (std::cos(M_PI * (treg - time) / treg) + 1.0);
    }
    reg_interp = reg_a + (reg_b - reg_a) * reg_factor;  // 0 <= reg_factor <= 1

    switch (select)
    {
    case 1:
        //
        // Given value at both sides
        //
        bc0 = reg_interp * bc0;
        bc1 = reg_interp * bc1;
        break;
    case 2:
        //
        // given sine function at left boundary
        //
        if (time < treg) 
        {
            bc0 = reg_interp;
        }
        else
        {
            bc0 = -std::cos(M_PI * (time) / treg);
        }
        bc1 = reg_interp * bc1;
        break;
    }
}
