//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2025-02-13
//   copyright © 2025 Mooiman
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Boundary conditions for the 1D wave equation
//

#include "boundary_condition.h"

void boundary_condition(double& bc0, double& bc1, double& bc0_in, double& bc1_in, double& time, double& treg, int select)
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

    switch (select)
    {
    case 1:
        //
        // Given value at both sides
        //
        bc0 = reg_interp * bc0_in;
        bc1 = reg_interp * bc1_in;
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
        bc1 = reg_interp * bc1_in;
        break;
    }
}
