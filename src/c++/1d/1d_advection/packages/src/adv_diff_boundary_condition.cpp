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

#include "adv_diff_init_concentration.h"
#include "adv_diff_boundary_condition.h"

void adv_diff_boundary_condition(double& bc0_out, double& bc1_out, double& bc0_in, double& bc1_in, double& time, double& treg, BND_TYPE bnd_type)
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

    switch (bnd_type)
    {
    case BND_TYPE::Constant:
        //
        // Given value at both sides
        //
        bc0_out = reg_interp * bc0_in;
        bc1_out = reg_interp * bc1_in;
        break;
    case BND_TYPE::Sine:
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
        break;
    }
}
