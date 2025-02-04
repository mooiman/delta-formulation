//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       $date$
//   version:    $version$
//   copyright © 2025 Mooiman
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Righthand side 0f  the Advection-Diffusion equation
//
//   L_h C = d/dx (u.C - eps.d(C)/dx) =q
//

#include "adv_diff_source.h"

void adv_diff_source(std::vector<double> & q, int select)
{
    switch (select)
    {
    case 1:
        //
        // Homogeneous righthand side
        //
        for (int i = 0; i < q.size(); ++i)
        {
            q[i] = 0.0;
        }
        break;
    }
}
