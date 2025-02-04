//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       $date$
//   version:    $version$
//   copyright © 2025 Mooiman
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Initial concentration for the Advection-Diffusion equation
//
//   L_h C = d/dx (u.C - eps.d(C)/dx)
//

#include "adv_diff_init_velocity.h"

void adv_diff_init_velocity(std::vector<double>& u, const double u_in,  const std::vector<double>& x, int select)
{
    switch (select)
    {
        case 1:
        case 2:
        {
            for (int i = 0; i < x.size(); ++i)
            {
                u[i] = u_in;
             }
            break;
        }
        case 3:
        {
            double umax = 9.01;
            for (int i = 0; i < x.size(); ++i)
            {
                if (x[i] < 0.4) {
                    u[i] = x[i] / 0.4 * umax;
                }
                else if (x[i] >= 0.4 && x[i] < 0.6) {
                    u[i] = umax;
                }
                else {
                    u[i] = (1.0 - x[i]) / 0.4 * umax;
                }
                u[i] = 0.50;
            }
            break;
        }
        default:
        {
            break;
        }
    }
}
