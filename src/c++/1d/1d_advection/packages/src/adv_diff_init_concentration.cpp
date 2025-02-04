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

#include "adv_diff_init_concentration.h"

void adv_diff_init_concentration(std::vector<double>& d, int select)
{
    switch (select)
    {
    case 1 or 2:
        for (int i = 0; i < d.size(); ++i)
        {
            d[i] = 0.0;
        }
        d[2 * d.size() / 4] = 0.0;
        break;
    default:
        break;
    }
}
