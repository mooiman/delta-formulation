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

void adv_diff_init_concentration(double time, double u_const, double eps, std::vector<double> & x, double Lx, SHAPE_CONC shape, std::vector<double>& u_out)
{
    double L_envelope;
    double fcent;
    double shift;
    size_t refine = 1;

    size_t nx = x.size();
    std::vector<double> x_ana(refine * (nx - 1) + 1, 0.0);
    std::vector<double> u_ana(refine * (nx - 1) + 1, 0.0);
    std::vector<double> cv(nx, 0.0);

    double dx = x[1] - x[0];
    
    for (size_t i = 0; i < refine * (nx - 1) + 1; ++i)
    {
        x_ana[i] = (double(i) * dx / double(refine) - dx);
    }
    switch (shape)
    {
    case SHAPE_CONC::Constant:
        for (size_t i = 0; i < nx; ++i)
        {
            u_out[i] = 0.0;
        }
        break;
    case SHAPE_CONC::Envelope:
        // Special function, as supplied by Mart Borsboom
        L_envelope = 0.25 * Lx;
        shift = 0.125 * Lx;
        fcent = 1. / 2. * L_envelope + shift;

        L_envelope = 0.25 * Lx;
        shift = 0.125 * Lx + u_const*time;
        fcent = 1. / 2. * L_envelope + shift;
        u_out[0] = 0.0;
        for (int i = 1; i < refine * (nx - 1) + 1; ++i)
        {
            u_out[i] = 0.0;
            if (x_ana[i] > 0 + shift and x_ana[i] < Lx / 4. + shift)
            {
                u_out[i] = (0.5 + 0.5 * cos(2.0 * M_PI * 5. * (x_ana[i] - fcent) / L_envelope)) *
                    (0.5 + 0.5 * cos(2.0 * M_PI * 1. * (x_ana[i] - fcent) / L_envelope));
            }
        }
        break;
    case SHAPE_CONC::Boundary_layer:
        {
            // Analytic solution outflow boundary layer
            double Pe = u_const * Lx/ eps;
            double a = 0.0;  // boundary value west
            double b = 1.0;  // boundary value east

            u_out[0] = 0.0;
            for (int i = 1; i < refine * (nx - 1); ++i)
            {
                u_out[i] = a + (b - a) *( std::exp(((x_ana[i] - Lx)/Lx) * Pe) - std::exp(-Pe) ) / (1.0 - std::exp(-Pe)); //eq 4.16 Wesseling
            }
            u_out[refine * (nx - 1)] = 2. * u_out[refine * (nx - 1) - 1] - u_out[refine * (nx - 1) - 2];
        }
        break;
    default:
        break;
    }
}
