//---------------------------------------------------------------
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
//    Solving the 1D advection/diffusion equation, fully implicit with delta-formuation and Modified Newton iteration 
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
//   Initial concentration for the Advection-Diffusion equation
//

#include "compatible_function.h"
#include "adv_diff_init_concentration.h"

void adv_diff_init_concentration(std::vector<double>& mass, std::vector<double> & x, double Lx, SHAPE_CONC shape, std::vector<double>& u_out)
{
    double L_envelope;
    double fcent;
    double shift;
    size_t refine = 100;

    size_t nx = x.size();
    std::vector<double> x_ana(refine * (nx - 1) + 1, 0.0);
    std::vector<double> u_ana(refine * (nx - 1) + 1, 0.0);
    std::vector<double> cv(nx, 0.0);

    double dx = x[1] - x[0];
    
    for (size_t i = 0; i < refine * (nx - 1) + 1; ++i)
    {
        x_ana[i] = double(i) * dx / double(refine) - dx;
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
        shift = 0.125 * Lx;
        fcent = 1. / 2. * L_envelope + shift;
        for (int i = 1; i < refine * (nx - 1) + 1; ++i)
        {
            if (x_ana[i] > 0 + shift and x_ana[i] < Lx / 4. + shift)
            {
                u_ana[i] = (0.5 + 0.5 * cos(2.0 * M_PI * 5. * (x_ana[i] - fcent) / L_envelope)) *
                    (0.5 + 0.5 * cos(2.0 * M_PI * 1. * (x_ana[i] - fcent) / L_envelope));
            }
        }
        (void) control_volumes(u_ana, cv, dx, refine);
        (void) compatible_function(mass, cv, u_ana, u_out, dx, refine);
        break;
    case SHAPE_CONC::EnvelopePhi:
        // Special function, as supplied by Mart Borsboom
        L_envelope = 0.25 * Lx;
        shift = 0.125 * Lx;
        fcent = 1. / 2. * L_envelope + shift;

        L_envelope = 0.25 * Lx;
        shift = 0.125 * Lx;
        fcent = 1. / 2. * L_envelope + shift;
        for (int i = 1; i < refine * (nx - 1) + 1; ++i)
        {
            if (x_ana[i] > 0 + shift and x_ana[i] < Lx / 4. + shift)
            {
                u_ana[i] = (0.5 + 0.5 * cos(2.0 * M_PI * 5. * (x_ana[i] - fcent) / L_envelope)) *
                    (0.5 + 0.5 * cos(2.0 * M_PI * 1. * (x_ana[i] - fcent) / L_envelope));
            }
        }
        (void)control_volumes(u_ana, cv, dx, refine);
        (void)compatible_function(mass, cv, u_ana, u_out, dx, refine);

        for (int i = 0; i < nx; ++i)
        {
            if (u_out[i] < 1.e-12)
            {
                u_out[i] = -25.;  // $ln(1e-12) \pm -27.0$
            }
            else
            {
                u_out[i] = std::log(u_out[i]);
            }
        }

        break;
    default:
        break;
    }
}
void control_volumes(std::vector<double>& u_ana, std::vector<double>& cv, double dx, size_t refine)
{
    size_t k;
    size_t nx = cv.size();
    for (size_t i = 1; i < nx - 1; ++i)
    {
        cv[i] = 0.;
        for (size_t j = 0; j < refine; ++j)
        {
            k = i * refine + j - size_t(refine / 2);
            cv[i] += dx / double(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
        }
    }
    size_t i = 0;
    cv[i] = 0.;
    for (size_t j = 0; j < (size_t)refine / 2; ++j)
    {
        k = i * refine + j;
        cv[i] += dx / double(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
    }
    i = nx - 1;
    cv[i] = 0.;
    for (int j = 0; j < (size_t)refine / 2; ++j)
    {
        k = i * refine + j - size_t(refine / 2);
        cv[i] += dx / double(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
    }

    return;
}
void compatible_function(std::vector<double>& mass, std::vector<double>& cv, std::vector<double>& u_ana, std::vector<double>& u_out, 
    double dx, size_t refine)
{
    size_t nx = cv.size();
    size_t nx_ana = u_ana.size();

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);               // solution vector 
    Eigen::VectorXd rhs(nx);                // RHS vector

    for (size_t i = 1; i < nx - 1; ++i)
    {
        A.coeffRef(i, i - 1) = dx * mass[0];
        A.coeffRef(i, i) = dx * mass[1];
        A.coeffRef(i, i + 1) = dx * mass[2];
        rhs[i] = cv[i];
    }
    size_t i = 0;
    A.coeffRef(i, i) = 1. / 12.;
    A.coeffRef(i, i + 1) = 10. / 12.;
    A.coeffRef(i, i + 2) = 1. / 12.;
    rhs[i] = u_ana[refine];
    i = nx - 1;
    A.coeffRef(i, i - 2) = 1. / 12.;
    A.coeffRef(i, i - 1) = 10. / 12.;
    A.coeffRef(i, i) = 1. / 12.;
    rhs[i] = u_ana[nx_ana - refine - 1];

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);
    for (size_t i = 0; i < nx; ++i)
    {
        u_out[i] = solution[i];
    }
}

