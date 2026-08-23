//
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
//   Initial concentration for the Burgers equation

#include "adv_diff_init_velocity.h"

int adv_diff_init_velocity(std::vector<double>& u, const std::vector<double>& ini_u_bnd, const std::vector<double>& x, std::string ini_var)
{
    int status = 1;
    size_t nx = x.size();

    bool compatible = false;
    size_t refine = 2;
    size_t refine_bnd = refine/2;  // boundary at controlvolumen edge
    refine_bnd = refine;  // boundary at node
    double Lx = x[nx - 2] - x[1];
    double dx = Lx / double(nx - 3);  

    std::vector<double> x_ana(refine * (nx - 1) + 1, 0.0);
    std::vector<double> u_ana(refine * (nx - 1) + 1, 0.0);

    for (size_t i = 0; i < refine*(nx - 1)+1; ++i)
    {
        x_ana[i] = double(i) / (double(refine * (nx - 3))) * Lx - dx;
    }
    if (ini_var == "constant")
    {
       compatible = true;  // fill u[i] with constant value ini_u_bnd[0]
       for (size_t i = 0; i < nx; ++i)
        {
            u[i] = ini_u_bnd[0];
        }
        status = 0;
    }
    else if (ini_var == "smoothstep")
    {
        compatible = false;
        u[0] = ini_u_bnd[0];
        for (size_t i = 1; i < nx - 1; ++i)
        {
            double xtmp= x[i]/Lx;
            double phi_x = std::exp(-1. / xtmp);
            u[i] = ini_u_bnd[0] + (ini_u_bnd[1] - ini_u_bnd[0]) * phi_x / (phi_x + std::exp(-1. / (1. - xtmp)));      
        }
        u[nx - 1] = ini_u_bnd[1];

        // smooth stepfunction between boundary values ini_u_bnd[0] and ini_u_bnd[1] at the left and right boundary of the domain x
        for (size_t i = 0; i < refine_bnd; ++i)
        {
            u_ana[i] = ini_u_bnd[0];
        }
        for (size_t i = refine_bnd; i < refine * (nx - 1) - refine_bnd + 1; ++i)
        {
            double xtmp= x_ana[i]/Lx;
            double phi_x = std::exp(-1. / xtmp);
            u_ana[i] = ini_u_bnd[0] + (ini_u_bnd[1] - ini_u_bnd[0]) * phi_x / (phi_x + std::exp(-1. / (1. - xtmp)));      
        }
        for (size_t i = refine * (nx - 1) - refine_bnd + 1; i < refine * (nx - 1) + 1; ++i)
        {
            u_ana[i] = ini_u_bnd[1];
        }
        //
        // u_init = a - (a+b)*phi(1/2+(x - x_c)/alpha) / ( phi(1/2+(x - x_c)/alpha) + phi(1/2 - (x - x_c)/alpha) )
        // 
        // Met phi de bekende functie phi(x) = exp(-1/x), x > 0; = 0, x <= 0, met a, b > 0, 
        // met alpha (in meters) de parameter waarmee je de steilheid van de overgang tussen a en -b regelt 
        // (typisch (een fractie van) de lengte van het domein), en met x_c de positie van die overgang.
    }
    else if (ini_var == "linear")
    {
        compatible = true;  // fill u[i] with linear function between boundary values ini_u_bnd[0] and ini_u_bnd[1]
        size_t shift = 1;
        for (size_t i = 0; i < shift; ++i)
        {
           u[i] = ini_u_bnd[0];
        }
        for (size_t i = shift; i < nx - shift; ++i)
        {
            double x_begin = x[shift];
            double x_end = x[nx - 1 - shift];
            double alpha = (x[i] - x_begin) / (x_end - x_begin);
            u[i] = (1. - alpha) * ini_u_bnd[0] + alpha * ini_u_bnd[1];
        }
        for (size_t i = nx - shift; i < nx; ++i)
        {
           u[i] = ini_u_bnd[1];
        }
        status = 0;
    }
    else if (ini_var == "colombo")
    {
        compatible = false;
        double mu = -0.5;
        double sigma = 0.05;
        double xtmp;
        // u(x,0) = \exp{x} + 0.3 \exp( -200 (x + 0.5)^2)
        for (size_t i = 0; i < nx; ++i)
        {
            u[i] = std::exp(x[i]) + 0.3 * std::exp( - (x[i] - mu) * (x[i] - mu) / (2. * sigma * sigma) );
        }

        for (size_t i = 0; i < refine_bnd; ++i)
        {
            xtmp = x[1]/Lx;  // boundary location at left side
            u_ana[i] = std::exp(xtmp) + 0.3 * std::exp( - (xtmp - mu/Lx) * (xtmp - mu/Lx) / (2. * sigma * sigma) );
        }
        for (size_t i = refine_bnd; i < refine * (nx - 1) - refine_bnd + 1; ++i)
        {
            double xtmp= x_ana[i]/Lx;
            u_ana[i] = std::exp(xtmp) + 0.3 * std::exp( - (xtmp - mu/Lx) * (xtmp - mu/Lx) / (2. * sigma * sigma) );
        }
        for (size_t i = refine * (nx - 1) - refine_bnd + 1; i < refine * (nx - 1) + 1; ++i)
        {
            xtmp = x[nx - 1]/Lx;  // boundary location at right side
            u_ana[i] = std::exp(xtmp) + 0.3 * std::exp( - (xtmp - mu/Lx) * (xtmp - mu/Lx) / (2. * sigma * sigma) );
        }
        status = 0;
    }
    else
    {
        status = 1;
    }
    //if (compatible)
    if (true)
    {
        return status;
    }
    //
    // Make the fuctions compatible with the numerical scheme
    //
    std::vector<double> cv_coarse(nx, 0.0);
    std::vector<double> cv_compatible(nx, 0.0);
    std::vector<double> cv_ana(nx, 0.0);

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector 
    Eigen::VectorXd rhs(nx);                // RHS vector
    // 
    // Integral over the control volumes (xi-1/2, xi+1/2)
    // Coarse solution
    for (size_t i = 1; i < nx - 1; ++i)
    {
        cv_coarse[i] = dx / 2.0 * 0.25 * (3.0 * u[i] + u[i - 1]) + dx / 2.0 * 0.25 * (3.0 * u[i] + u[i + 1]);
    }

    size_t i = 0;
    cv_coarse[i] = dx / 2.0 * 0.25 * (3.0 * u[i] + u[i + 1]);

    i = nx - 1;
    cv_coarse[i] = dx / 2.0 * 0.25 * (3.0 * u[i] + u[i - 1]);

    // Analytical solution
    for (size_t i = 1; i < nx - 1; ++i)
    {
        cv_ana[i] = 0.0;

        for (size_t j = 0; j < refine; ++j)
        {
            size_t k = i * refine + j - static_cast<int>(refine / 2);

            cv_ana[i] += dx / static_cast<double>(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
        }
    }

    i = 0;
    cv_ana[i] = 0.0;
    for (size_t j = 0; j < static_cast<int>(refine / 2); ++j)
    {
        size_t k = i * refine + j;

        cv_ana[i] += dx / static_cast<double>(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
    }
    i = nx - 1;
    cv_ana[i] = 0.0;

    for (size_t j = 0; j < static_cast<int>(refine / 2); ++j)
    {
        size_t k = i * refine + j - static_cast<int>(refine / 2);

        cv_ana[i] += dx / static_cast<double>(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
    }

    std::vector<double> m_mass;
    double m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);

    i = 0;
    A.coeffRef(i, i    ) = -1.0; 
    A.coeffRef(i, i + 1) =  1.0;
    A.coeffRef(i, i + 2) =  0.0;
    rhs[i] = 0.0;
    for (size_t i = 1; i < nx - 1; ++i)
    {
        A.coeffRef(i, i - 1) = dx * m_mass[0];
        A.coeffRef(i, i    ) = dx * m_mass[1];
        A.coeffRef(i, i + 1) = dx * m_mass[2];

        rhs[i] = cv_ana[i];
    }
    i = nx - 1;
    A.coeffRef(i, i - 2) =  0.0;
    A.coeffRef(i, i - 1) = -1.0;
    A.coeffRef(i, i    ) =  1.0;
    rhs[i] =  0.0;

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);

    for (size_t i = 0; i < nx; ++i)
    {
        u[i] = solution[i];
    }

    for (size_t i = 1; i < nx - 1; ++i)
        cv_compatible[i] =  dx/2. * 0.25 * (3. * u[i] +u[i - 1]) + dx/2. * 0.25 * (3. * u[i] + u[i + 1]);
    i = 0;
    cv_compatible[i] = dx/2. * 0.25 * (3. * u[i] + u[i + 1]);
    i = nx-1;
    cv_compatible[i] = dx/2. * 0.25 * (3. * u[i] + u[i - 1]);

    return status;
}
