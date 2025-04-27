//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2025-03-10
//   copyright © 2025 Mooiman
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
//   Initial function compatible with numerical scheme
//

#include "compatible_function.h"

void compatible_function(std::vector<double>& x, std::vector<double>& u_out, std::vector<double>& u_in, std::vector<double>& mass)
{
    int nx = (int) u_out.size();
    int refine = 4;
    double dx = x[1] - x[0];
    double Lx = x[nx - 2] - x[1];  // skip 2 virtual points
    
    std::vector<double> x_ana(refine*(nx-1) + 1, 0.0);  
    std::vector<double> u_ana(refine*(nx-1) + 1, 0.0);  
    std::vector<double> cv_ana(nx, 0.0);  
    
    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);               // solution vector 
    Eigen::VectorXd rhs(nx);                // RHS vector

    // compute the integral for the control volumes of fuction u_in.
    // integral of the control volumes (xi-1/2, xi+1/2)
    for (int i = 0; i < refine*(nx-1)+1; ++i)
    {
        x_ana[i] = double(i) * dx/double(refine) - x[0];
    }
    
    double L_envelope = 0.25 * Lx;
    double shift = 0.125 * Lx;
    double fcent = 1. / 2. * L_envelope + shift;
    for (int i = 1; i < refine*(nx-1)+1; ++i)
    {
        if (x_ana[i] > 0 + shift and x_ana[i] < Lx / 4. + shift)
        {
            u_ana[i] = (0.5 + 0.5 * cos(2.0 * M_PI * 5. * (x[i] - fcent) / L_envelope)) *
                       (0.5 + 0.5 * cos(2.0 * M_PI * 1. * (x[i] - fcent) / L_envelope));
        }
    }
    int k;
    for (int i = 1; i < nx-1; ++i)
    {
        cv_ana[i] = 0.;
        for (int j = 0; j < refine; ++j)
        {
            k = i * refine + j - int(refine / 2);
            cv_ana[i] += dx / double(refine) * 0.5 * (u_ana[k] + u_ana[k + 1]);
        }
    }
    int i = 0;
    cv_ana[i] = 0.;
    for (int j = 0; j < (int) refine/2; ++j)
    {
        k = i * refine + j;
        cv_ana[i] += dx/double(refine) * 0.5 * (u_ana[k] + u_ana[k+1]);
    }
    i = nx-1;
    cv_ana[i] = 0.;
    for (int j = 0; j < (int) refine/2; ++j)
    {
        k = i * refine + j - int(refine/2);
        cv_ana[i] += dx/double(refine) * 0.5 * (u_ana[k] + u_ana[k+1]);
    }

    for (int i = 0; i < nx; ++i)
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = 0.0;
    }
    for (int i = 1; i < nx-1; ++i)
    {
        A.coeffRef(i, i-1) = dx * mass[0];
        A.coeffRef(i, i  ) = dx * mass[1];
        A.coeffRef(i, i+1) = dx * mass[2];
        rhs[i] = cv_ana[i];
    }
    i = 0;
    A.coeffRef(i, i  ) = 1./12.;
    A.coeffRef(i, i+1) = 10./12.;
    A.coeffRef(i, i+2) = 1./12.;
    rhs[i] = u_ana[refine + 1];
    i = nx - 1;
    A.coeffRef(i, i-2) = 1./12.;
    A.coeffRef(i, i-1) = 10./12.;
    A.coeffRef(i, i  ) = 1./12.;
    rhs[i] = u_ana[i - refine];
    
    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);
    for (int i = 0; i < nx; ++i)
    {
        u_out[i] = solution[i];
    }
}