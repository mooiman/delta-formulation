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
//   Regularization of given function
//

#include "adv_diff_regularization.h"

void adv_diff_regularization(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& u_giv, double dx, double c_psi)
{
    int nx = (int) u_giv.size();
    int iter_max = 100;
    double c_error;
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    std::vector<double> err(nx, 0.);
    std::vector<double> u0(nx, 0.);
    std::vector<double> u1(nx, 0.);
    std::vector<double> u0_xx(nx, 0.);
    std::vector<double> tmp(nx, 0.);

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::SparseMatrix<double> B(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    std::ofstream log_file;
#if defined(DEBUG)
    log_file.open("../output/janm.log");
#endif

    for (int i = 0; i < nx; ++i)
    {
        u0[i] = u_giv[i];
    }
    for (int i = 0; i < nx; ++i)
    {
        A.coeffRef(i, i) = 1.0;
        B.coeffRef(i, i) = 1.0;
        rhs[i] = 0.0;
    }
    for (int iter = 0; iter < iter_max; ++iter)
    {
        //std::cout << "Iteration: " << iter << std::endl;
        //std::cout << "Compute second derivative" << std::endl;
        double u0xx_max = 0.0;
        for (int i = 1; i < nx - 1; ++i)
        {
            u0_xx[i] = (u0[i - 1] - 2. * u0[i] + u0[i + 1]);
            u0xx_max = std::max(u0xx_max, std::abs(u0_xx[i]));
        }
        int i = 0;
        u0_xx[i] = (u0[i] - 2. * u0[i] + u0[i + 1]);
        //u0_xx[i] = u0_xx[i + 1];
        i = nx - 1;
        u0_xx[i] = (u0[i - 1] - 2. * u0[i] + u0[i]);
        //u0_xx[i] = u0_xx[i - 1];

        //std::cout << "Initialization matrix A and rhs" << std::endl;
        c_error = c_psi;

        for (int i = 1; i < nx - 1; ++i)
        {
            A.coeffRef(i, i - 1) = 1. / 8. - c_error;
            A.coeffRef(i, i) = 6. / 8. + 2. * c_error;
            A.coeffRef(i, i + 1) = 1. / 8. - c_error;
        }
        //std::cout << std::endl;
        for (int i = 1; i < nx - 1; ++i)
        {
            tmp[i] = std::abs(u0_xx[i]);
            rhs[i] = tmp[i];
        }
        i = 0;
        A.coeffRef(i, i) = 1./12.;
        A.coeffRef(i, i + 1) = 10. / 12.;
        A.coeffRef(i, i + 2) = 1. / 12.;
        rhs[i] = rhs[i+1];

        i = nx - 1;
        A.coeffRef(i, i - 2) = 1. / 12.;
        A.coeffRef(i, i - 1) = 10./12.;
        A.coeffRef(i, i) = 1. / 12.;
        rhs[i] = rhs[i-1];

#if defined(DEBUG)
        if (false)
        {
            log_file << "=== Matrix A ==========================================" << std::endl;
            log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A) << std::endl;
            log_file << "=== RHS A =============================================" << std::endl;
            log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
        }
#endif

        //std::cout << "Copy previous iteration value u0 to solution vector for initial guess" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            solution[i] = u0[i];
        }

        //std::cout << "Initialize matrix for BiCGStab" << std::endl;
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
        //std::cout << "Set solver" << std::endl;
        solverA.compute(A);
        //std::cout << "Solve matrix A with BiCGStab" << std::endl;
        solution = solverA.solve(rhs);
        //solution = solverA.solveWithGuess(rhs, solution);

        double err_max = *std::max_element(solution.begin(), solution.end());
        // Compute artifical viscosity, psi
        //std::cout << "Compute artifical viscosity" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            err[i] = solution[i];
            psi[i] = c_psi * dx * dx  * err[i] / (err_max + 1e-13);  // to prevent division bij zero
        }

        //std::cout << "Defining matrix B and rhs" << std::endl;

        for (int i = 1; i < nx - 1; ++i)
        {
            double psi_im12 = 0.5 * (psi[i] + psi[i - 1]);
            double psi_ip12 = 0.5 * (psi[i + 1] + psi[i]);
            psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
            psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

            B.coeffRef(i, i - 1) = dx * 1./8. - psi_im12 / dx;
            B.coeffRef(i, i) = dx * 6./8. + psi_ip12 / dx + psi_im12 / dx;
            B.coeffRef(i, i + 1) = dx * 1./8. - psi_ip12 / dx;
        }
        for (int i = 1; i < nx - 1; ++i)
        {
            rhs[i] = dx * (1. / 8. * u_giv[i - 1] + 6. / 8. * u_giv[i] + 1. / 8. * u_giv[i + 1]);
        }
        //std::cout << std::endl;
        i = 0;
        B.coeffRef(i, i)     = 1. / 12.;  11. / 24.;
        B.coeffRef(i, i + 1) = 10. / 12.;  14. / 24.;
        B.coeffRef(i, i + 2) = 1. / 12.;  -1. / 24.;
        rhs[i] = u_giv[i+1];

        i = nx - 1;
        B.coeffRef(i, i - 2) = 1./12.;  -1. / 24.;
        B.coeffRef(i, i - 1) = 10. / 12.;  14. / 24.;
        B.coeffRef(i, i)     = 1. / 12.;  11. / 24.;
        rhs[i] = u_giv[i-1];

#if defined(DEBUG)
        if (false)
        {
            log_file << "=== Matrix B ==========================================" << std::endl;
            log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(B) << std::endl;
            log_file << "=== RHS B =============================================" << std::endl;
            log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
        }
#endif
        //std::cout << "Copy previous iteration value u0 to solution vector for initial guess" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            solution[i] = u0[i];
        }

        //std::cout << "Initialize matrix for BiCGStab" << std::endl;
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverB;
        //std::cout << "Set solver" << std::endl;
        solverB.compute(B);
        //std::cout << "Solve matrix B with BiCGStab" << std::endl;
        //solution = solverB.solve(rhs);
        solution = solverB.solveWithGuess(rhs, solution);
        //std::cout << "Copy solution to u1-vector" << std::endl;

        diff_max1 = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(solution[i] - u_giv[i]));
            u0[i] = solution[i];
        }
        if (std::abs(diff_max1 - diff_max0) > 1e-12)
        {
            //log_file << std::setprecision(8) << std::scientific << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            diff_max0 = diff_max1;
        }
        else
        {
            //log_file << std::setprecision(8) << std::scientific  << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            break;
        }
        //log_file << std::endl;
    }
#if defined(DEBUG)
    log_file.close();
#endif    
    for (int i = 0; i < nx; ++i)
    {
        u_out[i] = u0[i];
    }
}
