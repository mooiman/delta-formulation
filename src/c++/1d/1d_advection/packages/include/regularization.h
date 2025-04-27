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
#ifndef __ADV_DIFF_REGULARIZATION_H__
#define __ADV_DIFF_REGULARIZATION_H__

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

class REGULARIZATION
{
public:
    REGULARIZATION();
    REGULARIZATION(int);
    void old_version(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& err, std::vector<double>& u_giv,
        double dx, double c_psi, bool use_eq8);
    void given_function(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& eq8, std::vector<double>& u_giv_in,
        double dx, double c_psi, bool use_eq8);
    void first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx);
private:
    std::unique_ptr<std::vector<double>> solve_eq8(double c_psi, std::vector<double> u0, std::vector<double> u0_xixi);
    std::unique_ptr<std::vector<double>> solve_eq7(double dx, std::vector<double> psi, std::vector<double> u_giv);

    int m_iter_max;
};

#endif __ADV_DIFF_REGULARIZATION_H__
