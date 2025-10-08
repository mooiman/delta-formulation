//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formuation and Modified Newton iteration 
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
//------------------------------------------------------------------------------

#ifndef __2D_REGULARIZATION_H__
#define __2D_REGULARIZATION_H__

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <math.h>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

class REGULARIZATION
{
public:
    REGULARIZATION();
    REGULARIZATION(int iter_max, double g);

    void given_function(
        std::vector<double>& u_out, std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& eq8,
        std::vector<double>& u_giv,
        size_t nx, size_t ny, double dx, double dy, double c_psi, std::ofstream& log_file);

    void first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx);

private:
    std::unique_ptr<std::vector<double>> solve_eq7(size_t nx, size_t ny, double dx, double dy, 
        std::vector<double> psi_11, std::vector<double> psi_22, std::vector<double> u_giv, std::ofstream& log_file);
    std::unique_ptr<std::vector<double>> solve_eq8(size_t nx, size_t ny, double dx, double dy, double c_psi, std::vector<double> u0, 
        std::vector<double> u0_xixi, std::vector<double> u0_etaeta, std::ofstream& log_file);

    int m_iter_max;
    double m_g;

    double m_alpha;
    std::vector<double> m_mass;
    double eps_smooth;  // epsilon of regularization process
    double m_u0_is_smooth;  // is function u0 is smooth enough, 
    size_t p_index(size_t i, size_t j, size_t ny);
};

#endif __2D_REGULARIZATION_H__
