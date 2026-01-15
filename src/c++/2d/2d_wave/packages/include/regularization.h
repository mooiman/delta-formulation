//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formuation and Modified Newton iteration 
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

#include "matrix_assembly_psi_boundaries.h"
#include "matrix_assembly_psi_corners.h"
#include "matrix_assembly_psi_interior.h"
#include "matrix_assembly_utilde_boundaries.h"
#include "matrix_assembly_utilde_corners.h"
#include "matrix_assembly_utilde_interior.h"
#include "diffusion2d.h"
#include "print_matrix.h"

class REGULARIZATION
{
public:
    REGULARIZATION();
    REGULARIZATION(int iter_max, double g, std::string logging);

    void given_function(
        std::vector<double>& u_out, std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& eq8,
        std::vector<double>& x, std::vector<double>& y,
        std::vector<double>& u_giv,
        size_t nx, size_t ny, double c_psi, std::ofstream& log_file);

    void artificial_viscosity(std::vector<double>& psi, 
        std::vector<double>& h, std::vector<double>& q, std::vector<double>& r, std::vector<double>& zb, 
        std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny, 
        double c_psi_in, std::ofstream& log_file);

    void first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx);

private:
    std::unique_ptr<std::vector<double>> solve_eq7(size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y, 
        std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& u_giv, std::ofstream& log_file);
    std::unique_ptr<std::vector<double>> solve_eq8(size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y, 
        double c_psi, std::vector<double>& u0, std::vector<double>& u0_xixi, std::vector<double>& u0_etaeta, 
        std::ofstream& log_file);
    int reg_interior_rhs_psi( size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
        std::vector<double>& h, std::vector<double>& q, std::vector<double>& r,
        std::vector<double>& x, std::vector<double>& y,
        double c_psi, double g, size_t nx, size_t ny);

    inline double F1(std::vector<double> & u, std::vector<size_t>& p, std::vector<double> & x, std::vector<double> &y, size_t nx, size_t ny );
    inline double F2(std::vector<double> & u, std::vector<size_t>& p, std::vector<double> & x, std::vector<double> &y, size_t nx, size_t ny );
    inline double F3(std::vector<double> & u, std::vector<size_t>& p, std::vector<double> & x, std::vector<double> &y, size_t nx, size_t ny );

    inline double d2udxi2(std::vector<double> & u, std::vector<size_t>& p);
    inline double d2udxideta(std::vector<double> & u, std::vector<size_t>& p);
    inline double d2udeta2(std::vector<double> & u, std::vector<size_t>& p);

    int m_iter_max;
    double m_g;
    std::string m_logging;
    std::ofstream m_log_file;

    double m_alpha;
    std::vector<double> m_mass;
    double m_eps_smooth;  // epsilon of regularization process
    double m_u0_is_smooth;  // is function u0 is smooth enough, 
    size_t p_index(size_t i, size_t j, size_t ny);
};

#endif __2D_REGULARIZATION_H__
