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

#ifndef __MATRIX_ASSEMBLY_BOUNDARIES_H__
#define __MATRIX_ASSEMBLY_BOUNDARIES_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>


int boundary_north(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity,
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_NORTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess);

int boundary_east(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity,
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_EAST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess);

int boundary_south(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity,
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_SOUTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess);

int boundary_west(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity,
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_WEST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess);

// double natural_boundary_qp(std::vector<double>& hp, int p_b, int xi, int eta, std::vector<double>& w, int ny);

void molecule(std::vector<size_t>& p, size_t p_sw, size_t ny);

inline void add_value(double * values, size_t col, double data);
inline size_t ma_index(size_t i, size_t j, size_t ny_in);

#endif  // __MATRIX_ASSEMBLY_BOUNDARIES_H__
