//
// Programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
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

#ifndef __BEDSHEARSTRESS_H__
#define __BEDSHEARSTRESS_H__

#include <vector>
#include <Eigen/Sparse>

int bed_shear_stress_matrix_rhs(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double cf, double theta, double dx, double dy, int nx, int ny);
void bed_shear_stress_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double cf, int nx, int ny);

    inline int bed_shear_stress_p_index(int i, int j, int nx);
    inline double bed_shear_stress_scv(double& c0, double c1, double c2, double c3);
    inline double bed_shear_stress_J_10(double& h, double& q, double& r, double& cf);  // right hand side q-equation
    inline double bed_shear_stress_J_11(double& h, double& q, double& r, double& cf);
    inline double bed_shear_stress_J_12(double& h, double& q, double& r, double& cf);
    inline double bed_shear_stress_J_13(double& h, double& q, double& r, double& cf);
    inline double bed_shear_stress_J_20(double& h, double& q, double& r, double& cf);  // right hand side r-equation
    inline double bed_shear_stress_J_21(double& h, double& q, double& r, double& cf);
    inline double bed_shear_stress_J_22(double& h, double& q, double& r, double& cf);
    inline double bed_shear_stress_J_23(double& h, double& q, double& r, double& cf);
    inline double vecq(double& qp, double& rp);

#endif  // __BEDSHEARSTRESS_H__
