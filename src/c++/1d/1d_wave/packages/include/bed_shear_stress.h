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

    int bed_stress_matrix_rhs(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
        std::vector<double>& hp, std::vector<double>& qp,
        std::vector<double>& hn, std::vector<double>& qn,
        double cf, double dx, double theta);
    void bed_stress_rhs(std::vector<double>& rhs_q, std::vector<double> hn, std::vector<double>& qn);

    inline double bed_stress_scv(double c0, double c1);
    inline double bed_stress_J_10(double& h, double& q, double cf);  // right hand side q-equation
    inline double bed_stress_J_11(double& h, double& q, double cf);
    inline double bed_stress_J_12(double& h, double& q, double cf);
    inline double bed_stress_vecq(double& qp);

#endif  // __BEDSHEARSTRESS_H__
