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

#ifndef __CONVECTION_H__
#define __CONVECTION_H__

#include <vector>
#include <Eigen/Sparse>

#define UNUSED(x) (void)x;

int convection_matrix_rhs(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double theta, double dx, double dy, int nx, int ny);
int convection_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double dx, double dy, int nx, int ny );                // RHS vector [h, q, r]^{n}

inline double convection_J_10(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_20(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_11(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_21(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_12(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_22(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_13(double& h, double& q, double& r, double nxi, double neta);
inline double convection_J_23(double& h, double& q, double& r, double nxi, double neta);

inline double convection_scvf_xi (double c0, double c1, double c2, double c3);
inline double convection_scvf_eta(double c0, double c1, double c2, double c3);

inline int convection_p_index(int i, int j, int nx);
inline void set_value(double * values, int col, double data);

#endif  // __CONVECTION_H__
