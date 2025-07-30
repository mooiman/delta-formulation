//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
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

#include "convection.h"

//   nw - - - - - - - n - - - - - - - ne        nw - - - - - - - n - - - - - - - ne
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    | - - - - - x - | - x - - - - - |          | - - - - - 6 - | - 5 - - - - - |
//    |       |       |       |       |          |       |       |       |       |
//    |       x       |       x       |          |       7 scv_3 | scv_2 4       |
//    |       |       |       |       |          |       |       |       |       |
//    w - - - - - - - 0 - - - - - - - e          w - - - - - - - 0 - - - - - - - e
//    |       |       |       |       |          |       |       |       |       |
//    |       x       |       x       |          |       0 scv_0 | scv_1 3       |
//    |       |       |       |       |          |       |       |       |       |
//    | - - - - - x - | - x - - - - - |          | - - - - - 1 - | - 2 - - - - - |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//    |       |       |       |       |          |       |       |       |       |
//   sw - - - - - - - s - - - - - - - se        sw - - - - - - - s - - - - - - - se

int convection_matrix_rhs(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs,
    std::vector<double> htheta, std::vector<double> qtheta, std::vector<double> rtheta,
    double theta, double dx, double dy, int nx, int ny)
{
    double h;
    double q;
    double r;
    double scvf_fac;

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            int p_0  = convection_p_index(i, j, ny); // central point of control volume
            int p_sw = convection_p_index(i - 1, j - 1, ny);  
            int p_s  = convection_p_index(i    , j - 1, ny);  
            int p_se = convection_p_index(i + 1, j - 1, ny);  
            int p_w  = convection_p_index(i - 1, j    , ny);  
            int p_e  = convection_p_index(i + 1, j    , ny);  
            int p_nw = convection_p_index(i - 1, j + 1, ny);  
            int p_n  = convection_p_index(i    , j + 1, ny);  
            int p_ne = convection_p_index(i + 1, j + 1, ny);  

            // int c_eq = 3 * p_0; continuity equation, not applicable for convection term
            int q_eq = 3 * p_0 + 1;
            int r_eq = 3 * p_0 + 2;
            
            // sub control volume 0 ============================================
            // scv_0 face_0
            h = convection_scvf_xi(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
            q = convection_scvf_xi(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
            r = convection_scvf_xi(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
            scvf_fac = theta * -1. * 0.25 * dy * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_w ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_s ) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_w ) += scvf_fac * 3.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw) += scvf_fac * 1.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_s ) += scvf_fac * 1.* convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw + 1) += scvf_fac * 1.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 1) += scvf_fac * 1.* convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw + 2) += scvf_fac * 1.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 2) += scvf_fac * 1.* convection_J_23(h, q, r);

            scvf_fac = -1. * 0.25 * dy;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);

            // scv_0 face_1
            h = convection_scvf_eta(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
            q = convection_scvf_eta(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
            r = convection_scvf_eta(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
            scvf_fac = theta * -1. * 0.25 * dx * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_w ) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_s ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_w ) += scvf_fac * 1.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw) += scvf_fac * 1.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_s ) += scvf_fac * 3.* convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 1) += scvf_fac * 1.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw + 1) += scvf_fac * 1.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_sw + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 2) += scvf_fac * 1.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_sw + 2) += scvf_fac * 1.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);

            scvf_fac = -1. * 0.25 * dy;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);
 
            // sub control volume 1 ============================================
            // scv_1 face_2
            h = convection_scvf_eta(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
            q = convection_scvf_eta(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
            r = convection_scvf_eta(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
            scvf_fac = theta * -1. * 0.25 * dx * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_s ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_se) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_e ) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_s ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_se) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_e ) += scvf_fac * 1. * convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_se + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_se + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 1) += scvf_fac * 1. * convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_se + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_se + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 2) += scvf_fac * 1. * convection_J_23(h, q, r);

            scvf_fac = -1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);

            // scv_1 face_3
            h = convection_scvf_xi(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
            q = convection_scvf_xi(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
            r = convection_scvf_xi(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
            scvf_fac = theta * 1. * 0.25 * dy * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_s ) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_se) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_e ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_s ) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_se) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_e ) += scvf_fac * 3. * convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_se + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_se + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_s  + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_se + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_s  + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_se + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            scvf_fac = 1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);

            // sub control volume 2 ============================================
            // scv_2 face_4
            h = convection_scvf_xi(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
            q = convection_scvf_xi(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
            r = convection_scvf_xi(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
            scvf_fac = theta * 1. * 0.25 * dy * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_e ) += scvf_fac * 3.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_n ) += scvf_fac * 1.* convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_e ) += scvf_fac * 3.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne) += scvf_fac * 1.* convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_n ) += scvf_fac * 1.* convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 1) += scvf_fac * 3.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 1) += scvf_fac * 1.* convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 1) += scvf_fac * 3.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne + 1) += scvf_fac * 1.* convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 1) += scvf_fac * 1.* convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 2) += scvf_fac * 3.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 2) += scvf_fac * 1.* convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 2) += scvf_fac * 3.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne + 2) += scvf_fac * 1.* convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 2) += scvf_fac * 1.* convection_J_23(h, q, r);

            scvf_fac = 1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);


            // scv_2 face_5
            h = convection_scvf_eta(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
            q = convection_scvf_eta(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
            r = convection_scvf_eta(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
            scvf_fac = theta * 1. * 0.25 * dx * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_e ) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_n ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_e ) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_n ) += scvf_fac * 3. * convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_e  + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_ne + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_e  + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_ne + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);

            scvf_fac = 1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);

            // sub control volume 3 ============================================
            // scv_3 face_6
            h = convection_scvf_eta(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
            q = convection_scvf_eta(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
            r = convection_scvf_eta(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
            scvf_fac = theta * 1. * 0.25 * dx * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_n ) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_w ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_n ) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_w ) += scvf_fac * 3. * convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);

            scvf_fac = 1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);

            // scv_3 face_7
            h = convection_scvf_xi(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
            q = convection_scvf_xi(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
            r = convection_scvf_xi(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
            scvf_fac = theta * -1. * 0.25 * dy * 0.125;
            A.coeffRef(q_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_n ) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw) += scvf_fac * 1. * convection_J_11(h, q, r);
            A.coeffRef(q_eq, 3 * p_w ) += scvf_fac * 3. * convection_J_11(h, q, r);
            A.coeffRef(r_eq, 3 * p_0 ) += scvf_fac * 3. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_n ) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw) += scvf_fac * 1. * convection_J_21(h, q, r);
            A.coeffRef(r_eq, 3 * p_w ) += scvf_fac * 3. * convection_J_21(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw + 1) += scvf_fac * 1. * convection_J_12(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 1) += scvf_fac * 3. * convection_J_12(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw + 1) += scvf_fac * 1. * convection_J_22(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 1) += scvf_fac * 3. * convection_J_22(h, q, r);

            A.coeffRef(q_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_n  + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_nw + 2) += scvf_fac * 1. * convection_J_13(h, q, r);
            A.coeffRef(q_eq, 3 * p_w  + 2) += scvf_fac * 3. * convection_J_13(h, q, r);
            A.coeffRef(r_eq, 3 * p_0  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_n  + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_nw + 2) += scvf_fac * 1. * convection_J_23(h, q, r);
            A.coeffRef(r_eq, 3 * p_w  + 2) += scvf_fac * 3. * convection_J_23(h, q, r);

            scvf_fac = -1. * 0.25 * dx;
            rhs[q_eq] += -scvf_fac * convection_J_10(h, q, r);
            rhs[r_eq] += -scvf_fac * convection_J_20(h, q, r);
       }
    }
    return 0;
}
int cconvection_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    int nx, int ny)                          // RHS vector [h, q, r]^{n}
{
    // Convection for post processing; WITHOUT integration over the control volumes.
    double h;
    double q;
    double r;

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int p_0  = convection_p_index(i, j, ny); // central point of control volume
            h = hn[p_0];
            q = qn[p_0];
            r = rn[p_0];
            // q-momentum equation
            rhs_q[p_0] = -( convection_J_10(h, q, r) );

            // r-momentum equation
            rhs_r[p_0] = -( convection_J_20(h, q, r) );
        }
    }
    return 0;
}
inline int convection_p_index(int i, int j, int ny)
{
    return i * ny + j;
}
inline double convection_J_10(double& h, double& q, double& r)
{
    return q * q / h + q * r / h;
}
inline double convection_J_11(double& h, double& q, double& r)
{
    return -q * q / (h * h) - q * r / (h * h);  // d()/dh
}
inline double convection_J_12(double& h, double& q, double& r)
{
    return 2. * q / h + r / h;  // d()/dq
}
inline double convection_J_13(double& h, double& q, double& r)
{
    UNUSED(r);
    return q / h; // d()/dr
}
inline double convection_J_20(double& h, double& q, double& r)
{
    return r * q / h + r * r / h;
}
inline double convection_J_21(double& h, double& q, double& r)
{
    return -r * q / (h * h) - r * r / (h * h);  // d()/dh
}
inline double convection_J_22(double& h, double& q, double& r)
{
    UNUSED(q);
    return r / h; // d()/dq
}
inline double convection_J_23(double& h, double& q, double& r)
{
    return q / h + 2. * r / h; // d()/dr
}

inline double convection_scvf_xi(double c0, double c1, double c2, double c3)
{
    // c_{i+1/2, j+1/4}

//   3 - - - - - - - 2
//   |       |       | 
//   |       |       | 
//   |       |       | 
//   | - - - - - - - | 
//   |       |       | 
//   |       x       | 
//   |       |       | 
//   0 - - - - - - - 1 

    return 0.125 * (3. * c0 + 3. * c1 + 1. * c2 + 1. * c3);
}
inline double convection_scvf_eta(double c0, double c1, double c2, double c3)
{
    // c_{i+1/4, j+1/2}

//   3 - - - - - - - 2
//   |       |       | 
//   |       |       | 
//   |       |       | 
//   | - x - - - - - | 
//   |       |       | 
//   |       |       | 
//   |       |       | 
//   0 - - - - - - - 1 

    return 0.125 * (3. * c0 + 1. * c1 + 1. * c2 + 3. * c3);
}

