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

#include "convection.h"

// Just the convection term in xi-direction implemented

//------------------------------------------------------------------------------
CONVECTION::CONVECTION()
{
}
CONVECTION::~CONVECTION()
{
};
CONVECTION::CONVECTION(double theta_in, double dx_in, double dy_in, int nx_in, int ny_in)
{
    theta = theta_in;
    dx = dx_in;
    dy = dy_in;
    nx = nx_in;
    ny = ny_in;
}
int CONVECTION::matrix_2d_q_eq(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs,
    std::vector<double> hp, std::vector<double> qp, std::vector<double> rp)
{
    // First term qq/h
    double J_11 = 0.0;   // qq/h^2 dh
    double J_12 = 0.0;   // 2q/h dq
    double J_13 = 0.0;   // 0
    double c0, c1, c2, c3;

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            int p_0  = p_index(i, j, ny); // central point of control volume
            int p_sw = p_index(i - 1, j - 1, ny);  
            int p_s  = p_index(i    , j - 1, ny);  
            int p_se = p_index(i + 1, j - 1, ny);  
            int p_w  = p_index(i - 1, j    , ny);  
            int p_e  = p_index(i + 1, j    , ny);  
            int p_nw = p_index(i - 1, j + 1, ny);  
            int p_n  = p_index(i    , j + 1, ny);  
            int p_ne = p_index(i + 1, j + 1, ny);  

            int c_eq = 3 * p_0;
            int q_eq = 3 * p_0 + 1;
            int r_eq = 3 * p_0 + 2;
            
            // matrix columns for delta h:  -qq/h^2
            // scv 0
            c0 = -theta * (-qp[p_0]  * qp[p_0]  / (hp[p_0]  * hp[p_0] ));
            c1 = -theta * (-qp[p_w]  * qp[p_w]  / (hp[p_w]  * hp[p_w] ));
            c2 = -theta * (-qp[p_sw] * qp[p_sw] / (hp[p_sw] * hp[p_sw]));
            c3 = -theta * (-qp[p_s]  * qp[p_s]  / (hp[p_s]  * hp[p_s] ));
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_w, 3 * p_sw, 3 * p_s);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_w, 3 * p_sw, 3 * p_s);

            // scv 1
            c0 = -theta * (qp[p_0]  * qp[p_0]  / (hp[p_0]  * hp[p_0] ));
            c1 = -theta * (qp[p_e]  * qp[p_e]  / (hp[p_e]  * hp[p_e] ));
            c2 = -theta * (qp[p_se] * qp[p_se] / (hp[p_se] * hp[p_se]));
            c3 = -theta * (qp[p_s]  * qp[p_s]  / (hp[p_s]  * hp[p_s] ));
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_e, 3 * p_se, 3 * p_s);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_e, 3 * p_se, 3 * p_s);
            
            // scv 2
            c0 = -theta * (qp[p_0]  * qp[p_0]  / (hp[p_0]  * hp[p_0] ));
            c1 = -theta * (qp[p_e]  * qp[p_e]  / (hp[p_e]  * hp[p_e] ));
            c2 = -theta * (qp[p_ne] * qp[p_ne] / (hp[p_ne] * hp[p_ne]));
            c3 = -theta * (qp[p_n]  * qp[p_n]  / (hp[p_n]  * hp[p_n] ));
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_e, 3 * p_ne, 3 * p_n);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_e, 3 * p_ne, 3 * p_n);
          
            // scv 3
            c0 = -theta * (-qp[p_0]  * qp[p_0]  / (hp[p_0]  * hp[p_0] ));
            c1 = -theta * (-qp[p_w]  * qp[p_w]  / (hp[p_w]  * hp[p_w] ));
            c2 = -theta * (-qp[p_nw] * qp[p_nw] / (hp[p_nw] * hp[p_nw]));
            c3 = -theta * (-qp[p_n]  * qp[p_n]  / (hp[p_n]  * hp[p_n] ));
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_w, 3 * p_nw, 3 * p_n);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0, 3 * p_w, 3 * p_nw, 3 * p_n);
            
            // matrix columns for delta q: 2q/h
            // scv 0
            c0 = -theta * (2.0 * qp[p_0]  / hp[p_0]);
            c1 = -theta * (2.0 * qp[p_w] / hp[p_w]);
            c2 = -theta * (2.0 * qp[p_sw] / hp[p_sw]);
            c3 = -theta * (2.0 * qp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);

            // scv 1
            c0 = -theta * (2.0 * qp[p_0]  / hp[p_0]);
            c1 = -theta * (2.0 * qp[p_e] / hp[p_e]);
            c2 = -theta * (2.0 * qp[p_se] / hp[p_se]);
            c3 = -theta * (2.0 * qp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            
            // scv 2
            c0 = -theta * (2.0 * qp[p_0]  / hp[p_0] );
            c1 = -theta * (2.0 * qp[p_e]  / hp[p_e] );
            c2 = -theta * (2.0 * qp[p_ne] / hp[p_ne]);
            c3 = -theta * (2.0 * qp[p_n]  / hp[p_n] );
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
          
            // scv 3
            c0 = -theta * (2.0 * qp[p_0]  / hp[p_0]);
            c1 = -theta * (2.0 * qp[p_w] / hp[p_w]);
            c2 = -theta * (2.0 * qp[p_nw] / hp[p_nw]);
            c3 = -theta * (2.0 * qp[p_n] / hp[p_n]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);
            
             // matrix columns for: r/h dq + q/h dr
            // scv 0
            c0 = -theta * (rp[p_0]  / hp[p_0]);
            c1 = -theta * (rp[p_w] / hp[p_w]);
            c2 = -theta * (rp[p_sw] / hp[p_sw]);
            c3 = -theta * (rp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);

            // scv 1
            c0 = -theta * (rp[p_0]  / hp[p_0]);
            c1 = -theta * (rp[p_e] / hp[p_e]);
            c2 = -theta * (rp[p_se] / hp[p_se]);
            c3 = -theta * (rp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            
            // scv 2
            c0 = -theta * (rp[p_0]  / hp[p_0]);
            c1 = -theta * (rp[p_e] / hp[p_e]);
            c2 = -theta * (rp[p_ne] / hp[p_ne]);
            c3 = -theta * (rp[p_n] / hp[p_n]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
          
            // scv 3
            c0 = -theta * (rp[p_0]  / hp[p_0]);
            c1 = -theta * (rp[p_w] / hp[p_w]);
            c2 = -theta * (rp[p_nw] / hp[p_nw]);
            c3 = -theta * (rp[p_n] / hp[p_n]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);


            // scv 0
            c0 = -theta * (qp[p_0]  / hp[p_0]);
            c1 = -theta * (qp[p_w] / hp[p_w]);
            c2 = -theta * (qp[p_sw] / hp[p_sw]);
            c3 = -theta * (qp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_sw + 1, 3 * p_s + 1);

            // scv 1
            c0 = -theta * (rp[p_0]  / hp[p_0]);
            c1 = -theta * (rp[p_e] / hp[p_e]);
            c2 = -theta * (rp[p_se] / hp[p_se]);
            c3 = -theta * (rp[p_s] / hp[p_s]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_se + 1, 3 * p_s + 1);
            
            // scv 2
            c0 = -theta * (qp[p_0]  / hp[p_0]);
            c1 = -theta * (qp[p_e] / hp[p_e]);
            c2 = -theta * (qp[p_ne] / hp[p_ne]);
            c3 = -theta * (qp[p_n] / hp[p_n]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_e + 1, 3 * p_ne + 1, 3 * p_n + 1);
          
            // scv 3
            c0 = -theta * (qp[p_0]  / hp[p_0]);
            c1 = -theta * (qp[p_w] / hp[p_w]);
            c2 = -theta * (qp[p_nw] / hp[p_nw]);
            c3 = -theta * (qp[p_n] / hp[p_n]);
            c_scvf_xi (A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);
            c_scvf_eta(A, rhs, q_eq, c0, c1, c2, c3, 3 * p_0 + 1, 3 * p_w + 1, 3 * p_nw + 1, 3 * p_n + 1);
            
       }
    }


    return 0;
}
int CONVECTION::rhs_2d_q_eq(Eigen::VectorXd rhs)                          // RHS vector [h, q, r]^{n}
{
    double qqh;   // rhs: qq/h
    double qrh;   // rhs: qr/h


    for (int i = 0; i < rhs.size(); ++i)
    {
        qqh = 0.0;   // rhs: qq/h
        qrh = 0.0;   // rhs: qr/h

        rhs[i] = qqh + qrh;
    }
    return 0;
}
inline int CONVECTION::p_index(int i, int j, int ny)
{
    return i * ny + j;
}
void CONVECTION::c_scvf_xi(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs, int eq, 
    double c0, double c1, double c2, double c3, int p0, int p1, int p2, int p3)
{
    // c_{i+1/2, j+1/4}

    //    3 - - - - 2
    //    |    |    |
    //    | - - - - |
    //    |    x    |
    //    0 - - - - 1

    std::vector<double> a{ 0.125 * 3., 0.125 * 3., 0.125 * 1., 0.125 * 1. };
    
    A.coeffRef(eq, p0) += a[0];
    A.coeffRef(eq, p1) += a[1];
    A.coeffRef(eq, p2) += a[2];
    A.coeffRef(eq, p3) += a[3];
    
    rhs[p0] += a[0] * c0 + a[1] * c1 + a[2] * c2 + a[3] * c3;
    return ;
}
void CONVECTION::c_scvf_eta(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs, int eq, 
    double c0, double c1, double c2, double c3, int p0, int p1, int p2, int p3)
{
    // c_{i+1/2, j+1/4}

    //    3 - - - - 2
    //    |    |    |
    //    |  x  - - |
    //    |    |    |
    //    0 - - - - 1

    std::vector<double> a{ 0.125 * 3., 0.125 * 3., 0.125 * 1., 0.125 * 1. };
    
    A.coeffRef(eq, p0) += a[0];
    A.coeffRef(eq, p1) += a[3];
    A.coeffRef(eq, p2) += a[2];
    A.coeffRef(eq, p3) += a[1];
    
    rhs[p0] += a[0] * c0 + a[1] * c1 + a[2] * c2 + a[3] * c3;
    return ;
}

