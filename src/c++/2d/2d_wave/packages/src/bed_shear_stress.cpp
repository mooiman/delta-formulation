//
// programmer: Jan Mooiman
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

#include "bed_Shear_stress.h"

//------------------------------------------------------------------------------
BEDSHEARSTRESS::BEDSHEARSTRESS()
{
}
BEDSHEARSTRESS::~BEDSHEARSTRESS()
{
};
BEDSHEARSTRESS::BEDSHEARSTRESS(double theta_in, double dx_in, double dy_in, double cf_in, double eps_in, int nx_in, int ny_in)
{
    theta = theta_in;
    dx = dx_in;
    dy = dy_in;
    nx = nx_in;
    ny = ny_in;
    nxny = nx * ny;
    cf = cf_in;
    eps = eps_in;
    cv_area = dx * dy;

    J_11.resize(nxny);  // q_eq: d()/dh
    J_12.resize(nxny);  // q_eq: d()/dq
    J_13.resize(nxny);  // q_eq: d()/dr
    rhs_q.resize(nxny);  // q_eq: right hand side
    std::fill(J_11.begin(), J_11.end(), 0.0);
    std::fill(J_12.begin(), J_12.end(), 0.0);
    std::fill(J_13.begin(), J_13.end(), 0.0);
    std::fill(rhs_q.begin(), rhs_q.end(), 0.0);
    // r-momentum equation
    J_21.resize(nxny);  // r_eq: d()/dh
    J_22.resize(nxny);  // r_eq: d()/dq
    J_23.resize(nxny);  // r_eq: d()/dr
    rhs_r.resize(nxny);  // r_eq: right hand side
    std::fill(J_21.begin(), J_21.end(), 0.0);
    std::fill(J_22.begin(), J_22.end(), 0.0);
    std::fill(J_23.begin(), J_23.end(), 0.0);
    std::fill(rhs_r.begin(), rhs_r.end(), 0.0);

    //std::vector<int> vecOfInts;
    //vecOfInts.resize(10);
    //std::fill(vecOfInts.begin(), vecOfInts.end(), 0);


}
int BEDSHEARSTRESS::matrix_2d_qr_eq(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
    std::vector<double> hp, std::vector<double>& qp, std::vector<double>& rp)
{
    //   nw - - - - - - - n - - - - - - - ne        nw - - - - - - - n - - - - - - - ne
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |       |       |       |          |       |       |       |       |
    //    | - - - - - - - | - - - - - - - |          | - - - - - - - | - - - - - - - |
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |   x   |   x   |       |          |       | scv_3 | scv_2 |       |
    //    |       |       |       |       |          |       |       |       |       |
    //    w - - - - - - - 0 - - - - - - - e          w - - - - - - - 0 - - - - - - - e
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |   x   |   x   |       |          |       | scv_0 | scv_1 |       |
    //    |       |       |       |       |          |       |       |       |       |
    //    | - - - - - - - | - - - - - - - |          | - - - - - - - | - - - - - - - |
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |       |       |       |          |       |       |       |       |
    //    |       |       |       |       |          |       |       |       |       |
    //   sw - - - - - - - s - - - - - - - se        sw - - - - - - - s - - - - - - - se




    // Loop over the nodes first and store the Jacobians
    // q-momentum equation

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int p0  = p_index(i, j, ny); // central point of control volume
            // q-momentum equation
            J_11[p0] = -2. * cf * qp[p0] * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0] * hp[p0]);
            J_12[p0] = cf * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0]) + cf * std::pow(qp[p0], 4.0)/(hp[p0] * hp[p0] * std::pow(vecq(qp[p0], rp[p0]), 3.));
            J_13[p0] = cf * std::pow(qp[p0], 3.0) * rp[p0]/(hp[p0] * hp[p0] * vecq(qp[p0], rp[p0]));

            rhs_q[p0] = -( cf * qp[p0] * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0]) );

            // r-momentum equation
            J_21[p0] = -2. * cf * rp[p0] * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0] * hp[p0]);
            J_22[p0] = cf * std::pow(rp[p0], 3.0) * qp[p0]/(hp[p0] * hp[p0] * vecq(qp[p0], rp[p0]));
            J_23[p0] = cf * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0]) + cf * std::pow(rp[p0], 4.0)/(hp[p0] * hp[p0] * std::pow(vecq(qp[p0], rp[p0]), 3.));

            rhs_r[p0] = -( cf * rp[p0] * vecq(qp[p0], rp[p0]) / (hp[p0] * hp[p0]) );
        }
    }

    // The term s are added to the matrix coeeficients, tehy already contain contributions from other terms in momentum equation
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            int ph_0  = p_index(i, j, ny); // central point of control volume
            int ph_sw = p_index(i - 1, j - 1, ny);  
            int ph_s  = p_index(i    , j - 1, ny);  
            int ph_se = p_index(i + 1, j - 1, ny);  
            int ph_w  = p_index(i - 1, j    , ny);  
            int ph_e  = p_index(i + 1, j    , ny);  
            int ph_nw = p_index(i - 1, j + 1, ny);  
            int ph_n  = p_index(i    , j + 1, ny);  
            int ph_ne = p_index(i + 1, j + 1, ny);  

            //int c_eq = 3 * ph_0;
            int q_eq = 3 * ph_0 + 1;
            int r_eq = 3 * ph_0 + 2;
            //
            // q-momentum equation
            //
            // scv_0
            double scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.*(J_11[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_w ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 1.*(J_11[ph_sw]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.*(J_11[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_s ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 1.*(J_11[ph_se]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.*(J_11[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_e ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 1.*(J_11[ph_ne]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.*(J_11[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_n ]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 1.*(J_11[ph_nw]);
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 3.*(J_11[ph_w ]);

            // scv_0
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_12[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_w ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_12[ph_sw]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_12[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_s ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_12[ph_se]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_12[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_e ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_12[ph_ne]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_12[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_n ]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_12[ph_nw]);
            A.coeffRef(q_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_12[ph_w ]);

            // scv_0
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0  + 2 ) += scv_fac * 9.*(J_13[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_w  + 2 ) += scv_fac * 3.*(J_13[ph_w ]);
            A.coeffRef(q_eq, 3 * ph_sw + 2 ) += scv_fac * 1.*(J_13[ph_sw]);
            A.coeffRef(q_eq, 3 * ph_s  + 2 ) += scv_fac * 3.*(J_13[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_13[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_s  + 2) += scv_fac * 3.*(J_13[ph_s ]);
            A.coeffRef(q_eq, 3 * ph_se + 2) += scv_fac * 1.*(J_13[ph_se]);
            A.coeffRef(q_eq, 3 * ph_e  + 2) += scv_fac * 3.*(J_13[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_13[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_e  + 2) += scv_fac * 3.*(J_13[ph_e ]);
            A.coeffRef(q_eq, 3 * ph_ne + 2) += scv_fac * 1.*(J_13[ph_ne]);
            A.coeffRef(q_eq, 3 * ph_n  + 2) += scv_fac * 3.*(J_13[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_13[ph_0 ]);
            A.coeffRef(q_eq, 3 * ph_n  + 2) += scv_fac * 3.*(J_13[ph_n ]);
            A.coeffRef(q_eq, 3 * ph_nw + 2) += scv_fac * 1.*(J_13[ph_nw]);
            A.coeffRef(q_eq, 3 * ph_w  + 2) += scv_fac * 3.*(J_13[ph_w ]);

            scv_fac = 25 * dx * dy;
            rhs[q_eq] = scv_fac * scv(rhs_q[ph_0], rhs_q[ph_w], rhs_q[ph_sw], rhs_q[ph_s])  // scv_0
                + scv_fac * scv(rhs_q[ph_0], rhs_q[ph_s], rhs_q[ph_se], rhs_q[ph_e])  // scv_1
                + scv_fac * scv(rhs_q[ph_0], rhs_q[ph_e], rhs_q[ph_ne], rhs_q[ph_n])  // scv_2
                + scv_fac * scv(rhs_q[ph_0], rhs_q[ph_n], rhs_q[ph_nw], rhs_q[ph_w]);  //scv_3
            //
            // r-momentum equation
            //
            // scv_0
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.*(J_21[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_w ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 1.*(J_21[ph_sw]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.*(J_21[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_s ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 1.*(J_21[ph_se]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.*(J_21[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_e ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 1.*(J_21[ph_ne]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.*(J_21[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_n ]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 1.*(J_21[ph_nw]);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 3.*(J_21[ph_w ]);

            // scv_0
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_22[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_w ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_22[ph_sw]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_22[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_s ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_22[ph_se]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_22[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_e ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_22[ph_ne]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 9.*(J_22[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_n ]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 1.*(J_22[ph_nw]);
            A.coeffRef(r_eq, 3 * ph_0 + 1) += scv_fac * 3.*(J_22[ph_w ]);

            // scv_0
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0  + 2 ) += scv_fac * 9.*(J_23[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_w  + 2 ) += scv_fac * 3.*(J_23[ph_w ]);
            A.coeffRef(r_eq, 3 * ph_sw + 2 ) += scv_fac * 1.*(J_23[ph_sw]);
            A.coeffRef(r_eq, 3 * ph_s  + 2 ) += scv_fac * 3.*(J_23[ph_s ]);
            //scv_1
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_23[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_s  + 2) += scv_fac * 3.*(J_23[ph_s ]);
            A.coeffRef(r_eq, 3 * ph_se + 2) += scv_fac * 1.*(J_23[ph_se]);
            A.coeffRef(r_eq, 3 * ph_e  + 2) += scv_fac * 3.*(J_23[ph_e ]);
            //scv_2
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_23[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_e  + 2) += scv_fac * 3.*(J_23[ph_e ]);
            A.coeffRef(r_eq, 3 * ph_ne + 2) += scv_fac * 1.*(J_23[ph_ne]);
            A.coeffRef(r_eq, 3 * ph_n  + 2) += scv_fac * 3.*(J_23[ph_n ]);
            //scv_3
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.*(J_23[ph_0 ]);
            A.coeffRef(r_eq, 3 * ph_n  + 2) += scv_fac * 3.*(J_23[ph_n ]);
            A.coeffRef(r_eq, 3 * ph_nw + 2) += scv_fac * 1.*(J_23[ph_nw]);
            A.coeffRef(r_eq, 3 * ph_w  + 2) += scv_fac * 3.*(J_23[ph_w ]);

            scv_fac = 25 * dx * dy;
            rhs[r_eq] = scv_fac * scv(rhs_r[ph_0], rhs_r[ph_w], rhs_r[ph_sw], rhs_r[ph_s])  // scv_0
                + scv_fac * scv(rhs_r[ph_0], rhs_r[ph_s], rhs_r[ph_se], rhs_r[ph_e])  // scv_1
                + scv_fac * scv(rhs_r[ph_0], rhs_r[ph_e], rhs_r[ph_ne], rhs_r[ph_n])  // scv_2
                + scv_fac * scv(rhs_r[ph_0], rhs_r[ph_n], rhs_r[ph_nw], rhs_r[ph_w]);  //scv_3

        }
    }
    return 0;
}
inline int BEDSHEARSTRESS::p_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}
inline double BEDSHEARSTRESS::scv(double& c0, double c1, double c2, double c3)
{
    // value at the quadrature point of a sub control volume
    double value = 0.0625 * (9.0 * c0 + 3.0 * c1 + 1.0 * c2 + 3.0 * c3);
    return value;
}
inline double BEDSHEARSTRESS::vecq(double& qp, double& rp)
{
    double tilde_abs = std::pow(qp*qp*qp*qp + rp*rp*rp*rp + eps*eps*eps*eps,0.25);
    return tilde_abs;
}
