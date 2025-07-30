//
// programmer: Jan Mooiman
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

#include "bed_shear_stress.h"

//------------------------------------------------------------------------------
int bed_shear_stress_matrix_rhs(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    double cf, double theta, double dx, double dy, int nx, int ny)
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

    double h;
    double q;
    double r;
    double scv_fac;
    
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            int ph_0  = bed_shear_stress_p_index(i, j, ny); // central point of control volume
            int ph_sw = bed_shear_stress_p_index(i - 1, j - 1, ny);  
            int ph_s  = bed_shear_stress_p_index(i    , j - 1, ny);  
            int ph_se = bed_shear_stress_p_index(i + 1, j - 1, ny);  
            int ph_w  = bed_shear_stress_p_index(i - 1, j    , ny);  
            int ph_e  = bed_shear_stress_p_index(i + 1, j    , ny);  
            int ph_nw = bed_shear_stress_p_index(i - 1, j + 1, ny);  
            int ph_n  = bed_shear_stress_p_index(i    , j + 1, ny);  
            int ph_ne = bed_shear_stress_p_index(i + 1, j + 1, ny);  

            //int c_eq = 3 * ph_0;
            int q_eq = 3 * ph_0 + 1;
            int r_eq = 3 * ph_0 + 2;

            //
            // q- and r-momentum equation
            //
            // scv_0
            h = bed_shear_stress_scv(htheta[ph_0], htheta[ph_w], htheta[ph_sw], htheta[ph_s]);
            q = bed_shear_stress_scv(qtheta[ph_0], qtheta[ph_w], qtheta[ph_sw], qtheta[ph_s]);
            r = bed_shear_stress_scv(rtheta[ph_0], rtheta[ph_w], rtheta[ph_sw], rtheta[ph_s]);
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_sw) += scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_sw) += scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);
            
            A.coeffRef(q_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_sw + 1) += scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_sw + 1) += scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);
            
            A.coeffRef(q_eq, 3 * ph_0  + 2 ) += scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w  + 2 ) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_sw + 2 ) += scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s  + 2 ) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 2 ) += scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w  + 2 ) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_sw + 2 ) += scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s  + 2 ) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);

            scv_fac = 0.25 * dx * dy;
            rhs[q_eq] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
            rhs[r_eq] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
            //
            //scv_1
            h = bed_shear_stress_scv(htheta[ph_0], htheta[ph_s], htheta[ph_se], htheta[ph_e]);
            q = bed_shear_stress_scv(qtheta[ph_0], qtheta[ph_s], qtheta[ph_se], qtheta[ph_e]);
            r = bed_shear_stress_scv(rtheta[ph_0], rtheta[ph_s], rtheta[ph_se], rtheta[ph_e]);
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_se) += scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_se) += scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_se + 1) += scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_se + 1) += scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_s  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_se + 2) += scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_s  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_se + 2) += scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);

            scv_fac = 0.25 * dx * dy;
            rhs[q_eq] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
            rhs[r_eq] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
            //
            //scv_2
            h = bed_shear_stress_scv(htheta[ph_0], htheta[ph_e], htheta[ph_ne], htheta[ph_n]);
            q = bed_shear_stress_scv(qtheta[ph_0], qtheta[ph_e], qtheta[ph_ne], qtheta[ph_n]);
            r = bed_shear_stress_scv(rtheta[ph_0], rtheta[ph_e], rtheta[ph_ne], rtheta[ph_n]);
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_ne) += scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_ne) += scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_ne + 1) += scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_ne + 1) += scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_e  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_ne + 2) += scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_e  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_ne + 2) += scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);

            scv_fac = 0.25 * dx * dy ;
            rhs[q_eq] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
            rhs[r_eq] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
            //
            //scv_3
            h = bed_shear_stress_scv(htheta[ph_0], htheta[ph_n], htheta[ph_nw], htheta[ph_w]);
            q = bed_shear_stress_scv(qtheta[ph_0], qtheta[ph_n], qtheta[ph_nw], qtheta[ph_w]);
            r = bed_shear_stress_scv(rtheta[ph_0], rtheta[ph_n], rtheta[ph_nw], rtheta[ph_w]);
            scv_fac = theta * 0.25 * dx * dy * 0.0625;
            A.coeffRef(q_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_nw) += scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w ) += scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0 ) += scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_nw) += scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w ) += scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_nw + 1) += scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w  + 1) += scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 1) += scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_nw + 1) += scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w  + 1) += scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf);

            A.coeffRef(q_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_n  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_nw + 2) += scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(q_eq, 3 * ph_w  + 2) += scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_0  + 2) += scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_n  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_nw + 2) += scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf);
            A.coeffRef(r_eq, 3 * ph_w  + 2) += scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf);

            scv_fac = 0.25 * dx * dy;
            rhs[q_eq] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
            rhs[r_eq] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
        }
    }
    return 0;
}
void bed_shear_stress_rhs(std::vector<double>& rhs_q, std::vector<double>& rhs_r, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    double cf, int nx, int ny)
{
    // Bed shear stress for post processing; WITHOUT integration over the control volumes.
    double h;
    double q;
    double r;

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int p0  = bed_shear_stress_p_index(i, j, ny); // central point of control volume
            h = hn[p0];
            q = qn[p0];
            r = rn[p0];
            // q-momentum equation
            rhs_q[p0] = -( bed_shear_stress_J_10(h, q, r, cf) );

            // r-momentum equation
            rhs_r[p0] = -( bed_shear_stress_J_20(h, q, r, cf) );
        }
    }
}
inline int bed_shear_stress_p_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}
inline double bed_shear_stress_scv(double& c0, double c1, double c2, double c3)
{
    // value at the quadrature point of a sub control volume
    double value = 0.0625 * (9.0 * c0 + 3.0 * c1 + 1.0 * c2 + 3.0 * c3);
    return value;
}
inline double bed_shear_stress_J_10(double& h, double& q, double& r, double&  cf)
{
    return cf * q * vecq(q, r) / (h * h);
}
inline double bed_shear_stress_J_11(double& h, double& q, double& r, double&  cf)
{
    return -2. * cf * q * vecq(q, r) / (h * h * h);
}
inline double bed_shear_stress_J_12(double& h, double& q, double& r, double&  cf)
{
    return cf * vecq(q, r) / (h * h) + cf * std::pow(q, 4.0)/(h * h * std::pow(vecq(q, r), 3.0));
}
inline double bed_shear_stress_J_13(double& h, double& q, double& r, double&  cf)
{
    return cf * q * std::pow(r, 3.0) / (h * h * std::pow(vecq(q, r), 3.0));
}
inline double bed_shear_stress_J_20(double& h, double& q, double& r, double&  cf)
{
    return cf * r * vecq(q, r) / (h * h);
}
inline double bed_shear_stress_J_21(double& h, double& q, double& r, double&  cf)
{
    return -2. * cf * r * vecq(q, r) / (h * h * h);
}
inline double bed_shear_stress_J_22(double& h, double& q, double& r, double&  cf)
{
    return cf * r * std::pow(q, 3.0)/(h * h * std::pow(vecq(q, r), 3.0));
}
inline double bed_shear_stress_J_23(double& h, double& q, double& r, double&  cf)
{
    return cf * vecq(q, r) / (h * h) + cf * std::pow(r, 4.0)/(h * h * std::pow(vecq(q, r), 3.0));
}
inline double vecq(double& qp, double& rp)
{
    double eps = 0.01;
    double tilde_abs = std::pow(qp*qp*qp*qp + rp*rp*rp*rp + eps*eps*eps*eps,0.25);
    return tilde_abs;
}
