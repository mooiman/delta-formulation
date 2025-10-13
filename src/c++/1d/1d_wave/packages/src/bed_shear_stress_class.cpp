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

#include "bed_Shear_stress.h"

//------------------------------------------------------------------------------
BEDSHEARSTRESS::BEDSHEARSTRESS()
{
}
BEDSHEARSTRESS::~BEDSHEARSTRESS()
{
};
BEDSHEARSTRESS::BEDSHEARSTRESS(double theta_in, double dx_in, double cf_in, double eps_in, int nx_in)
{
    theta = theta_in;
    dx = dx_in;
    nx = nx_in;
    cf = cf_in;
    eps = eps_in;
    cv_area = dx;
}
int BEDSHEARSTRESS::matrix_rhs(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
    std::vector<double>& hp, std::vector<double>& qp,
    std::vector<double>& hn, std::vector<double>& qn)
{
    //
    //    w - - - | - x - 0 - x - | - - - e          w - - - | scv_0 | scv_1 | - - - e
    //


    std::vector<double> htheta;
    std::vector<double> qtheta;
    htheta.reserve(hn.size());
    qtheta.reserve(hn.size());

    double h;
    double q;
    double scv_fac;

    // Loop over the nodes first and store the Jacobians
    for (int k = 0; k < hn.size(); ++k)
    {
        htheta[k] = theta * hp[k] + (1.0 - theta) * hn[k];
        qtheta[k] = theta * qp[k] + (1.0 - theta) * qn[k];
    }
    
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    for (int i = 1; i < nx-1; ++i)
    {
            int ph_0  = i; // central point of control volume
            int ph_w  = i - 1;  
            int ph_e  = i + 1;  
            
            //int h_eq = 2 * ph_0;
            int q_eq = 2 * ph_0 + 1;

            // scv0
            h = scv(htheta[ph_0], htheta[ph_w]);
            q = scv(qtheta[ph_0], qtheta[ph_w]);
            scv_fac = 0.5 * dx;
            rhs[q_eq] += -scv_fac * J_10(h, q);   // J_10 === cf * q * vecq(q) / (h * h);
            // scv1
            h = scv(htheta[ph_0], htheta[ph_e]);
            q = scv(qtheta[ph_0], qtheta[ph_e]);
            scv_fac = 0.5 * dx;
            rhs[q_eq] += -scv_fac * J_10(h, q);


            // scv0
            h = scv(htheta[ph_0], htheta[ph_w]);
            q = scv(qtheta[ph_0], qtheta[ph_w]);
            scv_fac = theta * 0.5 * dx * 0.25;
            A.coeffRef(q_eq, 2 * ph_0 ) += scv_fac * 3.* J_11(h, q);   // J_11 === -2. * cf * q * vecq(q) / (h * h * h);
            A.coeffRef(q_eq, 2 * ph_w ) += scv_fac * 1.* J_11(h, q);
            // scv1
            h = scv(htheta[ph_0], htheta[ph_e]);
            q = scv(qtheta[ph_0], qtheta[ph_e]);
            scv_fac = theta * 0.5 * dx * 0.25;
            A.coeffRef(q_eq, 2 * ph_0 ) += scv_fac * 3.* J_11(h, q);
            A.coeffRef(q_eq, 2 * ph_e ) += scv_fac * 1.* J_11(h, q);


            // scv0
            h = scv(htheta[ph_0], htheta[ph_w]);
            q = scv(qtheta[ph_0], qtheta[ph_w]);
            scv_fac = theta * 0.5 * dx * 0.25;
            A.coeffRef(q_eq, 2 * ph_0 + 1) += scv_fac * 3.* J_12(h, q);   // J_12 === cf * vecq(q) / (h * h) + cf * std::pow(q, 4.0)/(h * h * std::pow(vecq(q), 3.0));
            A.coeffRef(q_eq, 2 * ph_w + 1) += scv_fac * 1.* J_12(h, q);
            // scv1
            h = scv(htheta[ph_0], htheta[ph_e]);
            q = scv(qtheta[ph_0], qtheta[ph_e]);
            scv_fac = theta * 0.5 * dx * 0.25;
            A.coeffRef(q_eq, 2 * ph_0 + 1) += scv_fac * 3.* J_12(h, q);
            A.coeffRef(q_eq, 2 * ph_e + 1) += scv_fac * 1.* J_12(h, q);
    }
    return 0;
}
void BEDSHEARSTRESS::rhs(std::vector<double>& rhs_q, std::vector<double> hn, std::vector<double>& qn)
{
    // Bed shear stress for post processing; WITHOUT integration over the control volumes.
    double h;
    double q;

    for (int i = 1; i < nx - 1; ++i)
    {
        int p0 = i; // central point of control volume
        h = hn[p0];
        q = qn[p0];
        // q-momentum equation
        rhs_q[p0] = -( J_10(h, q) );
    }
}
inline double BEDSHEARSTRESS::scv(double c0, double c1)
{
    // value at the quadrature point of a sub control volume
    double value = 0.25 * (3.0 * c0 + c1);
    return value;
}
inline double BEDSHEARSTRESS::J_10(double& h, double& q)
{
    // rhs
    return cf * q * vecq(q) / (h * h);
}
inline double BEDSHEARSTRESS::J_11(double& h, double& q)
{
    // d()dh
    return -2. * cf * q * vecq(q) / (h * h * h);
}
inline double BEDSHEARSTRESS::J_12(double& h, double& q)
{
    // d()dq
    return cf * vecq(q) / (h * h) + cf * std::pow(q, 4.0)/(h * h * std::pow(vecq(q), 3.0));
}
inline double BEDSHEARSTRESS::vecq(double& qp)
{
    double tilde_abs = std::pow(qp*qp*qp*qp +  eps*eps*eps*eps, 0.25);
    return tilde_abs;
}
