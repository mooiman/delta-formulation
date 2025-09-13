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

#include "jacobians.h"

//==============================================================================
// bed shear stress
//==============================================================================

double bed_shear_stress_J_10(double& h, double& q, double& r, double&  cf)
{
    return cf * q * abs_vecq(q, r, 1.0) / (h * h);
}
double bed_shear_stress_J_11(double& h, double& q, double& r, double&  cf)
{
    return -2. * cf * q * abs_vecq(q, r, 1.0) / (h * h * h);
}
double bed_shear_stress_J_12(double& h, double& q, double& r, double&  cf)
{
    return cf * abs_vecq(q, r, 1.0) / (h * h) + cf * q * q * (q * q + r * r)/(h * h * abs_vecq(q, r, 3.0));
}
double bed_shear_stress_J_13(double& h, double& q, double& r, double&  cf)
{
    return cf * q * r * (q * q + r * r) / (h * h * abs_vecq(q, r, 3.0));
}
double bed_shear_stress_J_20(double& h, double& q, double& r, double&  cf)
{
    return cf * r * abs_vecq(q, r, 1.0) / (h * h);
}
double bed_shear_stress_J_21(double& h, double& q, double& r, double&  cf)
{
    return -2. * cf * r * abs_vecq(q, r, 1.0) / (h * h * h);
}
double bed_shear_stress_J_22(double& h, double& q, double& r, double&  cf)
{
    return cf * r * q * (q * q + r * r)/(h * h * abs_vecq(q, r, 3.0));
}
double bed_shear_stress_J_23(double& h, double& q, double& r, double&  cf)
{
    return cf * abs_vecq(q, r, 1.0) / (h * h) + cf * r * r * (q * q + r * r)/(h * h * abs_vecq(q, r, 3.0));
}
double abs_vecq(double& q_qp, double& r_qp, double a)
{
    double eps = 0.01;
    double tilde_abs = std::pow((q_qp*q_qp + r_qp*r_qp)*(q_qp*q_qp + r_qp*r_qp) + eps*eps*eps*eps,0.25 * a);
    return tilde_abs;
}

//==============================================================================
// convection
//==============================================================================

double convection_J_10(double& h, double& q, double& r, double nxi, double neta)
{
    return q * q / h * nxi + q * r / h * neta;
}
double convection_J_11(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta h
    return -q * q / (h * h) * nxi - q * r / (h * h) * neta; 
}
double convection_J_12(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta q
    return 2. * q / h * nxi + r / h * neta;
}
double convection_J_13(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta r
    UNUSED(r);
    UNUSED(nxi);
    return q / h * neta;
}
double convection_J_20(double& h, double& q, double& r, double nxi, double neta)
{
    return r * q / h * nxi + r * r / h * neta;
}
double convection_J_21(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta h
    return -r * q / (h * h) * nxi - r * r / (h * h) * neta; 
}
double convection_J_22(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta q
    UNUSED(q);
    UNUSED(neta);
    return r / h * nxi;
}
double convection_J_23(double& h, double& q, double& r, double nxi, double neta)
{
    // Contribution Delta r
    return q / h * nxi + 2. * r / h * neta;
}

