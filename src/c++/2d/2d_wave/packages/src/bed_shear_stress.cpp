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
#include "jacobians.h"

//------------------------------------------------------------------------------
int bed_shear_stress_matrix_rhs(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs,
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
    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    int p_0 = c_eq/(3*27);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (std::fmod(p_0, ny) == 0) { return 1; }  // south boundary
    if (std::fmod(p_0 + 1, ny) == 0) { return 2; }  // north boundary

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

    //--------------------------------------------------------------------------
    // c-equation
    // 
    // No contribution to the continuity equation
    //

    //--------------------------------------------------------------------------
    // q-momentum equation
    // 
    int col_sw = q_eq;
    int col_w  = q_eq + 3;
    int col_nw = q_eq + 6;
    int col_s  = q_eq + 9;
    int col_0  = q_eq + 12;
    int col_n  = q_eq + 15;
    int col_se = q_eq + 18;
    int col_e  = q_eq + 21;
    int col_ne = q_eq + 24;
    //
    // scv_0
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_w , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_sw, scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_s , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));
            
    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_w  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_sw + 1, scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_s  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));
            
    set_value(values, col_0  + 2 , scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_w  + 2 , scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_sw + 2 , scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_s  + 2 , scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 1] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
    //
    //scv_1
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_s , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_se, scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_e , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_s  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_se + 1, scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_e  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_s  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_se + 2, scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_e  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 1] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
    //
    //scv_2
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_e , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_ne, scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_n , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_e  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_ne + 1, scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_n  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_e  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_ne + 2, scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_n  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));

    scv_fac = 0.25 * dx * dy ;
    rhs[row + 1] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
    //
    //scv_3
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_n , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_nw, scv_fac * 1.* bed_shear_stress_J_11(h, q, r, cf));
    set_value(values, col_w , scv_fac * 3.* bed_shear_stress_J_11(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_n  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_nw + 1, scv_fac * 1.* bed_shear_stress_J_12(h, q, r, cf));
    set_value(values, col_w  + 1, scv_fac * 3.* bed_shear_stress_J_12(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_n  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_nw + 2, scv_fac * 1.* bed_shear_stress_J_13(h, q, r, cf));
    set_value(values, col_w  + 2, scv_fac * 3.* bed_shear_stress_J_13(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 1] += -scv_fac * bed_shear_stress_J_10(h, q, r, cf);
    //--------------------------------------------------------------------------
    // r-momentum equation
    // 
    col_sw = r_eq;
    col_w  = r_eq + 3;
    col_nw = r_eq + 6;
    col_s  = r_eq + 9;
    col_0  = r_eq + 12;
    col_n  = r_eq + 15;
    col_se = r_eq + 18;
    col_e  = r_eq + 21;
    col_ne = r_eq + 24;
    //
    // scv_0
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_w], htheta[p_sw], htheta[p_s]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_w], qtheta[p_sw], qtheta[p_s]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_w], rtheta[p_sw], rtheta[p_s]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_w , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_sw, scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_s , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));
            
    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_w  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_sw + 1, scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_s  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));
            
    set_value(values, col_0  + 2 , scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_w  + 2 , scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_sw + 2 , scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_s  + 2 , scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 2] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
    //
    //scv_1
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_s], htheta[p_se], htheta[p_e]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_s], qtheta[p_se], qtheta[p_e]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_s], rtheta[p_se], rtheta[p_e]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_s , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_se, scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_e , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_s  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_se + 1, scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_e  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_s  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_se + 2, scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_e  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 2] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
    //
    //scv_2
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_e], htheta[p_ne], htheta[p_n]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_e], qtheta[p_ne], qtheta[p_n]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_e], rtheta[p_ne], rtheta[p_n]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_e , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_ne, scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_n , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_e  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_ne + 1, scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_n  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_e  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_ne + 2, scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_n  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));

    scv_fac = 0.25 * dx * dy ;
    rhs[row + 2] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);
    //
    //scv_3
    h = bed_shear_stress_scv(htheta[p_0], htheta[p_n], htheta[p_nw], htheta[p_w]);
    q = bed_shear_stress_scv(qtheta[p_0], qtheta[p_n], qtheta[p_nw], qtheta[p_w]);
    r = bed_shear_stress_scv(rtheta[p_0], rtheta[p_n], rtheta[p_nw], rtheta[p_w]);
    scv_fac = theta * 0.25 * dx * dy * 0.0625;
    set_value(values, col_0 , scv_fac * 9.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_n , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_nw, scv_fac * 1.* bed_shear_stress_J_21(h, q, r, cf));
    set_value(values, col_w , scv_fac * 3.* bed_shear_stress_J_21(h, q, r, cf));

    set_value(values, col_0  + 1, scv_fac * 9.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_n  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_nw + 1, scv_fac * 1.* bed_shear_stress_J_22(h, q, r, cf));
    set_value(values, col_w  + 1, scv_fac * 3.* bed_shear_stress_J_22(h, q, r, cf));

    set_value(values, col_0  + 2, scv_fac * 9.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_n  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_nw + 2, scv_fac * 1.* bed_shear_stress_J_23(h, q, r, cf));
    set_value(values, col_w  + 2, scv_fac * 3.* bed_shear_stress_J_23(h, q, r, cf));

    scv_fac = 0.25 * dx * dy;
    rhs[row + 2] += -scv_fac * bed_shear_stress_J_20(h, q, r, cf);

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

    memset(rhs_q.data(), 0, rhs_q.size() * sizeof(double));
    memset(rhs_r.data(), 0, rhs_r.size() * sizeof(double));

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
inline void set_value(double * values, int col, double data){ 
    values[col] += data; 
}

