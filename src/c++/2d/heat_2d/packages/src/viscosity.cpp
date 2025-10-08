//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formulation and Modified Newton iteration 
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

#include "viscosity.h"
#include "interpolations.h"

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

int viscosity_matrix_and_rhs(double* values, size_t row, int c_eq, Eigen::VectorXd& rhs,
    std::vector<double>& Ttheta,
    std::vector<double>& visc, double theta, double dx, double dy, size_t nx, size_t ny)
{
    double n_xi = 0.0;  // xi-component of the outward normal vector
    double n_eta = 0.0;  // eta-component of the outward normal vector
    double scvf_fac;
    double dxinv = 1./dx;
    double dyinv = 1./dy;

    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    int p_0 = c_eq/(9);  // node number;  // centre of discretization molecule
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
    int col_sw = c_eq;
    int col_w  = c_eq + 1;
    int col_nw = c_eq + 2;
    int col_s  = c_eq + 3;
    int col_0  = c_eq + 4;
    int col_n  = c_eq + 5;
    int col_se = c_eq + 6;
    int col_e  = c_eq + 7;
    int col_ne = c_eq + 8;

    double psi_11 = -visc[p_0];
    double psi_22 = -visc[p_0];
    double psi_12 = 0.0;
    double psi_21 = 0.0;
    //
    // sub control volume 0 ============================================
    // scv_0 face_0
    n_xi = -1.0;
    n_eta =  0.0;

    double Ttheta_0 = scvf_xi(Ttheta[p_0], Ttheta[p_w], Ttheta[p_sw], Ttheta[p_s]);
    double dTtheta_0 = dcdx_scvf_n(Ttheta[p_0], Ttheta[p_w], Ttheta[p_s], Ttheta[p_sw]);

    scvf_fac = 0.5 * dy * n_xi;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_sw, theta * scvf_fac * ( - psi_11 * 0.25 * 1. * dxinv) );
    add_value(values, col_s , theta * scvf_fac * ( + psi_11 * 0.25 * 1. * dxinv) );

    rhs[row] += - ( scvf_fac * psi_11 * dxinv * dTtheta_0 );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    Ttheta_0 = scvf_eta(Ttheta[p_0], Ttheta[p_s], Ttheta[p_sw], Ttheta[p_w]);
    dTtheta_0 = dcdy_scvf_n(Ttheta[p_0], Ttheta[p_s], Ttheta[p_w], Ttheta[p_sw]);

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_22 * 0.25 * 3. * dxinv) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_22 * 0.25 * 3. * dxinv) );
    add_value(values, col_sw, theta * scvf_fac * ( - psi_22 * 0.25 * 1. * dxinv) );
    add_value(values, col_w , theta * scvf_fac * ( + psi_22 * 0.25 * 1. * dxinv) );
        
    rhs[row] += - ( scvf_fac * psi_22 * dyinv * dTtheta_0 );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    Ttheta_0 = scvf_eta(Ttheta[p_0], Ttheta[p_s], Ttheta[p_se], Ttheta[p_e]);
    dTtheta_0 = dcdy_scvf_n(Ttheta[p_0], Ttheta[p_s], Ttheta[p_e], Ttheta[p_se]);

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_22 * 0.25 * 3. * dxinv) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_22 * 0.25 * 3. * dxinv) );
    add_value(values, col_se, theta * scvf_fac * ( - psi_22 * 0.25 * 1. * dxinv) );
    add_value(values, col_e , theta * scvf_fac * ( + psi_22 * 0.25 * 1. * dxinv) );
        
    rhs[row] += - ( scvf_fac * psi_22 * dyinv * dTtheta_0 );

    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    Ttheta_0 = scvf_xi(Ttheta[p_0], Ttheta[p_e], Ttheta[p_se], Ttheta[p_s]);
    dTtheta_0 = dcdx_scvf_n(Ttheta[p_e], Ttheta[p_0], Ttheta[p_se], Ttheta[p_s]);

    scvf_fac = 0.5 * dy * n_xi;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_e , theta * scvf_fac * ( + psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_se, theta * scvf_fac * ( + psi_11 * 0.25 * 1. * dxinv) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_11 * 0.25 * 1. * dxinv) );

    rhs[row] += - scvf_fac * ( psi_11 * dxinv * dTtheta_0 );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    Ttheta_0 = scvf_xi(Ttheta[p_0], Ttheta[p_e], Ttheta[p_ne], Ttheta[p_n]);
    dTtheta_0 = dcdx_scvf_n(Ttheta[p_e], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_n]);

    scvf_fac = 0.5 * dy * n_xi;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_e , theta * scvf_fac * ( + psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_ne, theta * scvf_fac * ( + psi_11 * 0.25 * 1. * dxinv) );
    add_value(values, col_n , theta * scvf_fac * ( - psi_11 * 0.25 * 1. * dxinv) );

    rhs[row] += - ( scvf_fac * psi_11 * dxinv * dTtheta_0 );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    Ttheta_0 = scvf_eta(Ttheta[p_0], Ttheta[p_n], Ttheta[p_ne], Ttheta[p_e]);
    dTtheta_0 = dcdy_scvf_n(Ttheta[p_n], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_e]);

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_22 * 0.25 * 3. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_22 * 0.25 * 3. * dyinv) );
    add_value(values, col_ne, theta * scvf_fac * ( + psi_22 * 0.25 * 1. * dyinv) );
    add_value(values, col_e , theta * scvf_fac * ( - psi_22 * 0.25 * 1. * dyinv) );
    
    rhs[row] += - ( scvf_fac * psi_22 * dyinv * dTtheta_0 );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    Ttheta_0 = scvf_eta(Ttheta[p_0], Ttheta[p_n], Ttheta[p_nw], Ttheta[p_w]);
    dTtheta_0 = dcdy_scvf_n(Ttheta[p_n], Ttheta[p_0], Ttheta[p_nw], Ttheta[p_w]);

    scvf_fac = 0.5 * dx * n_eta;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_22 * 0.25 * 3. * dyinv) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_22 * 0.25 * 3. * dyinv) );
    add_value(values, col_nw, theta * scvf_fac * ( + psi_22 * 0.25 * 1. * dyinv) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_22 * 0.25 * 1. * dyinv) );
    
    rhs[row] += - ( scvf_fac * psi_22 * dyinv * dTtheta_0 );

    // scv_3 face_7
    n_xi =  -1.0;
    n_eta =  0.0;

    Ttheta_0 = scvf_xi(Ttheta[p_0], Ttheta[p_w], Ttheta[p_nw], Ttheta[p_n]);
    dTtheta_0 = dcdx_scvf_n(Ttheta[p_0], Ttheta[p_w], Ttheta[p_n], Ttheta[p_nw]);

    scvf_fac = 0.5 * dy * n_xi;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_11 * 0.25 * 3. * dxinv) );
    add_value(values, col_nw, theta * scvf_fac * ( - psi_11 * 0.25 * 1. * dxinv) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_11 * 0.25 * 1. * dxinv) );

    rhs[row] += - ( scvf_fac * psi_11 * dxinv * dTtheta_0 );

    // source at (0.0, 0.0)
    if (row == nx * ny/2)
    {
        rhs[p_0] += 0.0;
    }

    return 0;
}
int viscosity_post_rhs(std::vector<double>& rhs_q,
    std::vector<double>& Tn,
    std::vector<double>& visc, double dx, double dy, size_t nx, size_t ny)                          // RHS vector [h, q, r]^{n}
{
    // Viscosity for post processing; WITHOUT integration over the control volumes, just the line-integral.
    double nxi;
    double neta;
    double nxi_dl;
    double neta_dl;
    double dxinv = 1./dx;
    double dyinv = 1./dy;

    memset(rhs_q.data(), 0, rhs_q.size() * sizeof(double));

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int p_0  = viscosity_idx(i    , j    , ny); // central point of control volume
            int p_sw = viscosity_idx(i - 1, j - 1, ny);  
            int p_s  = viscosity_idx(i    , j - 1, ny);  
            int p_se = viscosity_idx(i + 1, j - 1, ny);  
            int p_w  = viscosity_idx(i - 1, j    , ny);  
            int p_e  = viscosity_idx(i + 1, j    , ny);  
            int p_nw = viscosity_idx(i - 1, j + 1, ny);  
            int p_n  = viscosity_idx(i    , j + 1, ny);  
            int p_ne = viscosity_idx(i + 1, j + 1, ny);  

            double psi_11 = visc[p_0];
            double psi_22 = visc[p_0];
            double psi_12 = 0.0;
            double psi_21 = 0.0;

            // scv_0 face_0
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = 0.5 * dy;
            neta_dl = 0.0;

            double Tn_0 = scvf_xi(Tn[p_0], Tn[p_w], Tn[p_sw], Tn[p_s]);
            double dTn_0 = dcdx_scvf_n(Tn[p_0], Tn[p_e], Tn[p_s], Tn[p_sw]);

            rhs_q[p_0] += nxi_dl * ( psi_11 * dxinv * dTn_0 );

            // scv_0 face_1
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = neta * 0.5 * dx;

            Tn_0 = scvf_eta(Tn[p_0], Tn[p_s], Tn[p_sw], Tn[p_w]);
            dTn_0 = dcdy_scvf_n(Tn[p_0], Tn[p_s], Tn[p_w], Tn[p_sw]);

            rhs_q[p_0] += nxi_dl * ( psi_11 * dxinv * dTn_0 );

            // sub control volume 1 ============================================
            // scv_1 face_2
            nxi = 0.0;
            neta = -1.0;
            nxi_dl = 0.0;
            neta_dl = neta * 0.5 * dx;

            Tn_0 = scvf_eta(Tn[p_0], Tn[p_s], Tn[p_e], Tn[p_se]);
            dTn_0 = dcdy_scvf_n(Tn[p_0], Tn[p_s], Tn[p_e], Tn[p_se]);

            rhs_q[p_0] += neta_dl * ( psi_22 * dyinv * dTn_0 );

            // scv_1 face_3
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = nxi * 0.5 * dy;
            neta_dl = 0.0;

            Tn_0 = scvf_xi(Tn[p_0], Tn[p_e], Tn[p_se], Tn[p_s]);
            dTn_0 = dcdx_scvf_n(Tn[p_e], Tn[p_0], Tn[p_se], Tn[p_s]);

            rhs_q[p_0] += nxi_dl * ( psi_11 * dxinv * dTn_0 );

            // sub control volume 2 ============================================
            // scv_2 face_4
            nxi = 1.0;
            neta = 0.0;
            nxi_dl = nxi * 0.5 * dy;
            neta_dl = 0.0;

            Tn_0 = scvf_xi(Tn[p_0], Tn[p_e], Tn[p_ne], Tn[p_n]);
            dTn_0 = dcdx_scvf_n(Tn[p_e], Tn[p_0], Tn[p_ne], Tn[p_n]);

            rhs_q[p_0] += nxi_dl * ( psi_11 * dxinv * dTn_0 );

            // scv_2 face_5
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = neta * 0.5 * dx;

            Tn_0 = scvf_eta(Tn[p_0], Tn[p_n], Tn[p_ne], Tn[p_e]);
            dTn_0 = dcdy_scvf_n(Tn[p_n], Tn[p_0], Tn[p_ne], Tn[p_e]);

            rhs_q[p_0] += neta_dl * ( psi_22 * dyinv * dTn_0 );

            // sub control volume 3 ============================================
            // scv_3 face_6
            nxi = 0.0;
            neta = 1.0;
            nxi_dl = 0.0;
            neta_dl = neta * 0.5 * dx;

            Tn_0 = scvf_eta(Tn[p_0], Tn[p_n], Tn[p_nw], Tn[p_w]);
            dTn_0 = dcdy_scvf_n(Tn[p_n], Tn[p_0], Tn[p_nw], Tn[p_w]);

            rhs_q[p_0] += neta_dl * ( psi_22 * dyinv * dTn_0 );

            // scv_3 face_7
            nxi = -1.0;
            neta = 0.0;
            nxi_dl = nxi * 0.5 * dy;
            neta_dl = 0.0;

            Tn_0 = scvf_xi(Tn[p_0], Tn[p_w], Tn[p_nw], Tn[p_n]);
            dTn_0 = dcdx_scvf_n(Tn[p_0], Tn[p_w], Tn[p_n], Tn[p_nw]);

            rhs_q[p_0] += nxi_dl * ( psi_11 * dxinv * dTn_0 );
        }
    }
    return 0;
}
inline size_t viscosity_idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}

inline void add_value(double * values, size_t col, double data){ 
    values[col] += data; 
}


