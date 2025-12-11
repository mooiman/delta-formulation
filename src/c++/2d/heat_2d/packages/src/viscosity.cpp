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
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& Ttheta,
    std::vector<double>& visc, double theta, size_t nx, size_t ny)
{
    double n_xi = 0.0;  // xi-component of the outward normal vector
    double n_eta = 0.0;  // eta-component of the outward normal vector
    double scvf_fac;

    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    int p_0 = c_eq/(9);  // node number;  // centre of discretization molecule
    // if node number is south or north boundary point, exit the function
    if (std::fmod(p_0, ny) == 0) { return 1; }  // south boundary
    if (std::fmod(p_0 + 1, ny) == 0) { return 2; }  // north boundary

    //std::fill_n(values + c_eq, 9, 0.0);  // set all coefficients for one row of Delta T to zero

    int p_sw = p_0 - ny - 1;
    int p_w  = p_0 - ny;
    int p_nw = p_0 - ny + 1;
    int p_s  = p_0 - 1; 
    int p_n  = p_0 + 1;
    int p_se = p_0 + ny - 1;
    int p_e  = p_0 + ny;
    int p_ne = p_0 + ny + 1;

    std::vector<double> x_pol = scv_nodes(0, x[p_0], x[p_w], x[p_sw], x[p_s]);
    std::vector<double> y_pol = scv_nodes(0, y[p_0], y[p_w], y[p_sw], y[p_s]);
    double scv_area_0 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(1, x[p_0], x[p_s], x[p_se], x[p_e]);
    y_pol = scv_nodes(1, y[p_0], y[p_s], y[p_se], y[p_e]);
    double scv_area_1 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(2, x[p_0], x[p_e], x[p_ne], x[p_n]);
    y_pol = scv_nodes(2, y[p_0], y[p_e], y[p_ne], y[p_n]);
    double scv_area_2 = polygon_area(x_pol, y_pol);

    x_pol = scv_nodes(3, x[p_0], x[p_n], x[p_nw], x[p_w]);
    y_pol = scv_nodes(3, y[p_0], y[p_n], y[p_nw], y[p_w]);
    double scv_area_3 = polygon_area(x_pol, y_pol);

    double cv_area_0 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_1 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_2 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;
    double cv_area_3 = scv_area_0 + scv_area_1 + scv_area_2 + scv_area_3;

    // at face 0
    double dy_dxi_f0  = dcdx_scvf_n(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_dxi_f0  = dcdx_scvf_n(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_deta_f0 = dcdy_scvf_t(y[p_0], y[p_s], y[p_w], y[p_sw]);
    double dx_deta_f0 = dcdy_scvf_t(x[p_0], x[p_s], x[p_w], x[p_sw]);
    // at face 1
    double dy_dxi_f1  = dcdx_scvf_t(y[p_0], y[p_w], y[p_s], y[p_sw]);
    double dx_dxi_f1  = dcdx_scvf_t(x[p_0], x[p_w], x[p_s], x[p_sw]);
    double dy_deta_f1 = dcdy_scvf_n(y[p_0], y[p_s], y[p_w], y[p_sw]);
    double dx_deta_f1 = dcdy_scvf_n(x[p_0], x[p_s], x[p_w], x[p_sw]);
    // at face 2
    double dy_dxi_f2  = dcdx_scvf_t(y[p_e], y[p_0], y[p_se], y[p_s]);
    double dx_dxi_f2  = dcdx_scvf_t(x[p_e], x[p_0], x[p_se], x[p_s]);
    double dy_deta_f2 = dcdy_scvf_n(y[p_0], y[p_s], y[p_e], y[p_se]);
    double dx_deta_f2 = dcdy_scvf_n(x[p_0], x[p_s], x[p_e], x[p_se]);
    // at face 3
    double dy_dxi_f3  = dcdx_scvf_n(y[p_e], y[p_0], y[p_se], y[p_s]);
    double dx_dxi_f3  = dcdx_scvf_n(x[p_e], x[p_0], x[p_se], x[p_s]);
    double dy_deta_f3 = dcdy_scvf_t(y[p_0], y[p_s], y[p_e], y[p_se]);
    double dx_deta_f3 = dcdy_scvf_t(x[p_0], x[p_s], x[p_e], x[p_se]);
    // at face 4
    double dy_dxi_f4  = dcdx_scvf_n(y[p_e], y[p_0], y[p_ne], y[p_n]);
    double dx_dxi_f4  = dcdx_scvf_n(x[p_e], x[p_0], x[p_ne], x[p_n]);
    double dy_deta_f4 = dcdy_scvf_t(y[p_n], y[p_0], y[p_ne], y[p_e]);
    double dx_deta_f4 = dcdy_scvf_t(x[p_n], x[p_0], x[p_ne], x[p_e]);
    // at face 5
    double dy_dxi_f5  = dcdx_scvf_t(y[p_e], y[p_0], y[p_ne], y[p_n]);
    double dx_dxi_f5  = dcdx_scvf_t(x[p_e], x[p_0], x[p_ne], x[p_n]);
    double dy_deta_f5 = dcdy_scvf_n(y[p_n], y[p_0], y[p_ne], y[p_e]);
    double dx_deta_f5 = dcdy_scvf_n(x[p_n], x[p_0], x[p_ne], x[p_e]);
    // at face 6
    double dy_dxi_f6  = dcdx_scvf_t(y[p_0], y[p_w], y[p_n], y[p_nw]);
    double dx_dxi_f6  = dcdx_scvf_t(x[p_0], x[p_w], x[p_n], x[p_nw]);
    double dy_deta_f6 = dcdy_scvf_n(y[p_n], y[p_0], y[p_nw], y[p_w]);
    double dx_deta_f6 = dcdy_scvf_n(x[p_n], x[p_0], x[p_nw], x[p_w]);
    // at face 7
    double dy_dxi_f7  = dcdx_scvf_n(y[p_0], y[p_w], y[p_n], y[p_nw]);
    double dx_dxi_f7  = dcdx_scvf_n(x[p_0], x[p_w], x[p_n], x[p_nw]);
    double dy_deta_f7 = dcdy_scvf_t(y[p_n], y[p_0], y[p_nw], y[p_w]);
    double dx_deta_f7 = dcdy_scvf_t(x[p_n], x[p_0], x[p_nw], x[p_w]);

    //--------------------------------------------------------------------------
    // heat-equation
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

    double dT_dxi_0 = dcdx_scvf_n(Ttheta[p_0], Ttheta[p_w], Ttheta[p_s], Ttheta[p_sw]);
    double dT_deta_0 = dcdy_scvf_t(Ttheta[p_0], Ttheta[p_s], Ttheta[p_w], Ttheta[p_sw]);

    scvf_fac = 0.5 * dy_deta_f0 * n_xi ;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_11 *  dy_deta_f0 * 3./4. ) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_11 *  dy_deta_f0 * 3./4. ) );
    add_value(values, col_s , theta * scvf_fac * ( + psi_11 *  dy_deta_f0 * 1./4. ) );
    add_value(values, col_sw, theta * scvf_fac * ( - psi_11 *  dy_deta_f0 * 1./4. ) );

    rhs[row] += - ( 0.5 *  dy_deta_f0 / cv_area_0 * n_xi * ( psi_11 * ( dy_deta_f0 * dT_dxi_0 ) )
                  );

    // scv_0 face_1
    n_xi =  0.0;
    n_eta = -1.0;

    dT_deta_0 = dcdy_scvf_n(Ttheta[p_0], Ttheta[p_s], Ttheta[p_w], Ttheta[p_sw]);
    dT_dxi_0 = dcdx_scvf_t(Ttheta[p_0], Ttheta[p_w], Ttheta[p_s], Ttheta[p_sw]);

    scvf_fac = 0.5 * dx_dxi_f1 * n_eta ;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_22 *  dx_dxi_f1 * 3./4. ) );
    add_value(values, col_w , theta * scvf_fac * ( + psi_22 *  dx_dxi_f1 * 1./4. ) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_22 *  dx_dxi_f1 * 3./4. ) );
    add_value(values, col_sw, theta * scvf_fac * ( - psi_22 *  dx_dxi_f1 * 1./4. ) );
        
    rhs[row] += - (
                  + 0.5 *  dx_dxi_f1 / cv_area_0 * n_eta * ( psi_22 * (dx_dxi_f1 * dT_deta_0) )
                  );
 
    // sub control volume 1 ============================================
    // scv_1 face_2
    n_xi =  0.0;
    n_eta = -1.0;

    dT_deta_0 = dcdy_scvf_n(Ttheta[p_0], Ttheta[p_s], Ttheta[p_e], Ttheta[p_se]);
    dT_dxi_0 = dcdx_scvf_t(Ttheta[p_e], Ttheta[p_0], Ttheta[p_se], Ttheta[p_s]);

    scvf_fac = 0.5 * dx_dxi_f2 * n_eta ;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_22 * dx_dxi_f2 * 3./4. ) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_22 * dx_dxi_f2 * 3./4. ) );
    add_value(values, col_e , theta * scvf_fac * ( + psi_22 * dx_dxi_f2 * 1./4. ) );
    add_value(values, col_se, theta * scvf_fac * ( - psi_22 * dx_dxi_f2 * 1./4. ) );

    rhs[row] += - (
                  + 0.5 *  dx_dxi_f2 / cv_area_1 * n_eta * ( psi_22 * (dx_dxi_f2 * dT_deta_0) )
                  );
 
    // scv_1 face_3
    n_xi =  1.0;
    n_eta =  0.0;

    dT_dxi_0 = dcdx_scvf_n(Ttheta[p_e], Ttheta[p_0], Ttheta[p_se], Ttheta[p_s]);
    dT_deta_0 = dcdy_scvf_t(Ttheta[p_0], Ttheta[p_s], Ttheta[p_e], Ttheta[p_se]);

    scvf_fac = 0.5 * dy_deta_f3 * n_xi ;
    add_value(values, col_e , theta * scvf_fac * ( + psi_11 * dy_deta_f3 * 3./4. ) );
    add_value(values, col_0 , theta * scvf_fac * ( - psi_11 * dy_deta_f3 * 3./4. ) );
    add_value(values, col_se, theta * scvf_fac * ( + psi_11 * dy_deta_f3 * 1./4. ) );
    add_value(values, col_s , theta * scvf_fac * ( - psi_11 * dy_deta_f3 * 1./4. ) );

    rhs[row] += - ( 0.5 *  dy_deta_f3 / cv_area_1 * n_xi * ( psi_11 * ( dy_deta_f3 * dT_dxi_0 ) )
                  );

    // sub control volume 2 ============================================
    // scv_2 face_4
    n_xi =  1.0;
    n_eta =  0.0;

    dT_dxi_0 = dcdx_scvf_n(Ttheta[p_e], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_n]);
    dT_deta_0 = dcdy_scvf_t(Ttheta[p_n], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_e]);

    scvf_fac = 0.5 * dy_deta_f4 * n_xi ;
    add_value(values, col_e , theta * scvf_fac * ( + psi_11 * 3./4. * dy_deta_f4 ) );
    add_value(values, col_0 , theta * scvf_fac * ( - psi_11 * 3./4. * dy_deta_f4 ) );
    add_value(values, col_ne, theta * scvf_fac * ( + psi_11 * 1./4. * dy_deta_f4 ) );
    add_value(values, col_n , theta * scvf_fac * ( - psi_11 * 1./4. * dy_deta_f4 ) );

    rhs[row] += - (  0.5 *  dy_deta_f4 / cv_area_2 * n_xi * ( psi_11 * ( dy_deta_f4 * dT_dxi_0 ) )
                  );

    // scv_2 face_5
    n_xi =  0.0;
    n_eta =  1.0;

    dT_deta_0 = dcdy_scvf_n(Ttheta[p_n], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_e]);
    dT_dxi_0 = dcdx_scvf_t(Ttheta[p_e], Ttheta[p_0], Ttheta[p_ne], Ttheta[p_n]);

    scvf_fac = 0.5 * dx_dxi_f5 * n_eta ;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_22 * dx_dxi_f5 * 3./4. ) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_22 * dx_dxi_f5 * 3./4. ) );
    add_value(values, col_e , theta * scvf_fac * ( - psi_22 * dx_dxi_f5 * 1./4. ) );
    add_value(values, col_ne, theta * scvf_fac * ( + psi_22 * dx_dxi_f5 * 1./4. ) );

    rhs[row] += - ( 
                  + 0.5 *  dx_dxi_f5 / cv_area_2 * n_eta * ( psi_22 * (dx_dxi_f5 * dT_deta_0) )
                  );

    // sub control volume 3 ============================================
    // scv_3 face_6
    n_xi =  0.0;
    n_eta =  1.0;

    dT_deta_0 = dcdy_scvf_n(Ttheta[p_n], Ttheta[p_0], Ttheta[p_nw], Ttheta[p_w]);
    dT_dxi_0 = dcdx_scvf_t(Ttheta[p_0], Ttheta[p_w], Ttheta[p_n], Ttheta[p_nw]);

    scvf_fac = 0.5 * dx_dxi_f6 * n_eta ;
    add_value(values, col_0 , theta * scvf_fac * ( - psi_22 * dx_dxi_f6 * 3./4. ) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_22 * dx_dxi_f6 * 3./4. ) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_22 * dx_dxi_f6 * 1./4. ) );
    add_value(values, col_nw, theta * scvf_fac * ( + psi_22 * dx_dxi_f6 * 1./4. ) );

    rhs[row] += - (
                  + 0.5 *  dx_dxi_f6 / cv_area_3 * n_eta * ( psi_22 * (dx_dxi_f6 * dT_deta_0) )
                  );
    
    // scv_3 face_7
    n_xi =  -1.0;
    n_eta =  0.0;

    dT_dxi_0 = dcdx_scvf_n(Ttheta[p_0], Ttheta[p_w], Ttheta[p_n], Ttheta[p_nw]);
    dT_deta_0 = dcdy_scvf_t(Ttheta[p_n], Ttheta[p_0], Ttheta[p_nw], Ttheta[p_w]);

    scvf_fac = 0.5 * dy_deta_f7 * n_xi ;
    add_value(values, col_0 , theta * scvf_fac * ( + psi_11 * dy_deta_f7 * 3./4. ) );
    add_value(values, col_w , theta * scvf_fac * ( - psi_11 * dy_deta_f7 * 3./4. ) );
    add_value(values, col_n , theta * scvf_fac * ( + psi_11 * dy_deta_f7 * 1./4. ) );
    add_value(values, col_nw, theta * scvf_fac * ( - psi_11 * dy_deta_f7 * 1./4. ) );

    rhs[row] += - ( 0.5 *  dy_deta_f7 / cv_area_3 * n_xi * ( psi_11 * ( dy_deta_f7 * dT_dxi_0 ) )
                  );

    // source at (0.0, 0.0)
    if (row == nx * ny/2)
    {
        rhs[p_0] += 0.0;
    }

    return 0;
}
//------------------------------------------------------------------------------
int viscosity_post_rhs(std::vector<double>& rhs_q,
    std::vector<double>& Tn, std::vector<double>& visc, 
    size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y)
{
    double nxi;
    double neta;
    double nxi_dl;
    double neta_dl;

    double dx = 1.0;
    double dy = 1.0;
    double dxinv = 1.0/dx;
    double dyinv = 1.0/dy;

    std::fill_n(rhs_q.data(), rhs_q.size(), 0.0);

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


