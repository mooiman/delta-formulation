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

#include "interpolations.h"
#include "jacobians.h"
#include "matrix_assembly_boundaries.h"

//------------------------------------------------------------------------------
//   nw - - - - - - - n - - - - - - - ne          2 - - - - - - - 5 - - - - - - - 8
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    w - - - - - - - c - - - - - - - e           1 - - - - - - - 4 - - - - - - - 7
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//    |               |               |           |               |               |
//   sw - - - - - - - s - - - - - - - se          0 - - - - - - - 3 - - - - - - - 6

//
// boundary nodes
//
//==============================================================================
int boundary_north(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_NORTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    std::fill_n(values + c_eq, 3 * 27, 0.0);  // set all coefficients for one row of Delta c-, Delta q- and Delta r-equation to zero

    size_t p_5 = c_eq/(3*27);  // node number of boundary point, ie north point of molecule
    size_t p_4 = p_5 - 1;
    size_t p_3 = p_5 - 2;
    size_t p_0 = p_5 - ny - 2;
    size_t p_1 = p_5 - ny - 1;
    size_t p_2 = p_5 - ny;
    size_t p_6 = p_5 + ny - 2;
    size_t p_7 = p_5 + ny - 1;
    size_t p_8 = p_5 + ny;

    double htheta_b = w_nat[0] * htheta[p_5] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_3];
    double zb_b = w_nat[0] * zb[p_5] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_3];

    double sign = 1.0;
    double h_given = bc[BC_NORTH] - zb_b;
    double h_infty = -zb_b;
    double c_wave = std::sqrt(g * htheta_b);
    double con_fac = c_wave;

    // nnorth
    double dy_deta_w = y[p_2] - y[p_1];
    double dy_deta_b = y[p_5] - y[p_4];
    double dy_deta_e = y[p_8] - y[p_7];
    double dx_dxi_0 = 0.5 * (x[p_5] - x[p_2]) + 0.5 * (x[p_4] - x[p_1]);
    double dx_dxi_1 = 0.5 * (x[p_8] - x[p_5]) + 0.5 * (x[p_7] - x[p_4]);
    double dy_deta_0 = 0.25 * (3. * (y[p_5] - y[p_4]) + 1. * (y[p_2] - y[p_1]));
    double dy_deta_1 = 0.25 * (3. * (y[p_5] - y[p_4]) + 1. * (y[p_8] - y[p_7]));

    if (bc_type[BC_NORTH] == "no_slip" || bc_type[BC_NORTH] == "free_slip")
    {
        if (bc_type[BC_NORTH] == "no_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 5 * 3;
            size_t col_s  = c_eq + 4 * 3;
            size_t col_ss = c_eq + 3 * 3;
            values[col_b ] =  1.0 * theta;
            values[col_s ] = -1.0 * theta;
            values[col_ss] = 0.0;
            rhs[row] = -(htheta[p_5] + zb[p_5] - htheta[p_4] - zb[p_4]);

            // Contribution Delta q
            col_b  = q_eq + 5 * 3;
            col_s  = q_eq + 4 * 3;
            col_ss = q_eq + 3 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_s  + 1] = 0.5 * theta;
            values[col_ss + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 5 * 3;
            col_s  = r_eq + 4 * 3;
            col_ss = r_eq + 3 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_s  + 2]  = 0.5 * theta;
            values[col_ss + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
        else if (bc_type[BC_NORTH] == "free_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 5 * 3;
            size_t col_s  = c_eq + 4 * 3;
            size_t col_ss = c_eq + 3 * 3;
            values[col_b ] =  1.0 * theta;
            values[col_s ] = -1.0 * theta;
            values[col_ss] = 0.0;
            rhs[row] = -(htheta[p_5] + zb[p_5] - htheta[p_4] - zb[p_4]);

            // Contribution Delta q
            col_b  = q_eq + 5 * 3;
            col_s  = q_eq + 4 * 3;
            col_ss = q_eq + 3 * 3;
            values[col_b  + 1] = 1.0 * theta;
            values[col_s  + 1] =-1.0 * theta;
            values[col_ss + 1] = 0.0;
            rhs[row + 1] = -(qtheta[p_5] - qtheta[p_4]);

            // Contribution Delta r
            col_b  = r_eq + 5 * 3;
            col_s  = r_eq + 4 * 3;
            col_ss = r_eq + 3 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_s  + 2]  = 0.5 * theta;
            values[col_ss + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
    }
    else
    {
        if (bc_type[BC_NORTH] == "mooiman")
        {
            // essential boundary condition
            // ----------------------------
            // first row
            size_t col_wb  = c_eq + 2 * 3;
            size_t col_ws  = c_eq + 1 * 3;
            size_t col_wss = c_eq + 0 * 3;
            size_t col_b   = c_eq + 5 * 3;
            size_t col_s   = c_eq + 4 * 3;
            size_t col_ss  = c_eq + 3 * 3;
            size_t col_eb  = c_eq + 8 * 3;
            size_t col_es  = c_eq + 7 * 3;
            size_t col_ess = c_eq + 6 * 3;
            //
            // face_0 
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_wb , face_fac * 1./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_ws , face_fac * 1./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_wss, face_fac * 1./4. * theta * w_ess[2] * -c_wave);
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_s  , face_fac * 3./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_ss , face_fac * 3./4. * theta * w_ess[2] * -c_wave);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_s  , face_fac * 3./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_ss , face_fac * 3./4. * theta * w_ess[2] * -c_wave);
            add_value(values, col_eb , face_fac * 1./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_es , face_fac * 1./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_ess, face_fac * 1./4. * theta * w_ess[2] * -c_wave);

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_s  + 1, 0.0);
            add_value(values, col_ss + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_wb  + 2, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_ws  + 2, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_wss + 2, face_fac * 1./4. * theta * w_ess[2]);
            add_value(values, col_b   + 2, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_s   + 2, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ss  + 2, face_fac * 3./4. * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_s   + 2, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ss  + 2, face_fac * 3./4. * theta * w_ess[2]);
            add_value(values, col_eb  + 2, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_es  + 2, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_ess + 2, face_fac * 1./4. * theta * w_ess[2]);
            //
            // Right hand side
            face_fac = 0.5;
            double htheta_b   = 0.25 * (3.0 * htheta[p_5] + 1.0 * htheta[p_2]);
            double htheta_jm1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_1]);
            double htheta_jm2 = 0.25 * (3.0 * htheta[p_3] + 1.0 * htheta[p_0]);
            double ht_jm12 = face_fac * ( w_ess[0] * htheta_b + w_ess[1] * htheta_jm1 + w_ess[2] * htheta_jm2);
            htheta_b   = 0.25 * (3.0 * htheta[p_5] + 1.0 * htheta[p_8]);
            htheta_jm1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_7]);
            htheta_jm2 = 0.25 * (3.0 * htheta[p_3] + 1.0 * htheta[p_6]);
            ht_jm12 += face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_jm1 + w_ess[2] * htheta_jm2);
            
            double rtheta_b   = 0.25 * (3.0 * rtheta[p_5] + 1.0 * rtheta[p_2]);
            double rtheta_jm1 = 0.25 * (3.0 * rtheta[p_4] + 1.0 * rtheta[p_1]);
            double rtheta_jm2 = 0.25 * (3.0 * rtheta[p_3] + 1.0 * rtheta[p_0]);
            double rt_jm12 = face_fac * (w_ess[0] * rtheta_b + w_ess[1] * rtheta_jm1 + w_ess[2] * rtheta_jm2);
            rtheta_b   = 0.25 * (3.0 * rtheta[p_5] + 1.0 * rtheta[p_8]);
            rtheta_jm1 = 0.25 * (3.0 * rtheta[p_4] + 1.0 * rtheta[p_7]);
            rtheta_jm2 = 0.25 * (3.0 * rtheta[p_3] + 1.0 * rtheta[p_6]);
            rt_jm12 += face_fac * (w_ess[0] * rtheta_b + w_ess[1] * rtheta_jm1 + w_ess[2] * rtheta_jm2);

            rhs[row] = -0.5 * (dx_dxi_0 + dx_dxi_1) * ( rt_jm12 - c_wave * (ht_jm12 - h_given) );
            
            if (bc_vars[BC_NORTH] == "zeta")
            {
                rhs[row] += - 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * c_wave * bc[BC_NORTH];
            }
            if (bc_vars[BC_NORTH] == "q")
            {
                rhs[row] += - 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * bc[BC_NORTH];
            }
        }
        else if (bc_type[BC_NORTH] == "borsboom")
        {
            //----------------------------------------------------------------------
            // Essential boundary condition
            //----------------------------------------------------------------------
            double hp_w = w_ess[0] * hp[p_2] + w_ess[1] * hp[p_1] + w_ess[2] * hp[p_0];
            double hp_b = w_ess[0] * hp[p_5] + w_ess[1] * hp[p_4] + w_ess[2] * hp[p_3];
            double hp_e = w_ess[0] * hp[p_8] + w_ess[1] * hp[p_7] + w_ess[2] * hp[p_6];
        
            double hn_w = w_ess[0] * hn[p_2] + w_ess[1] * hn[p_1] + w_ess[2] * hn[p_0];
            double hn_b = w_ess[0] * hn[p_5] + w_ess[1] * hn[p_4] + w_ess[2] * hn[p_3];
            double hn_e = w_ess[0] * hn[p_8] + w_ess[1] * hn[p_7] + w_ess[2] * hn[p_6];

            double htheta_w = w_ess[0] * htheta[p_2] + w_ess[1] * htheta[p_1] + w_ess[2] * htheta[p_0];  // west of boundary location
            double htheta_b = w_ess[0] * htheta[p_5] + w_ess[1] * htheta[p_4] + w_ess[2] * htheta[p_3];
            double htheta_e = w_ess[0] * htheta[p_8] + w_ess[1] * htheta[p_7] + w_ess[2] * htheta[p_6];  // east of boundary location

            double rp_w = w_ess[0] * rp[p_2] + w_ess[1] * rp[p_1] + w_ess[2] * rp[p_0];
            double rp_b = w_ess[0] * rp[p_5] + w_ess[1] * rp[p_4] + w_ess[2] * rp[p_3];
            double rp_e = w_ess[0] * rp[p_8] + w_ess[1] * rp[p_7] + w_ess[2] * rp[p_6];
        
            double rn_w = w_ess[0] * rn[p_2] + w_ess[1] * rn[p_1] + w_ess[2] * rn[p_0];
            double rn_b = w_ess[0] * rn[p_5] + w_ess[1] * rn[p_4] + w_ess[2] * rn[p_3];
            double rn_e = w_ess[0] * rn[p_8] + w_ess[1] * rn[p_7] + w_ess[2] * rn[p_6];

            double rtheta_w = w_ess[0] * rtheta[p_2] + w_ess[1] * rtheta[p_1] + w_ess[2] * rtheta[p_0];  // west of boundary location
            double rtheta_b = w_ess[0] * rtheta[p_5] + w_ess[1] * rtheta[p_4] + w_ess[2] * rtheta[p_3];
            double rtheta_e = w_ess[0] * rtheta[p_8] + w_ess[1] * rtheta[p_7] + w_ess[2] * rtheta[p_6];  // east of boundary location
            //------------------------------------------------------------------
            // first row
            //
            size_t col_wb  = c_eq + 2 * 3;
            size_t col_ws  = c_eq + 1 * 3;
            size_t col_wss = c_eq + 0 * 3;
            size_t col_b   = c_eq + 5 * 3;
            size_t col_s   = c_eq + 4 * 3;
            size_t col_ss  = c_eq + 3 * 3;
            size_t col_eb  = c_eq + 8 * 3;
            size_t col_es  = c_eq + 7 * 3;
            size_t col_ess = c_eq + 6 * 3;
            //
            if (do_convection) { con_fac = c_wave + rp_b / hp_b; }
            //
            // Contribution Delta h
            // face 0 
            double face_fac = 0.5 * dy_deta_0;
            add_value(values, col_wb , face_fac * 1./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_ws , face_fac * 1./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_wss, face_fac * 1./4. * dtinv * -con_fac * w_ess[2]);
            add_value(values, col_b  , face_fac * 3./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_s  , face_fac * 3./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_ss , face_fac * 3./4. * dtinv * -con_fac * w_ess[2]);

            // face_1 
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b  , face_fac * 3./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_s  , face_fac * 3./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_ss , face_fac * 3./4. * dtinv * -con_fac * w_ess[2]);
            add_value(values, col_eb , face_fac * 1./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_es , face_fac * 1./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_ess, face_fac * 1./4. * dtinv * -con_fac * w_ess[2]);

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_s  + 1, 0.0);
            add_value(values, col_ss + 1, 0.0);
            //
            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dy_deta_0;
            add_value(values, col_wb  + 2, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_ws  + 2, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_wss + 2, face_fac * 1./4. * dtinv * w_ess[2]);
            add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_s   + 2, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ss  + 2, face_fac * 3./4. * dtinv * w_ess[2]);

            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_s   + 2, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ss  + 2, face_fac * 3./4. * dtinv * w_ess[2]);
            add_value(values, col_eb  + 2, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_es  + 2, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_ess + 2, face_fac * 1./4. * dtinv * w_ess[2]);
            //
            double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_w);
            double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_e);
            double dhdt_w = dtinv * (hp_w - hn_w);
            double dhdt_b = dtinv * (hp_b - hn_b);
            double dhdt_e = dtinv * (hp_e - hn_e);
            double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_w);
            double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_e);
            double drdt_w = dtinv * (rp_w - rn_w);
            double drdt_b = dtinv * (rp_b - rn_b);
            double drdt_e = dtinv * (rp_e - rn_e);
            double drdt_0 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_w);
            double drdt_1 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_e);
            rhs[row] = - ( 0.5 * dy_deta_0 * (drdt_0 - con_fac * dhdt_0)
                         + 0.5 * dy_deta_1 * (drdt_1 - con_fac * dhdt_1) 
                         );

            double corr_term = 0.0;
            if (bc_vars[BC_NORTH] == "zeta")
            {
                // face 0
                face_fac = 0.5 * dy_deta_0;
                add_value(values, col_wb , face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_ws , face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_wss, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b  , face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_s  , face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ss , face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dy_deta_1;
                add_value(values, col_b  , face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_s  , face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ss , face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_eb , face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_es , face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ess, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                double zb_w = w_ess[0] * zb[p_2] + w_ess[1] * zb[p_1] + w_ess[2] * zb[p_0];
                double zb_b = w_ess[0] * zb[p_5] + w_ess[1] * zb[p_4] + w_ess[2] * zb[p_3];
                double zb_e = w_ess[0] * zb[p_8] + w_ess[1] * zb[p_7] + w_ess[2] * zb[p_6];
                double zb_0 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_w);
                double zb_1 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_e);

                corr_term = - ( 0.5 * dy_deta_0 * ( dhdt_0 + ( eps_bc_corr * ((bc[BC_NORTH] - zb_0) - htheta_0)) ) 
                              + 0.5 * dy_deta_1 * ( dhdt_1 + ( eps_bc_corr * ((bc[BC_NORTH] - zb_1) - htheta_1)) ) 
                              );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_NORTH] == "q")
            {
                // face 0
                face_fac = 0.5 * dy_deta_0;
                add_value(values, col_wb  + 2, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_ws  + 2, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_wss + 2, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b   + 2, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_s   + 2, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ss  + 2, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                //face 1
                face_fac = 0.5 * dy_deta_1;
                add_value(values, col_b   + 2, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_s   + 2, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ss  + 2, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_eb  + 2, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_es  + 2, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ess + 2, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                double rtheta_0 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_w);
                double rtheta_1 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_e);

                corr_term = ( 0.5 * dy_deta_0 * (- drdt_0 + sign * eps_bc_corr * (bc[BC_NORTH]- rtheta_0))  
                            + 0.5 * dy_deta_1 * (- drdt_1 + sign * eps_bc_corr * (bc[BC_NORTH]- rtheta_1))
                            );
                rhs[row] += corr_term;
            }
        }
        //----------------------------------------------------------------------
        // natural boundary condition
        // ---------------------------------------------------------------------
        double hp_w = w_nat[0] * hp[p_2] + w_nat[1] * hp[p_1] + w_nat[2] * hp[p_0];
        double hp_b = w_nat[0] * hp[p_5] + w_nat[1] * hp[p_4] + w_nat[2] * hp[p_3];
        double hp_e = w_nat[0] * hp[p_8] + w_nat[1] * hp[p_7] + w_nat[2] * hp[p_6];
        
        double hn_w = w_nat[0] * hn[p_2] + w_nat[1] * hn[p_1] + w_nat[2] * hn[p_0];
        double hn_b = w_nat[0] * hn[p_5] + w_nat[1] * hn[p_4] + w_nat[2] * hn[p_3];
        double hn_e = w_nat[0] * hn[p_8] + w_nat[1] * hn[p_7] + w_nat[2] * hn[p_6];

        double drdy_w = rtheta[p_2] - rtheta[p_1];
        double drdy_b = rtheta[p_5] - rtheta[p_4];
        double drdy_e = rtheta[p_8] - rtheta[p_7];

        double rp_w = w_nat[0] * rp[p_2] + w_nat[1] * rp[p_1] + w_nat[2] * rp[p_0];
        double rp_b = w_nat[0] * rp[p_5] + w_nat[1] * rp[p_4] + w_nat[2] * rp[p_3];
        double rp_e = w_nat[0] * rp[p_8] + w_nat[1] * rp[p_7] + w_nat[2] * rp[p_6];
        
        double rn_w = w_nat[0] * rn[p_2] + w_nat[1] * rn[p_1] + w_nat[2] * rn[p_0];
        double rn_b = w_nat[0] * rn[p_5] + w_nat[1] * rn[p_4] + w_nat[2] * rn[p_3];
        double rn_e = w_nat[0] * rn[p_8] + w_nat[1] * rn[p_7] + w_nat[2] * rn[p_6];
        
        double htheta_w = w_nat[0] * htheta[p_2] + w_nat[1] * htheta[p_1] + w_nat[2] * htheta[p_0];  // west of boundary location
        double htheta_b = w_nat[0] * htheta[p_5] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_3];
        double htheta_e = w_nat[0] * htheta[p_8] + w_nat[1] * htheta[p_7] + w_nat[2] * htheta[p_6];  // east of boundary location

        double qtheta_w = w_nat[0] * qtheta[p_2] + w_nat[1] * qtheta[p_1] + w_nat[2] * qtheta[p_0];  // west of boundary location
        double qtheta_b = w_nat[0] * qtheta[p_5] + w_nat[1] * qtheta[p_4] + w_nat[2] * qtheta[p_3];
        double qtheta_e = w_nat[0] * qtheta[p_8] + w_nat[1] * qtheta[p_7] + w_nat[2] * qtheta[p_6];  // east of boundary location

        double rtheta_w = w_nat[0] * rtheta[p_2] + w_nat[1] * rtheta[p_1] + w_nat[2] * rtheta[p_0];  // west of boundary location
        double rtheta_b = w_nat[0] * rtheta[p_5] + w_nat[1] * rtheta[p_4] + w_nat[2] * rtheta[p_3];
        double rtheta_e = w_nat[0] * rtheta[p_8] + w_nat[1] * rtheta[p_7] + w_nat[2] * rtheta[p_6];  // east of boundary location

        double dzetady_w = htheta[p_2] + zb[p_2] - htheta[p_1] - zb[p_1];
        double dzetady_b = htheta[p_5] + zb[p_5] - htheta[p_4] - zb[p_4];
        double dzetady_e = htheta[p_8] + zb[p_8] - htheta[p_7] - zb[p_7];
        // ---------------------------------------------------------------------
        if (do_convection) { con_fac = c_wave - rp_b / hp_b; }
        // ---------------------------------------------------------------------
        // second row
        // q-momentum (tangential equation, q == 0)
        //
        size_t col_wb  = q_eq + 2 * 3;
        size_t col_ws  = q_eq + 1 * 3;
        size_t col_wss = q_eq + 0 * 3;
        size_t col_b   = q_eq + 5 * 3;
        size_t col_s   = q_eq + 4 * 3;
        size_t col_ss  = q_eq + 3 * 3;
        size_t col_eb  = q_eq + 8 * 3;
        size_t col_es  = q_eq + 7 * 3;
        size_t col_ess = q_eq + 6 * 3;
        //
        add_value(values, col_b , 0.0);
        add_value(values, col_s , 0.0);
        add_value(values, col_ss, 0.0);
        //
        add_value(values, col_b  + 1, 1.0);
        add_value(values, col_s  + 1, 0.0);
        add_value(values, col_ss + 1, 0.0);
        //
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_s  + 2, 0.0);
        add_value(values, col_ss + 2, 0.0);
        //
        rhs[row + 1] = 0.0;
        //
        // -----------------------------------------------------------------
        // third row
        // momentum part dr/dt + gh d(zeta)/dy
        //
        col_wb  = r_eq + 2 * 3;
        col_ws  = r_eq + 1 * 3;
        col_wss = r_eq + 0 * 3;
        col_b   = r_eq + 5 * 3;
        col_s   = r_eq + 4 * 3;
        col_ss  = r_eq + 3 * 3;
        col_eb  = r_eq + 8 * 3;
        col_es  = r_eq + 7 * 3;
        col_ess = r_eq + 6 * 3;
        //
        // Contribution Delta h
        // face 0
        double face_fac = 0.5;
        add_value(values, col_wb , face_fac * 1./4. * w_nat[0] * theta * g * dzetady_e +  face_fac * 1./4. * theta * g * htheta_e);
        add_value(values, col_ws , face_fac * 1./4. * w_nat[1] * theta * g * dzetady_e -  face_fac * 1./4. * theta * g * htheta_e);
        add_value(values, col_wss, face_fac * 1./4. * w_nat[2] * theta * g * dzetady_e);
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetady_b +  face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_s  , face_fac * 3./4. * w_nat[1] * theta * g * dzetady_b -  face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ss , face_fac * 3./4. * w_nat[2] * theta * g * dzetady_b);
        
        //face 1
        face_fac = 0.5;
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetady_b +  face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_s  , face_fac * 3./4. * w_nat[1] * theta * g * dzetady_b -  face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ss , face_fac * 3./4. * w_nat[2] * theta * g * dzetady_b);
        add_value(values, col_eb , face_fac * 1./4. * w_nat[0] * theta * g * dzetady_w +  face_fac * 1./4. * theta * g * htheta_w);
        add_value(values, col_es , face_fac * 1./4. * w_nat[1] * theta * g * dzetady_w -  face_fac * 1./4. * theta * g * htheta_w);
        add_value(values, col_ess, face_fac * 1./4. * w_nat[2] * theta * g * dzetady_w);

        // Contribution Delta q
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_s  + 1, 0.0);
        add_value(values, col_ss + 1, 0.0);

        // Contribution Delta r
        // face 0
        face_fac = 0.5 * dx_dxi_0;
        add_value(values, col_wb  + 2, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_ws  + 2, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_wss + 2, face_fac * 1./4. * dtinv * w_nat[2]);
        add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_s   + 2, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ss  + 2, face_fac * 3./4. * dtinv * w_nat[2]);
        
        //face 1
        face_fac = 0.5 * dx_dxi_1;
        add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_s   + 2, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ss  + 2, face_fac * 3./4. * dtinv * w_nat[2]);
        add_value(values, col_eb  + 2, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_es  + 2, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_ess + 2, face_fac * 1./4. * dtinv * w_nat[2]);
        //
        double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_w);
        double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_e);
        double rtheta_0 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_w);
        double rtheta_1 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_e);
        double dzetady_0 = 0.25 * (3.0 * dzetady_b + 1.0 * dzetady_w);
        double dzetady_1 = 0.25 * (3.0 * dzetady_b + 1.0 * dzetady_e);
        double drdt_w = dtinv * (rp_w - rn_w);
        double drdt_b = dtinv * (rp_b - rn_b);
        double drdt_e = dtinv * (rp_e - rn_e);
        double drdt_0 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_w);
        double drdt_1 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_e);

        rhs[row + 2] = - ( 0.5 * ( dx_dxi_0 * drdt_0 + g * htheta_0 * dzetady_0 ) 
                         + 0.5 * ( dx_dxi_1 * drdt_1 + g * htheta_1 * dzetady_1 ) );
        //
        if (do_convection)
        {
            // North boundary convection
            
            double dh_deta_0 = 0.25 * (3.0 * (htheta[p_5] - htheta[p_4]) + 1.0 * (htheta[p_2] - htheta[p_1]));
            double dr_deta_0 = 0.25 * (3.0 * (rtheta[p_5] - rtheta[p_4]) + 1.0 * (rtheta[p_2] - rtheta[p_1]));
            double dh_deta_1 = 0.25 * (3.0 * (htheta[p_5] - htheta[p_4]) + 1.0 * (htheta[p_8] - htheta[p_7]));
            double dr_deta_1 = 0.25 * (3.0 * (rtheta[p_5] - rtheta[p_4]) + 1.0 * (rtheta[p_8] - rtheta[p_7]));
            
            double dh_dxi_0 = 0.5 * ((htheta[p_4] - htheta[p_1]) + (htheta[p_5] - htheta[p_2]));
            double dr_dxi_0 = 0.5 * ((rtheta[p_4] - rtheta[p_1]) + (rtheta[p_5] - rtheta[p_2]));
            double dh_dxi_1 = 0.5 * ((htheta[p_4] - htheta[p_1]) + (htheta[p_6] - htheta[p_3]));
            double dr_dxi_1 = 0.5 * ((rtheta[p_4] - rtheta[p_1]) + (rtheta[p_6] - rtheta[p_3]));

            double dx_deta_0 = dcdy_scvf_n(y[p_5], y[p_4], y[p_2], y[p_1]);
            double dx_deta_1 = dcdy_scvf_n(y[p_5], y[p_4], y[p_8], y[p_7]);

            double aa_0 = - 2. * rtheta_0 / (htheta_0 * htheta_0) * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) + 2. * (rtheta_0 * rtheta_0) / (htheta_0 * htheta_0 * htheta_0) * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0);
            double aa_1 = - 2. * rtheta_1 / (htheta_1 * htheta_1) * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) + 2. * (rtheta_1 * rtheta_1) / (htheta_1 * htheta_1 * htheta_1) * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1);
            
            double bb_0 = 2. / htheta_0 * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) - 2. * rtheta_0 / (htheta_0 * htheta_0) * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0);
            double bb_1 = 2. / htheta_1 * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) - 2. * rtheta_1 / (htheta_1 * htheta_1) * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1);
            
            double cc_0 = - rtheta_0 * rtheta_0 / (htheta_0 * htheta_0);
            double cc_1 = - rtheta_1 * rtheta_1 / (htheta_1 * htheta_1);

            double dd_0 = 2. * rtheta_0 / htheta_0;
            double dd_1 = 2. * rtheta_1 / htheta_1;
            
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * theta / dy_deta_0;
            add_value(values, col_wb , face_fac * 1./4. * (aa_0 * w_nat[0] + dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ws , face_fac * 1./4. * (aa_0 * w_nat[1] - dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_wss, face_fac * 1./4. * (aa_0 * w_nat[2]));
            add_value(values, col_b  , face_fac * 3./4. * (aa_0 * w_nat[0] + dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_s  , face_fac * 3./4. * (aa_0 * w_nat[1] - dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ss , face_fac * 3./4. * (aa_0 * w_nat[2]));

            //face 1
            face_fac = 0.5 * theta / dy_deta_1;
            add_value(values, col_b  , face_fac * 3./4. * (aa_1 * w_nat[0] + dx_dxi_1 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_s  , face_fac * 3./4. * (aa_1 * w_nat[1] - dx_dxi_1 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ss , face_fac * 3./4. * (aa_1 * w_nat[2]));
            add_value(values, col_eb , face_fac * 1./4. * (aa_1 * w_nat[0] + dx_dxi_1 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_es , face_fac * 1./4. * (aa_1 * w_nat[1] - dx_dxi_1 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ess, face_fac * 1./4. * (aa_1 * w_nat[2]));
            
            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_s  + 1, 0.0);
            add_value(values, col_ss + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * theta / dy_deta_0;
            add_value(values, col_wb  + 2, face_fac * 1./4. * (bb_0 * w_nat[0] + dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ws  + 2, face_fac * 1./4. * (bb_0 * w_nat[1] - dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_wss + 2, face_fac * 1./4. * (bb_0 * w_nat[2]));
            add_value(values, col_b   + 2, face_fac * 3./4. * (bb_0 * w_nat[0] + dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_s   + 2, face_fac * 3./4. * (bb_0 * w_nat[1] - dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_ss  + 2, face_fac * 3./4. * (bb_0 * w_nat[2]));

            // face 1
            face_fac = 0.5 * theta / dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * (bb_1 * w_nat[0] + dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_s   + 2, face_fac * 3./4. * (bb_1 * w_nat[1] - dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_ss  + 2, face_fac * 3./4. * (bb_1 * w_nat[2]));
            add_value(values, col_eb  + 2, face_fac * 1./4. * (bb_1 * w_nat[0] + dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_es  + 2, face_fac * 1./4. * (bb_1 * w_nat[1] - dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_ess + 2, face_fac * 1./4. * (bb_1 * w_nat[2]));

            rhs[row + 2] += - (
                0.5 * dd_0 * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) / dy_deta_0
              + 0.5 * cc_0 * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0) / dy_deta_0
              + 0.5 * dd_1 * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) / dy_deta_1
              + 0.5 * cc_1 * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1) / dy_deta_1
                );
        }
        if (do_bed_shear_stress)
        {
            // North boundary bed shear stress

            double cf_w = cf;
            double cf_b = cf;
            double cf_e = cf;
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dy_deta_0;
            add_value(values, col_wb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_ws , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wss, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_s  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ss , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_s  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ss , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_eb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_es , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_ess, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_s  + 1, 0.0);
            add_value(values, col_ss + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dy_deta_0;
            add_value(values, col_wb  + 2, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_ws  + 2, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wss + 2, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );
            add_value(values, col_b   + 2, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_s   + 2, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ss  + 2, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_s   + 2, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ss  + 2, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_eb  + 2, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_es  + 2, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_ess + 2, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    

            // right hand side
            double qtheta_0 = 0.0;  // no tangential discharge q
            double qtheta_1 = 0.0;
            double abs_qtheta_0 = abs_vecq(qtheta_0, rtheta_0, 1.0);
            double abs_qtheta_1 = abs_vecq(qtheta_1, rtheta_1, 1.0);

            double cf_0 = 0.25 * (3.0 * cf_b + 1.0 * cf_w);
            double cf_1 = 0.25 * (3.0 * cf_b + 1.0 * cf_e);

            rhs[row + 2] += -(
                  0.5 * dy_deta_0 * cf_0 * rtheta_0 * abs_qtheta_0 / (htheta_0 * htheta_0)
                + 0.5 * dy_deta_1 * cf_1 * rtheta_1 * abs_qtheta_1 / (htheta_1 * htheta_1)
            );
        }
        do_viscosity = false;  // north boundary
        if (do_viscosity)
        {
            // North boundary viscosity
        }
        //
        // continuity part +c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        // face 0
        face_fac = 0.5 * dy_deta_0;
        add_value(values, col_wb , face_fac * 1./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_ws , face_fac * 1./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_wss, face_fac * 1./4. * con_fac * dtinv * w_nat[2]);
        add_value(values, col_b  , face_fac * 3./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_s  , face_fac * 3./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_ss , face_fac * 3./4. * con_fac * dtinv * w_nat[2]);

        // face 1
        face_fac = 0.5 * dy_deta_1;
        add_value(values, col_b  , face_fac * 3./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_s  , face_fac * 3./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_ss , face_fac * 3./4. * con_fac * dtinv * w_nat[2]);
        add_value(values, col_eb , face_fac * 1./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_es , face_fac * 1./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_ess, face_fac * 1./4. * con_fac * dtinv * w_nat[2]);

        // Contribution Delta q
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_s  + 1, 0.0);
        add_value(values, col_ss + 1, 0.0);

        // Contribution Delta r
        // face 0
        face_fac = 0.5;
        add_value(values, col_wb  + 2, face_fac * 1./4. * con_fac *  theta);
        add_value(values, col_ws  + 2, face_fac * 1./4. * con_fac * -theta);
        add_value(values, col_wss + 2, face_fac * 1./4. * 0.0);
        add_value(values, col_b   + 2, face_fac * 3./4. * con_fac *  theta);
        add_value(values, col_s   + 2, face_fac * 3./4. * con_fac * -theta);
        add_value(values, col_ss  + 2, face_fac * 3./4. * 0.0);

        // face 1
        face_fac = 0.5;
        add_value(values, col_b   + 2, face_fac * 3./4. * con_fac *  theta);
        add_value(values, col_s   + 2, face_fac * 3./4. * con_fac * -theta);
        add_value(values, col_ss  + 2, face_fac * 3./4. * 0.0);
        add_value(values, col_eb  + 2, face_fac * 1./4. * con_fac *  theta);
        add_value(values, col_es  + 2, face_fac * 1./4. * con_fac * -theta);
        add_value(values, col_ess + 2, face_fac * 1./4. * 0.0);

        //
        double dhdt_w = dtinv * (hp_w - hn_w);
        double dhdt_b = dtinv * (hp_b - hn_b);
        double dhdt_e = dtinv * (hp_e - hn_e);
        double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_w);
        double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_e);
        double drdy_0 = 0.25 * (3.0 * drdy_b + 1.0 * drdy_w);
        double drdy_1 = 0.25 * (3.0 * drdy_b + 1.0 * drdy_e);

        rhs[row + 2] += - ( 0.5 * con_fac * (dy_deta_0 * dhdt_0 + drdy_0) 
                          + 0.5 * con_fac * (dy_deta_1 * dhdt_1 + drdy_1) );
    }
    return 0;
}
//==============================================================================
int boundary_east(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_EAST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    std::fill_n(values + c_eq, 3 * 27, 0.0);  // set all coefficients for one row of Delta c-, Delta q- and Delta r-equation to zero

    int p_7 = c_eq/(3*27);  // node number of boundary point, ie east point of molecule
    int p_8 = p_7 + 1;
    int p_6 = p_7 - 1;
    int p_5 = p_7 - ny + 1;
    int p_4 = p_7 - ny;
    int p_3 = p_7 - ny - 1;
    int p_2 = p_7 - 2 * ny + 1;
    int p_1 = p_7 - 2 * ny;
    int p_0 = p_7 - 2 * ny - 1;

    double htheta_b = w_nat[0] * htheta[p_7] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_1];
    double zb_b = w_nat[0] * zb[p_7] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_1];

    double sign = 1.0;
    double h_given = bc[BC_EAST] - zb_b;
    double h_infty = -zb_b;
    double c_wave = std::sqrt(g * htheta_b);
    double con_fac = c_wave;

    // eeast
    double dx_dxi_n = x[p_8] - x[p_5];
    double dx_dxi_b = x[p_7] - x[p_4];
    double dx_dxi_s = x[p_6] - x[p_3];
    double dy_deta_0 = 0.5 * (y[p_5] - y[p_4]) + 0.5 * (y[p_8] - y[p_7]);
    double dy_deta_1 = 0.5 * (y[p_4] - y[p_3]) + 0.5 * (y[p_7] - y[p_6]);
    double dx_dxi_0 = 0.25 * (3. * (x[p_7] - x[p_4]) + 1. * (x[p_8] - x[p_5]));
    double dx_dxi_1 = 0.25 * (3. * (x[p_7] - x[p_4]) + 1. * (x[p_6] - x[p_3]));

    if (bc_type[BC_EAST] == "no_slip" || bc_type[BC_EAST] == "free_slip")
    {
        if (bc_type[BC_EAST] == "no_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 7 * 3;
            size_t col_w  = c_eq + 4 * 3;
            size_t col_ww = c_eq + 1 * 3;
            values[col_b ] =  1.0 * theta;
            values[col_w ] = -1.0 * theta;
            values[col_ww] = 0.0;
            rhs[row] = -(htheta[p_7] + zb[p_7] - htheta[p_4] - zb[p_4]);

            // Contribution Delta q
            col_b  = q_eq + 7 * 3;
            col_w  = q_eq + 4 * 3;
            col_ww = q_eq + 1 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_w  + 1] = 0.5 * theta;
            values[col_ww + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 7 * 3;
            col_w  = r_eq + 4 * 3;
            col_ww = r_eq + 1 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_w  + 2]  = 0.5 * theta;
            values[col_ww + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
        else if (bc_type[BC_EAST] == "free_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 7 * 3;
            size_t col_w  = c_eq + 4 * 3;
            size_t col_ww = c_eq + 1 * 3;
            values[col_b ] =  1.0 * theta;
            values[col_w ] = -1.0 * theta;
            values[col_ww] = 0.0;
            rhs[row] = -(htheta[p_7] + zb[p_7] - htheta[p_4] - zb[p_4]);

            // Contribution Delta q
            col_b  = q_eq + 7 * 3;
            col_w  = q_eq + 4 * 3;
            col_ww = q_eq + 1 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_w  + 1] = 0.5 * theta;
            values[col_ww + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 7 * 3;
            col_w  = r_eq + 4 * 3;
            col_ww = r_eq + 1 * 3;
            values[col_b  + 2]  =  1.0 * theta;
            values[col_w  + 2]  = -1.0 * theta;
            values[col_ww + 2]  = 0.0;
            rhs[row + 2] = -(rtheta[p_7] - rtheta[p_4]);
        }
    }
    else
    {
        if (bc_type[BC_EAST] == "mooiman")
        {
            // essential boundary condition
            // ----------------------------
            // first row
            size_t col_nb  = c_eq + 8 * 3;
            size_t col_nw  = c_eq + 5 * 3;
            size_t col_nww = c_eq + 2 * 3;
            size_t col_b  = c_eq + 7 * 3;
            size_t col_w  = c_eq + 4 * 3;
            size_t col_ww = c_eq + 1 * 3;
            size_t col_sb  = c_eq + 6 * 3;
            size_t col_sw  = c_eq + 3 * 3;
            size_t col_sww = c_eq + 0 * 3;
            //
            // face_0 
            double face_fac = 0.5 * dy_deta_0;
            add_value(values, col_nb , face_fac * 1./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_nw , face_fac * 1./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_nww, face_fac * 1./4. * theta * w_ess[2] * -c_wave);
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_w  , face_fac * 3./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_ww , face_fac * 3./4. * theta * w_ess[2] * -c_wave);
            
             // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_w  , face_fac * 3./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_ww , face_fac * 3./4. * theta * w_ess[2] * -c_wave);
            add_value(values, col_sb , face_fac * 1./4. * theta * w_ess[0] * -c_wave);
            add_value(values, col_sw , face_fac * 1./4. * theta * w_ess[1] * -c_wave);
            add_value(values, col_sww, face_fac * 1./4. * theta * w_ess[2] * -c_wave);
            
            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dy_deta_0;
            add_value(values, col_nb  + 1, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_nw  + 1, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_nww + 1, face_fac * 1./4. * theta * w_ess[2]);
            add_value(values, col_b   + 1, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_w   + 1, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ww  + 1, face_fac * 3./4. * theta * w_ess[2]);
            
            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_w   + 1, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ww  + 1, face_fac * 3./4. * theta * w_ess[2]);
            add_value(values, col_sb  + 1, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_sw  + 1, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_sww + 1, face_fac * 1./4. * theta * w_ess[2]);
            
            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_w  + 2, 0.0);
            add_value(values, col_ww + 2, 0.0);
            //
            // Right hand side
            face_fac = 0.5;
            double htheta_b   = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_8]);
            double htheta_im1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_5]);
            double htheta_im2 = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_2]);
            double ht_im12 = face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2);
            htheta_b   = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_6]);
            htheta_im1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_3]);
            htheta_im2 = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_0]);
            ht_im12 += face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2);
            
            double qtheta_b   = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_8]);
            double qtheta_im1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_5]);
            double qtheta_im2 = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_2]);
            double qt_im12 = face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_im1 + w_ess[2] * qtheta_im2);
            qtheta_b   = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_6]);
            qtheta_im1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_3]);
            qtheta_im2 = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_0]);
            qt_im12 += face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_im1 + w_ess[2] * qtheta_im2);
            
            rhs[row] = -0.5 * (dy_deta_0 + dy_deta_1) * (qt_im12 - c_wave * (ht_im12 - h_given));
            if (bc_vars[BC_EAST] == "zeta")
            {
                rhs[row] += -0.5 * (dy_deta_0 + dy_deta_1) * 2. * c_wave * bc[BC_EAST];
            }
            if (bc_vars[BC_EAST] == "q")
            {
                rhs[row] += -0.5 * ( dy_deta_0 + dy_deta_1) * 2. * bc[BC_EAST];
            }
        }
        else if (bc_type[BC_EAST] == "borsboom")
        {
            //----------------------------------------------------------------------
            // Essential boundary condition
            //----------------------------------------------------------------------
            double hp_n = w_ess[0] * hp[p_8] + w_ess[1] * hp[p_5] + w_ess[2] * hp[p_2];
            double hp_b = w_ess[0] * hp[p_7] + w_ess[1] * hp[p_4] + w_ess[2] * hp[p_1];
            double hp_s = w_ess[0] * hp[p_6] + w_ess[1] * hp[p_3] + w_ess[2] * hp[p_0];
        
            double hn_n = w_ess[0] * hn[p_8] + w_ess[1] * hn[p_5] + w_ess[2] * hn[p_2];
            double hn_b = w_ess[0] * hn[p_7] + w_ess[1] * hn[p_4] + w_ess[2] * hn[p_1];
            double hn_s = w_ess[0] * hn[p_6] + w_ess[1] * hn[p_3] + w_ess[2] * hn[p_0];

            double htheta_n = w_ess[0] * htheta[p_8] + w_ess[1] * htheta[p_5] + w_ess[2] * htheta[p_2];  // north of boundary location
            double htheta_b = w_ess[0] * htheta[p_7] + w_ess[1] * htheta[p_4] + w_ess[2] * htheta[p_1];
            double htheta_s = w_ess[0] * htheta[p_6] + w_ess[1] * htheta[p_3] + w_ess[2] * htheta[p_0];  // south of boundary location

            double qp_n = w_ess[0] * qp[p_8] + w_ess[1] * qp[p_5] + w_ess[2] * qp[p_2];
            double qp_b = w_ess[0] * qp[p_7] + w_ess[1] * qp[p_4] + w_ess[2] * qp[p_1];
            double qp_s = w_ess[0] * qp[p_6] + w_ess[1] * qp[p_3] + w_ess[2] * qp[p_0];
        
            double qn_n = w_ess[0] * qn[p_8] + w_ess[1] * qn[p_5] + w_ess[2] * qn[p_2];
            double qn_b = w_ess[0] * qn[p_7] + w_ess[1] * qn[p_4] + w_ess[2] * qn[p_1];
            double qn_s = w_ess[0] * qn[p_6] + w_ess[1] * qn[p_3] + w_ess[2] * qn[p_0];
        
            double qtheta_n = w_ess[0] * qtheta[p_8] + w_ess[1] * qtheta[p_5] + w_ess[2] * qtheta[p_2];  // north of boundary location
            double qtheta_b = w_ess[0] * qtheta[p_7] + w_ess[1] * qtheta[p_4] + w_ess[2] * qtheta[p_1];
            double qtheta_s = w_ess[0] * qtheta[p_6] + w_ess[1] * qtheta[p_3] + w_ess[2] * qtheta[p_0];  // south of boundary location
            //------------------------------------------------------------------
            // first row
            //
            size_t col_nb  = c_eq + 8 * 3;
            size_t col_nw  = c_eq + 5 * 3;
            size_t col_nww = c_eq + 2 * 3;
            size_t col_b   = c_eq + 7 * 3;
            size_t col_w   = c_eq + 4 * 3;
            size_t col_ww  = c_eq + 1 * 3;
            size_t col_sb  = c_eq + 6 * 3;
            size_t col_sw  = c_eq + 3 * 3;
            size_t col_sww = c_eq + 0 * 3;
            //
            if (do_convection) { con_fac = c_wave + qp_b / hp_b; }
            //
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_nb , face_fac * 1./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_nw , face_fac * 1./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_nww, face_fac * 1./4. * dtinv * -con_fac * w_ess[2]);
            add_value(values, col_b  , face_fac * 3./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_w  , face_fac * 3./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_ww , face_fac * 3./4. * dtinv * -con_fac * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_w  , face_fac * 3./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_ww , face_fac * 3./4. * dtinv * -con_fac * w_ess[2]);
            add_value(values, col_sb , face_fac * 1./4. * dtinv * -con_fac * w_ess[0]);
            add_value(values, col_sw , face_fac * 1./4. * dtinv * -con_fac * w_ess[1]);
            add_value(values, col_sww, face_fac * 1./4. * dtinv * -con_fac * w_ess[2]);

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_nb  + 1, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_nw  + 1, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_nww + 1, face_fac * 1./4. * dtinv * w_ess[2]);
            add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_w   + 1, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ww  + 1, face_fac * 3./4. * dtinv * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_w   + 1, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ww  + 1, face_fac * 3./4. * dtinv * w_ess[2]);
            add_value(values, col_sb  + 1, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_sw  + 1, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_sww + 1, face_fac * 1./4. * dtinv * w_ess[2]);

            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_w  + 2, 0.0);
            add_value(values, col_ww + 2, 0.0);
            //
            double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_n);
            double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_s);
            double dhdt_n = dtinv * (hp_n - hn_n);
            double dhdt_b = dtinv * (hp_b - hn_b);
            double dhdt_s = dtinv * (hp_s - hn_s);
            double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_n);
            double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_s);
            double dqdt_n = dtinv * (qp_n - qn_n);
            double dqdt_b = dtinv * (qp_b - qn_b);
            double dqdt_s = dtinv * (qp_s - qn_s);
            double dqdt_0 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_n);
            double dqdt_1 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_s);

            rhs[row] = - ( 0.5 * dx_dxi_0 * (dqdt_0 - con_fac * dhdt_0) 
                         + 0.5 * dx_dxi_1 * (dqdt_1 - con_fac * dhdt_1) 
                         );

            double corr_term = 0.0;
            if (bc_vars[BC_EAST] == "zeta")
            {
                // face 0
                face_fac = 0.5 * dx_dxi_0;
                add_value(values, col_nb , face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_nw , face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nww, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b  , face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_w  , face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ww , face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dx_dxi_1;
                add_value(values, col_b  , face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_w  , face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ww , face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_sb , face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_sw , face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_sww, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                double zb_n = w_ess[0] * zb[p_8] + w_ess[1] * zb[p_5] + w_ess[2] * zb[p_2];
                double zb_b = w_ess[0] * zb[p_7] + w_ess[1] * zb[p_4] + w_ess[2] * zb[p_1];
                double zb_s = w_ess[0] * zb[p_6] + w_ess[1] * zb[p_3] + w_ess[2] * zb[p_0];
                double zb_0 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_n);
                double zb_1 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_s);

                corr_term = - ( 0.5 * dx_dxi_0 * ( dhdt_0 + ( eps_bc_corr * ( (bc[BC_EAST] - zb_0) - htheta_0)) ) 
                              + 0.5 * dx_dxi_1 * ( dhdt_1 + ( eps_bc_corr * ( (bc[BC_EAST] - zb_1) - htheta_1)) ) 
                              );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_EAST] == "q")
            {
                // face 0
                face_fac = 0.5 * dx_dxi_0;
                add_value(values, col_nb  + 1, face_fac * 1./4. * (dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_nw  + 1, face_fac * 1./4. * (dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nww + 1, face_fac * 1./4. * (dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b   + 1, face_fac * 3./4. * (dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_w   + 1, face_fac * 3./4. * (dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ww  + 1, face_fac * 3./4. * (dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dx_dxi_1;
                add_value(values, col_b   + 1, face_fac * 3./4. * (dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_w   + 1, face_fac * 3./4. * (dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ww  + 1, face_fac * 3./4. * (dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_sb  + 1, face_fac * 1./4. * (dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_sw  + 1, face_fac * 1./4. * (dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_sww + 1, face_fac * 1./4. * (dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                double qtheta_0 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_n);
                double qtheta_1 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_s);

                corr_term = ( 0.5 * dx_dxi_0 * (- dqdt_0 + sign * eps_bc_corr * (bc[BC_EAST] - qtheta_0))
                            + 0.5 * dx_dxi_1 * (- dqdt_1 + sign * eps_bc_corr * (bc[BC_EAST] - qtheta_1))
                            );
                rhs[row] += corr_term;
            }
        }
        //----------------------------------------------------------------------
        // natural boundary condition
        // ---------------------------------------------------------------------
        double hp_n = w_nat[0] * hp[p_8] + w_nat[1] * hp[p_5] + w_nat[2] * hp[p_2];
        double hp_b = w_nat[0] * hp[p_7] + w_nat[1] * hp[p_4] + w_nat[2] * hp[p_1];
        double hp_s = w_nat[0] * hp[p_6] + w_nat[1] * hp[p_3] + w_nat[2] * hp[p_0];
        
        double hn_n = w_nat[0] * hn[p_8] + w_nat[1] * hn[p_5] + w_nat[2] * hn[p_2];
        double hn_b = w_nat[0] * hn[p_7] + w_nat[1] * hn[p_4] + w_nat[2] * hn[p_1];
        double hn_s = w_nat[0] * hn[p_6] + w_nat[1] * hn[p_3] + w_nat[2] * hn[p_0];

        double dqdx_n = qtheta[p_8] - qtheta[p_5];
        double dqdx_b = qtheta[p_7] - qtheta[p_4];
        double dqdx_s = qtheta[p_6] - qtheta[p_3];

        double qp_n = w_nat[0] * qp[p_8] + w_nat[1] * qp[p_5] + w_nat[2] * qp[p_2];
        double qp_b = w_nat[0] * qp[p_7] + w_nat[1] * qp[p_4] + w_nat[2] * qp[p_1];
        double qp_s = w_nat[0] * qp[p_6] + w_nat[1] * qp[p_3] + w_nat[2] * qp[p_0];
        
        double qn_n = w_nat[0] * qn[p_8] + w_nat[1] * qn[p_5] + w_nat[2] * qn[p_2];
        double qn_b = w_nat[0] * qn[p_7] + w_nat[1] * qn[p_4] + w_nat[2] * qn[p_1];
        double qn_s = w_nat[0] * qn[p_6] + w_nat[1] * qn[p_3] + w_nat[2] * qn[p_0];
        
        double htheta_n = w_nat[0] * htheta[p_8] + w_nat[1] * htheta[p_5] + w_nat[2] * htheta[p_2];  // north of boundary location
        double htheta_b = w_nat[0] * htheta[p_7] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_1];
        double htheta_s = w_nat[0] * htheta[p_6] + w_nat[1] * htheta[p_3] + w_nat[2] * htheta[p_0];  // south of boundary location

        double qtheta_n = w_nat[0] * qtheta[p_8] + w_nat[1] * qtheta[p_5] + w_nat[2] * qtheta[p_2];  // north of boundary location
        double qtheta_b = w_nat[0] * qtheta[p_7] + w_nat[1] * qtheta[p_4] + w_nat[2] * qtheta[p_1];
        double qtheta_s = w_nat[0] * qtheta[p_6] + w_nat[1] * qtheta[p_3] + w_nat[2] * qtheta[p_0];  // south of boundary location

        double rtheta_n = w_nat[0] * rtheta[p_8] + w_nat[1] * rtheta[p_5] + w_nat[2] * rtheta[p_2];  // north of boundary location
        double rtheta_b = w_nat[0] * rtheta[p_7] + w_nat[1] * rtheta[p_4] + w_nat[2] * rtheta[p_1];
        double rtheta_s = w_nat[0] * rtheta[p_6] + w_nat[1] * rtheta[p_3] + w_nat[2] * rtheta[p_0];  // south of boundary location

        double dzetadx_n = htheta[p_8] + zb[p_8] - htheta[p_5] - zb[p_5];
        double dzetadx_b = htheta[p_7] + zb[p_7] - htheta[p_4] - zb[p_4];
        double dzetadx_s = htheta[p_6] + zb[p_6] - htheta[p_3] - zb[p_3];
        // ---------------------------------------------------------------------
        if (do_convection) { con_fac = c_wave - qp_b / hp_b; }
        // ---------------------------------------------------------------------
        // second row
        // momentum part dq/dt + gh d(zeta)/dx
        //
        size_t col_nb  = q_eq + 8 * 3;
        size_t col_nw  = q_eq + 5 * 3;
        size_t col_nww = q_eq + 2 * 3;
        size_t col_b   = q_eq + 7 * 3;
        size_t col_w   = q_eq + 4 * 3;
        size_t col_ww  = q_eq + 1 * 3;
        size_t col_sb  = q_eq + 6 * 3;
        size_t col_sw  = q_eq + 3 * 3;
        size_t col_sww = q_eq + 0 * 3;
        //
        // Contribution Delta h
        // face 0
        double face_fac = 0.5;
        add_value(values, col_nb , face_fac * 1./4. * w_nat[0] * theta * g * dzetadx_n + face_fac * 1./4. * theta * g * htheta_n);
        add_value(values, col_nw , face_fac * 1./4. * w_nat[1] * theta * g * dzetadx_n - face_fac * 1./4. * theta * g * htheta_n);
        add_value(values, col_nww, face_fac * 1./4. * w_nat[2] * theta * g * dzetadx_n);
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetadx_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_w  , face_fac * 3./4. * w_nat[1] * theta * g * dzetadx_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ww , face_fac * 3./4. * w_nat[2] * theta * g * dzetadx_b);

        //face 1
        face_fac = 0.5;
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetadx_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_w  , face_fac * 3./4. * w_nat[1] * theta * g * dzetadx_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ww , face_fac * 3./4. * w_nat[2] * theta * g * dzetadx_b);
        add_value(values, col_sb , face_fac * 1./4. * w_nat[0] * theta * g * dzetadx_s + face_fac * 1./4. * theta * g * htheta_s);
        add_value(values, col_sw , face_fac * 1./4. * w_nat[1] * theta * g * dzetadx_s - face_fac * 1./4. * theta * g * htheta_s);
        add_value(values, col_sww, face_fac * 1./4. * w_nat[2] * theta * g * dzetadx_s);

        // Contribution Delta q
        // face 0
        face_fac = 0.5 * dy_deta_0;
        add_value(values, col_nb  + 1, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_nw  + 1, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_nww + 1, face_fac * 1./4. * dtinv * w_nat[2]);
        add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_w   + 1, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ww  + 1, face_fac * 3./4. * dtinv * w_nat[2]);
        
        //face 1
        face_fac = 0.5 * dy_deta_1;
        add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_w   + 1, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ww  + 1, face_fac * 3./4. * dtinv * w_nat[2]);
        add_value(values, col_sb  + 1, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_sw  + 1, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_sww + 1, face_fac * 1./4. * dtinv * w_nat[2]); 

        // Contribution Delta r
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_w  + 2, 0.0);
        add_value(values, col_ww + 2, 0.0);
        //
        double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_n);
        double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_s);
        double qtheta_0 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_n);
        double qtheta_1 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_s);
        double dzetadx_0 = 0.25 * (3.0 * dzetadx_b + 1.0 * dzetadx_n);
        double dzetadx_1 = 0.25 * (3.0 * dzetadx_b + 1.0 * dzetadx_s);
        double dqdt_n = dtinv * (qp_n - qn_n);
        double dqdt_b = dtinv * (qp_b - qn_b);
        double dqdt_s = dtinv * (qp_s - qn_s);
        double dqdt_0 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_n);
        double dqdt_1 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_s);

        rhs[row + 1] = - ( 0.5 * ( dy_deta_0 * dqdt_0 + g * htheta_0 * dzetadx_0 ) 
                         + 0.5 * ( dy_deta_1 * dqdt_1 + g * htheta_1 * dzetadx_1 ) );
        //
        if (do_convection)
        {
            // East boundary convection

            double dh_dxi_0 = 0.25 * (3.0 * (htheta[p_7] - htheta[p_4]) +  1.0 * (htheta[p_8] - htheta[p_5]));
            double dh_dxi_1 = 0.25 * (3.0 * (htheta[p_7] - htheta[p_4]) +  1.0 * (htheta[p_6] - htheta[p_3]));
            double dq_dxi_0 = 0.25 * (3.0 * (qtheta[p_7] - qtheta[p_4]) +  1.0 * (qtheta[p_8] - qtheta[p_5]));
            double dq_dxi_1 = 0.25 * (3.0 * (qtheta[p_7] - qtheta[p_4]) +  1.0 * (qtheta[p_6] - qtheta[p_3]));

            double dh_deta_0 = 0.5 * ((htheta[p_5] - htheta[p_4]) + (htheta[p_8] - htheta[p_7]));
            double dh_deta_1 = 0.5 * ((htheta[p_4] - htheta[p_3]) + (htheta[p_7] - htheta[p_6]));
            double dq_deta_0 = 0.5 * ((qtheta[p_5] - qtheta[p_4]) + (qtheta[p_8] - qtheta[p_7]));
            double dq_deta_1 = 0.5 * ((qtheta[p_4] - qtheta[p_3]) + (qtheta[p_7] - qtheta[p_6]));

            double dy_dxi_0 = dcdx_scvf_n(y[p_7], y[p_4], y[p_8], y[p_5]);
            double dy_dxi_1 = dcdx_scvf_n(y[p_7], y[p_4], y[p_6], y[p_3]);

            double aa_0 = - 2. * qtheta_0 / (htheta_0 * htheta_0) * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) + 2. * (qtheta_0 * qtheta_0) / (htheta_0 * htheta_0 * htheta_0) * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0);
            double aa_1 = - 2. * qtheta_1 / (htheta_1 * htheta_1) * (dy_deta_1 * dq_dxi_1 - dy_dxi_1 * dq_deta_1) + 2. * (qtheta_1 * qtheta_1) / (htheta_1 * htheta_1 * htheta_1) * (dy_deta_1 * dh_dxi_1 - dy_dxi_1 * dh_deta_1);

            double bb_0 = 2. / htheta_0 * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) - 2. * qtheta_0 / (htheta_0 * htheta_0) * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0);
            double bb_1 = 2. / htheta_1 * (dy_deta_1 * dq_dxi_1 - dy_dxi_1 * dq_deta_1) - 2. * qtheta_1 / (htheta_1 * htheta_1) * (dy_deta_1 * dh_dxi_1 - dy_dxi_1 * dh_deta_1);

            double cc_0 = - qtheta_0 * qtheta_0 / (htheta_0 * htheta_0);
            double cc_1 = - qtheta_1 * qtheta_1 / (htheta_1 * htheta_1);

            double dd_0 = 2. * qtheta_0 / htheta_0;
            double dd_1 = 2. * qtheta_1 / htheta_1;

            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * theta / dx_dxi_0;
            add_value(values, col_nb , face_fac * 1./4. * (aa_0 * w_nat[0] + dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_nw , face_fac * 1./4. * (aa_0 * w_nat[1] - dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_nww, face_fac * 1./4. * (aa_0 * w_nat[2]));
            add_value(values, col_b  , face_fac * 3./4. * (aa_0 * w_nat[0] + dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_w  , face_fac * 3./4. * (aa_0 * w_nat[1] - dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_ww , face_fac * 3./4. * (aa_0 * w_nat[2]));

            //face 1
            face_fac = 0.5 * theta / dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * (aa_1 * w_nat[0] + dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_w  , face_fac * 3./4. * (aa_1 * w_nat[1] - dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ww , face_fac * 3./4. * (aa_1 * w_nat[2]));
            add_value(values, col_sb , face_fac * 1./4. * (aa_1 * w_nat[0] + dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_sw , face_fac * 1./4. * (aa_1 * w_nat[1] - dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_sww, face_fac * 1./4. * (aa_1 * w_nat[2]));

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * theta / dx_dxi_0;
            add_value(values, col_nb  + 1, face_fac * 1./4. * (bb_0 * w_nat[0] + dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_nw  + 1, face_fac * 1./4. * (bb_0 * w_nat[1] - dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_nww + 1, face_fac * 1./4. * (bb_0 * w_nat[2]));
            add_value(values, col_b   + 1, face_fac * 3./4. * (bb_0 * w_nat[0] + dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_w   + 1, face_fac * 3./4. * (bb_0 * w_nat[1] - dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_ww  + 1, face_fac * 3./4. * (bb_0 * w_nat[2]));

            // face 1
            face_fac = 0.5 * theta / dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * (bb_1 * w_nat[0] + dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_w   + 1, face_fac * 3./4. * (bb_1 * w_nat[1] - dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ww  + 1, face_fac * 3./4. * (bb_1 * w_nat[2]));
            add_value(values, col_sb  + 1, face_fac * 1./4. * (bb_1 * w_nat[0] + dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_sw  + 1, face_fac * 1./4. * (bb_1 * w_nat[1] - dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_sww + 1, face_fac * 1./4. * (bb_1 * w_nat[2]));

            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_w  + 2, 0.0);
            add_value(values, col_ww + 2, 0.0);

            rhs[row + 1] += - (
                0.5 * dd_0 * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) / dx_dxi_0
              + 0.5 * cc_0 * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0) / dx_dxi_0
              + 0.5 * dd_1 * (dy_deta_1 * dq_dxi_1 - dy_dxi_1 * dq_deta_1) / dx_dxi_1
              + 0.5 * cc_1 * (dy_deta_1 * dh_dxi_1 - dy_dxi_1 * dh_deta_1) / dx_dxi_1
                );
        }
        if (do_bed_shear_stress)
        {
            // East boundary bed shear stress

            double cf_n = cf;
            double cf_b = cf;
            double cf_s = cf;
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_nb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nw , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nww, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_w  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ww , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_w  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ww , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_sb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_sw , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_sww, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_nb  + 1, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nw  + 1, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nww + 1, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );
            add_value(values, col_b   + 1, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_w   + 1, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ww  + 1, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_w   + 1, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ww  + 1, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_sb  + 1, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_sw  + 1, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_sww + 1, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    

            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_w  + 2, 0.0);
            add_value(values, col_ww + 2, 0.0);

            // right hand side
            double rtheta_0 = 0.0;  // no tangential discharge r
            double rtheta_1 = 0.0;
            double abs_qtheta_0 = abs_vecq(qtheta_0, rtheta_0, 1.0);
            double abs_qtheta_1 = abs_vecq(qtheta_1, rtheta_1, 1.0);

            double cf_0 = 0.25 * (3.0 * cf_b + 1.0 * cf_n);
            double cf_1 = 0.25 * (3.0 * cf_b + 1.0 * cf_s);

            rhs[row + 1] += - (
                  0.5 * dx_dxi_0 * cf_0 * qtheta_0 * abs_qtheta_0 / (htheta_0 * htheta_0)
                + 0.5 * dx_dxi_1 * cf_1 * qtheta_1 * abs_qtheta_1 / (htheta_1 * htheta_1)
                );
        }  
        do_viscosity = false;  // east boundary
        if (do_viscosity)
        {
            // East boundary viscosity
        }
        //
        // continuity part +c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        // face 0
        face_fac = 0.5 * dx_dxi_0;
        add_value(values, col_nb , face_fac * 1./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_nw , face_fac * 1./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_nww, face_fac * 1./4. * con_fac * dtinv * w_nat[2]);
        add_value(values, col_b  , face_fac * 3./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_w  , face_fac * 3./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_ww , face_fac * 3./4. * con_fac * dtinv * w_nat[2]);

        // face 1
        face_fac = 0.5 * dx_dxi_1;
        add_value(values, col_b  , face_fac * 3./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_w  , face_fac * 3./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_ww , face_fac * 3./4. * con_fac * dtinv * w_nat[2]);
        add_value(values, col_sb , face_fac * 1./4. * con_fac * dtinv * w_nat[0]);
        add_value(values, col_sw , face_fac * 1./4. * con_fac * dtinv * w_nat[1]);
        add_value(values, col_sww, face_fac * 1./4. * con_fac * dtinv * w_nat[2]);

        // Contribution Delta q
        // face 0
        face_fac = 0.5;
        add_value(values, col_nb  + 1, face_fac * 1./4. * con_fac *  theta);
        add_value(values, col_nw  + 1, face_fac * 1./4. * con_fac * -theta);
        add_value(values, col_nww + 1, face_fac * 1./4. * 0.0);
        add_value(values, col_b   + 1, face_fac * 3./4. * con_fac *  theta);
        add_value(values, col_w   + 1, face_fac * 3./4. * con_fac * -theta);
        add_value(values, col_ww  + 1, face_fac * 3./4. * 0.0);

        // face 1
        face_fac = 0.5;
        add_value(values, col_b   + 1, face_fac * 3./4. * con_fac *  theta);
        add_value(values, col_w   + 1, face_fac * 3./4. * con_fac * -theta);
        add_value(values, col_ww  + 1, face_fac * 3./4. * 0.0);
        add_value(values, col_sb  + 1, face_fac * 1./4. * con_fac *  theta);
        add_value(values, col_sw  + 1, face_fac * 1./4. * con_fac * -theta);
        add_value(values, col_sww + 1, face_fac * 1./4. * 0.0);

        // Contribution Delta r
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_w  + 2, 0.0);
        add_value(values, col_ww + 2, 0.0);
        //

        double dhdt_n = dtinv * (hp_n - hn_n);
        double dhdt_b = dtinv * (hp_b - hn_b);
        double dhdt_s = dtinv * (hp_s - hn_s);
        double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_n);
        double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_s);
        double dqdx_0 = 0.25 * (3.0 * dqdx_b + 1.0 * dqdx_n);
        double dqdx_1 = 0.25 * (3.0 * dqdx_b + 1.0 * dqdx_s);

        rhs[row + 1] += - ( 0.5 * con_fac * (dx_dxi_0 * dhdt_0 + dqdx_0) 
                          + 0.5 * con_fac * (dx_dxi_1 * dhdt_1 + dqdx_1) );
        //
        // -----------------------------------------------------------------
        // third row
        // r-momentum (tangential equation, r == 0)
        //
        col_nb  = r_eq + 8 * 3;
        col_nw  = r_eq + 5 * 3;
        col_nww = r_eq + 2 * 3;
        col_b   = r_eq + 7 * 3;
        col_w   = r_eq + 4 * 3;
        col_ww  = r_eq + 1 * 3;
        col_sb  = r_eq + 6 * 3;
        col_sw  = r_eq + 3 * 3;
        col_sww = r_eq + 0 * 3;
        // 
        // Contribution Delta h
        add_value(values, col_b , 0.0);
        add_value(values, col_w , 0.0);
        add_value(values, col_ww, 0.0);
        // Contribution Delta q
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_w  + 1, 0.0);
        add_value(values, col_ww + 1, 0.0);
        // Contribution Delta r
        add_value(values, col_b  + 2, 1.0);
        add_value(values, col_w  + 2, 0.0);
        add_value(values, col_ww + 2, 0.0);
        //
        rhs[row + 2] = 0.0;
    }
    return 0;
}
//==============================================================================
int boundary_south(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf, 
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_SOUTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    std::fill_n(values + c_eq, 3 * 27, 0.0);  // set all coefficients for one row of Delta c-, Delta q- and Delta r-equation to zero

    int p_3 = c_eq/(3*27);  // node number of boundary point, ie south point of molecule
    int p_4 = p_3 + 1;
    int p_5 = p_3 + 2;
    int p_0 = p_3 - ny;
    int p_1 = p_3 - ny + 1;
    int p_2 = p_3 - ny + 2;
    int p_6 = p_3 + ny;
    int p_7 = p_3 + ny + 1;
    int p_8 = p_3 + ny + 2;

    double htheta_b = w_nat[0] * htheta[p_3] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_5];
    double zb_b = w_nat[0] * zb[p_3] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_5];

    double h_given = bc[BC_SOUTH] - zb_b;
    double h_infty = -zb_b;
    double c_wave = std::sqrt(g * htheta_b);
    double con_fac = c_wave;

    // ssouth
    double dy_deta_e = y[p_7] - y[p_6];
    double dy_deta_b = y[p_4] - y[p_3];
    double dy_deta_w = y[p_1] - y[p_0];
    double dx_dxi_0 = 0.5 * (x[p_5] - x[p_2]) + 0.5 * (x[p_4] - x[p_1]);
    double dx_dxi_1 = 0.5 * (x[p_8] - x[p_5]) + 0.5 * (x[p_7] - x[p_4]);
    double dy_deta_0 = 0.25 * (3. * (y[p_4] - y[p_3]) + 1. * (y[p_7] - y[p_6]));
    double dy_deta_1 = 0.25 * (3. * (y[p_4] - y[p_3]) + 1. * (y[p_1] - y[p_0]));

    if (bc_type[BC_SOUTH] == "no_slip" || bc_type[BC_SOUTH] == "free_slip")
    {
        if (bc_type[BC_SOUTH] == "no_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 3 * 3;
            size_t col_n  = c_eq + 4 * 3;
            size_t col_nn = c_eq + 5 * 3;
            values[col_b ] = -1.0 * theta;
            values[col_n ] =  1.0 * theta;
            values[col_nn] = 0.0;
            rhs[row] = -(htheta[p_4] + zb[p_4] - htheta[p_3] - zb[p_3]);

            // Contribution Delta q
            col_b  = q_eq + 3 * 3;
            col_n  = q_eq + 4 * 3;
            col_nn = q_eq + 5 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_n  + 1] = 0.5 * theta;
            values[col_nn + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 3 * 3;
            col_n  = r_eq + 4 * 3;
            col_nn = r_eq + 5 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_n  + 2]  = 0.5 * theta;
            values[col_nn + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
        else if (bc_type[BC_SOUTH] == "free_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 3 * 3;
            size_t col_n  = c_eq + 4 * 3;
            size_t col_nn = c_eq + 5 * 3;
            values[col_b ] = -1.0 * theta;
            values[col_n ] =  1.0 * theta;
            values[col_nn] = 0.0;
            rhs[row] = -(htheta[p_4] + zb[p_4] - htheta[p_3] - zb[p_3]);

            // Contribution Delta q
            col_b  = q_eq + 3 * 3;
            col_n  = q_eq + 4 * 3;
            col_nn = q_eq + 5 * 3;
            values[col_b  + 1] = -1.0 * theta;
            values[col_n  + 1] =  1.0 * theta;
            values[col_nn + 1] = 0.0;
            rhs[row + 1] = -(qtheta[p_4] - qtheta[p_3]);

            // Contribution Delta r
            col_b  = r_eq + 3 * 3;
            col_n  = r_eq + 4 * 3;
            col_nn = r_eq + 5 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_n  + 2]  = 0.5 * theta;
            values[col_nn + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
    }
    else
    {
        if (bc_type[BC_SOUTH] == "mooiman")
        {
            // Essential boundary condition
            // ----------------------------
            // first row
            size_t col_wb  = c_eq + 0 * 3;
            size_t col_wn  = c_eq + 1 * 3;
            size_t col_wnn = c_eq + 2 * 3;
            size_t col_b   = c_eq + 3 * 3;
            size_t col_n   = c_eq + 4 * 3;
            size_t col_nn  = c_eq + 5 * 3;
            size_t col_eb  = c_eq + 6 * 3;
            size_t col_en  = c_eq + 7 * 3;
            size_t col_enn = c_eq + 8 * 3;
            //
            // face 0 
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_eb , face_fac * 1./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_en , face_fac * 1./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_enn, face_fac * 1./4. * theta * w_ess[2] * c_wave);
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_n  , face_fac * 3./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_nn , face_fac * 3./4. * theta * w_ess[2] * c_wave);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_n  , face_fac * 3./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_nn , face_fac * 3./4. * theta * w_ess[2] * c_wave);
            add_value(values, col_wb , face_fac * 1./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_wn , face_fac * 1./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_wnn, face_fac * 1./4. * theta * w_ess[2] * c_wave);

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_n  + 1, 0.0);
            add_value(values, col_nn + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_eb  + 2, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_en  + 2, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_enn + 2, face_fac * 1./4. * theta * w_ess[2]);
            add_value(values, col_b   + 2, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_n   + 2, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_nn  + 2, face_fac * 3./4. * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_n   + 2, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_nn  + 2, face_fac * 3./4. * theta * w_ess[2]);
            add_value(values, col_wb  + 2, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_wn  + 2, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_wnn + 2, face_fac * 1./4. * theta * w_ess[2]);
            //
            // Right hand side
            face_fac = 0.5;
            double htheta_b   = 0.25 * (3.0 * htheta[p_3] + 1.0 * htheta[p_0]);
            double htheta_jp1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_1]);
            double htheta_jp2 = 0.25 * (3.0 * htheta[p_5] + 1.0 * htheta[p_2]);
            double ht_jp12 = face_fac *  (w_ess[0] * htheta_b + w_ess[1] * htheta_jp1 + w_ess[2] * htheta_jp2);
            htheta_b   = 0.25 * (3.0 * htheta[p_3] + 1.0 * htheta[p_6]);
            htheta_jp1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_7]);
            htheta_jp2 = 0.25 * (3.0 * htheta[p_5] + 1.0 * htheta[p_8]);
            ht_jp12 += face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_jp1 + w_ess[2] * htheta_jp2);

            double rtheta_b   = 0.25 * (3.0 * rtheta[p_3] + 1.0 * rtheta[p_0]);
            double rtheta_jp1 = 0.25 * (3.0 * rtheta[p_4] + 1.0 * rtheta[p_1]);
            double rtheta_jp2 = 0.25 * (3.0 * rtheta[p_5] + 1.0 * rtheta[p_2]);
            double rt_jp12 = face_fac * (w_ess[0] * rtheta_b + w_ess[1] * rtheta_jp1 + w_ess[2] * rtheta_jp2);
            rtheta_b   = 0.25 * (3.0 * rtheta[p_3] + 1.0 * rtheta[p_6]);
            rtheta_jp1 = 0.25 * (3.0 * rtheta[p_4] + 1.0 * rtheta[p_7]);
            rtheta_jp2 = 0.25 * (3.0 * rtheta[p_5] + 1.0 * rtheta[p_8]);
            rt_jp12 += face_fac * (w_ess[0] * rtheta_b + w_ess[1] * rtheta_jp1 + w_ess[2] * rtheta_jp2);

            rhs[row] = - 0.5 * (dx_dxi_0 + dx_dxi_1) * ( rt_jp12 + c_wave * (ht_jp12 - h_given) );

            if (bc_vars[BC_SOUTH] == "zeta")
            {
                rhs[row] += 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * c_wave * bc[BC_SOUTH];
            }
            if (bc_vars[BC_SOUTH] == "q")
            {
                rhs[row] += 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * bc[BC_SOUTH];
            }
        }
        else if (bc_type[BC_SOUTH] == "borsboom")
        {
            //----------------------------------------------------------------------
            // Essential boundary condition
            //----------------------------------------------------------------------
            double hp_e = w_ess[0] * hp[p_6] + w_ess[1] * hp[p_7] + w_ess[2] * hp[p_8];
            double hp_b = w_ess[0] * hp[p_3] + w_ess[1] * hp[p_4] + w_ess[2] * hp[p_5];
            double hp_w = w_ess[0] * hp[p_0] + w_ess[1] * hp[p_1] + w_ess[2] * hp[p_2];
        
            double hn_e = w_ess[0] * hn[p_6] + w_ess[1] * hn[p_7] + w_ess[2] * hn[p_8];
            double hn_b = w_ess[0] * hn[p_3] + w_ess[1] * hn[p_4] + w_ess[2] * hn[p_5];
            double hn_w = w_ess[0] * hn[p_0] + w_ess[1] * hn[p_1] + w_ess[2] * hn[p_2];

            double htheta_e = w_ess[0] * htheta[p_6] + w_ess[1] * htheta[p_7] + w_ess[2] * htheta[p_8];  // east of boundary location
            double htheta_b = w_ess[0] * htheta[p_3] + w_ess[1] * htheta[p_4] + w_ess[2] * htheta[p_5];
            double htheta_w = w_ess[0] * htheta[p_0] + w_ess[1] * htheta[p_1] + w_ess[2] * htheta[p_2];  // west of boundary location

            double rp_e = w_ess[0] * rp[p_6] + w_ess[1] * rp[p_7] + w_ess[2] * rp[p_8];
            double rp_b = w_ess[0] * rp[p_3] + w_ess[1] * rp[p_4] + w_ess[2] * rp[p_5];
            double rp_w = w_ess[0] * rp[p_0] + w_ess[1] * rp[p_1] + w_ess[2] * rp[p_2];
        
            double rn_e = w_ess[0] * rn[p_6] + w_ess[1] * rn[p_7] + w_ess[2] * rn[p_8];
            double rn_b = w_ess[0] * rn[p_3] + w_ess[1] * rn[p_4] + w_ess[2] * rn[p_5];
            double rn_w = w_ess[0] * rn[p_0] + w_ess[1] * rn[p_1] + w_ess[2] * rn[p_2];

            double rtheta_e = w_ess[0] * rtheta[p_6] + w_ess[1] * rtheta[p_7] + w_ess[2] * rtheta[p_8];  // east of boundary location
            double rtheta_b = w_ess[0] * rtheta[p_3] + w_ess[1] * rtheta[p_4] + w_ess[2] * rtheta[p_5];
            double rtheta_w = w_ess[0] * rtheta[p_0] + w_ess[1] * rtheta[p_1] + w_ess[2] * rtheta[p_2];  // west of boundary location
            //------------------------------------------------------------------
            // first row
            //
            size_t col_wb  = c_eq + 0 * 3;
            size_t col_wn  = c_eq + 1 * 3;
            size_t col_wnn = c_eq + 2 * 3;
            size_t col_b   = c_eq + 3 * 3;
            size_t col_n   = c_eq + 4 * 3;
            size_t col_nn  = c_eq + 5 * 3;
            size_t col_eb  = c_eq + 6 * 3;
            size_t col_en  = c_eq + 7 * 3;
            size_t col_enn = c_eq + 8 * 3;
            //
            if (do_convection) { con_fac = c_wave - rp_b / hp_b; }
            //
            // Contribution h
            // face 0
            double face_fac = 0.5 * dy_deta_0;
            add_value(values, col_eb , face_fac * 1./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_en , face_fac * 1./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_enn, face_fac * 1./4. * dtinv * con_fac * w_ess[2]);
            add_value(values, col_b  , face_fac * 3./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_n  , face_fac * 3./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_nn , face_fac * 3./4. * dtinv * con_fac * w_ess[2]);

            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b  , face_fac * 3./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_n  , face_fac * 3./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_nn , face_fac * 3./4. * dtinv * con_fac * w_ess[2]);
            add_value(values, col_wb , face_fac * 1./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_wn , face_fac * 1./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_wnn, face_fac * 1./4. * dtinv * con_fac * w_ess[2]);

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_n  + 1, 0.0);
            add_value(values, col_nn + 1, 0.0);
            //
            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dy_deta_0;
            add_value(values, col_eb  + 2, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_en  + 2, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_enn + 2, face_fac * 1./4. * dtinv * w_ess[2]);
            add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_n   + 2, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_nn  + 2, face_fac * 3./4. * dtinv * w_ess[2]);

            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_n   + 2, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_nn  + 2, face_fac * 3./4. * dtinv * w_ess[2]);
            add_value(values, col_wb  + 2, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_wn  + 2, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_wnn + 2, face_fac * 1./4. * dtinv * w_ess[2]);
            //
            double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_e);
            double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_w);
            double dhdt_e = dtinv * (hp_e - hn_e);
            double dhdt_b = dtinv * (hp_b - hn_b);
            double dhdt_w = dtinv * (hp_w - hn_w);
            double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_e);
            double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_w);
            double drdt_e = dtinv * (rp_e - rn_e);
            double drdt_b = dtinv * (rp_b - rn_b);
            double drdt_w = dtinv * (rp_w - rn_w);
            double drdt_0 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_e);
            double drdt_1 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_w);
            rhs[row] = - ( 0.5 * dy_deta_0 * (drdt_0 + con_fac * dhdt_0) 
                         + 0.5 * dy_deta_1 * (drdt_1 + con_fac * dhdt_1) 
                         );

            double corr_term = 0.0;
            if (bc_vars[BC_SOUTH] == "zeta")
            {
                // face 0
                face_fac = 0.5 * dy_deta_0;
                add_value(values, col_eb , face_fac * 1./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_en , face_fac * 1./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_enn, face_fac * 1./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b  , face_fac * 3./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_n  , face_fac * 3./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nn , face_fac * 3./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dy_deta_1;
                add_value(values, col_b  , face_fac * 3./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_n  , face_fac * 3./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nn , face_fac * 3./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_wb , face_fac * 1./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_wn , face_fac * 1./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_wnn, face_fac * 1./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                double zb_e = w_ess[0] * zb[p_6] + w_ess[1] * zb[p_7] + w_ess[2] * zb[p_8];
                double zb_b = w_ess[0] * zb[p_3] + w_ess[1] * zb[p_4] + w_ess[2] * zb[p_5];
                double zb_w = w_ess[0] * zb[p_0] + w_ess[1] * zb[p_1] + w_ess[2] * zb[p_2];
                double zb_0 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_e);
                double zb_1 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_w);

                corr_term = ( 0.5 * dy_deta_0 * ( dhdt_0 + ( eps_bc_corr * ((bc[BC_SOUTH] - zb_0) - htheta_0)) )
                            + 0.5 * dy_deta_1 * ( dhdt_1 + ( eps_bc_corr * ((bc[BC_SOUTH] - zb_1) - htheta_1)) )
                            );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_SOUTH] == "q")
            {
                // face 0
                face_fac = 0.5 * dy_deta_0;
                add_value(values, col_eb  + 2, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_en  + 2, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_enn + 2, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b   + 2, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_n   + 2, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nn  + 2, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                face_fac = 0.5 * dy_deta_1;
                add_value(values, col_b   + 2, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_n   + 2, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nn  + 2, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_wb  + 2, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_wn  + 2, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_wnn + 2, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                double rtheta_0 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_e);
                double rtheta_1 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_w);

                corr_term = ( 0.5 * dy_deta_0 * (- drdt_0 + eps_bc_corr * (bc[BC_SOUTH] - rtheta_0)) +
                              0.5 * dy_deta_1 * (- drdt_1 + eps_bc_corr * (bc[BC_SOUTH] - rtheta_1))
                            );
                rhs[row] += corr_term;
            }
        }
        //------------------------------------------------------------------
        // natural boundary condition
        // -----------------------------------------------------------------
        double hp_e = w_nat[0] * hp[p_6] + w_nat[1] * hp[p_7] + w_nat[2] * hp[p_8];
        double hp_b = w_nat[0] * hp[p_3] + w_nat[1] * hp[p_4] + w_nat[2] * hp[p_5];
        double hp_w = w_nat[0] * hp[p_0] + w_nat[1] * hp[p_1] + w_nat[2] * hp[p_2];
        //
        double hn_e = w_nat[0] * hn[p_6] + w_nat[1] * hn[p_7] + w_nat[2] * hn[p_8];
        double hn_b = w_nat[0] * hn[p_3] + w_nat[1] * hn[p_4] + w_nat[2] * hn[p_5];
        double hn_w = w_nat[0] * hn[p_0] + w_nat[1] * hn[p_1] + w_nat[2] * hn[p_2];
        //
        double drdy_e = rtheta[p_7] - rtheta[p_6];
        double drdy_b = rtheta[p_4] - rtheta[p_3];
        double drdy_w = rtheta[p_1] - rtheta[p_0];
        //
        double rp_e = w_nat[0] * rp[p_6] + w_nat[1] * rp[p_7] + w_nat[2] * rp[p_8];
        double rp_b = w_nat[0] * rp[p_3] + w_nat[1] * rp[p_4] + w_nat[2] * rp[p_5];
        double rp_w = w_nat[0] * rp[p_0] + w_nat[1] * rp[p_1] + w_nat[2] * rp[p_2];
        //
        double rn_e = w_nat[0] * rn[p_6] + w_nat[1] * rn[p_7] + w_nat[2] * rn[p_8];
        double rn_b = w_nat[0] * rn[p_3] + w_nat[1] * rn[p_4] + w_nat[2] * rn[p_5];
        double rn_w = w_nat[0] * rn[p_0] + w_nat[1] * rn[p_1] + w_nat[2] * rn[p_2];
        //
        double htheta_e = w_nat[0] * htheta[p_6] + w_nat[1] * htheta[p_7] + w_nat[2] * htheta[p_8];  // east of boundary location
        double htheta_b = w_nat[0] * htheta[p_3] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_5];
        double htheta_w = w_nat[0] * htheta[p_0] + w_nat[1] * htheta[p_1] + w_nat[2] * htheta[p_2];  // west of boundary location
        //
        double qtheta_e = w_nat[0] * qtheta[p_6] + w_nat[1] * qtheta[p_7] + w_nat[2] * qtheta[p_8];  // north of boundary location
        double qtheta_b = w_nat[0] * qtheta[p_3] + w_nat[1] * qtheta[p_4] + w_nat[2] * qtheta[p_5];
        double qtheta_w = w_nat[0] * qtheta[p_0] + w_nat[1] * qtheta[p_1] + w_nat[2] * qtheta[p_2];  // south of boundary location
        //
        double rtheta_e = w_nat[0] * rtheta[p_6] + w_nat[1] * rtheta[p_7] + w_nat[2] * rtheta[p_8];  // east of boundary location
        double rtheta_b = w_nat[0] * rtheta[p_3] + w_nat[1] * rtheta[p_4] + w_nat[2] * rtheta[p_5];
        double rtheta_w = w_nat[0] * rtheta[p_0] + w_nat[1] * rtheta[p_1] + w_nat[2] * rtheta[p_2];  // west of boundary location
        //
        double dzetady_e = htheta[p_7] + zb[p_7] - htheta[p_6] - zb[p_6];
        double dzetady_b = htheta[p_4] + zb[p_4] - htheta[p_3] - zb[p_3];
        double dzetady_w = htheta[p_1] + zb[p_1] - htheta[p_0] - zb[p_0];
        // ---------------------------------------------------------------------
        if (do_convection) { con_fac = c_wave + rp_b / hp_b; }
        // ---------------------------------------------------------------------
        // second row
        // q-momentum (tangential equation, q == 0)
        //
        size_t col_wb  = q_eq + 0 * 3;
        size_t col_wn  = q_eq + 1 * 3;
        size_t col_wnn = q_eq + 2 * 3;
        size_t col_b   = q_eq + 3 * 3;
        size_t col_n   = q_eq + 4 * 3;
        size_t col_nn  = q_eq + 5 * 3;
        size_t col_eb  = q_eq + 6 * 3;
        size_t col_en  = q_eq + 7 * 3;
        size_t col_enn = q_eq + 8 * 3;
        //
        add_value(values, col_b, 0.0);
        add_value(values, col_n, 0.0);
        add_value(values, col_nn, 0.0);
        //
        add_value(values, col_b  + 1, 1.0);
        add_value(values, col_n  + 1, 0.0);
        add_value(values, col_nn + 1, 0.0);
        //
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_n  + 2, 0.0);
        add_value(values, col_nn + 2, 0.0);
        //
        rhs[row + 1] = 0.0;
        //
        // ---------------------------------------------------------------------
        // third row
        // momentum part dr/dt + gh d(zeta)/dy
        //
        col_wb  = r_eq + 0 * 3;
        col_wn  = r_eq + 1 * 3;
        col_wnn = r_eq + 2 * 3;
        col_b   = r_eq + 3 * 3;
        col_n   = r_eq + 4 * 3;
        col_nn  = r_eq + 5 * 3;
        col_eb  = r_eq + 6 * 3;
        col_en  = r_eq + 7 * 3;
        col_enn = r_eq + 8 * 3;
        //
        // Contribution Delta h
        // face 0
        double face_fac = 0.5;
        add_value(values, col_eb , face_fac * 1./4. * w_nat[0] * theta * g * dzetady_e - face_fac * 1./4. * theta * g * htheta_e);
        add_value(values, col_en , face_fac * 1./4. * w_nat[1] * theta * g * dzetady_e + face_fac * 1./4. * theta * g * htheta_e);
        add_value(values, col_enn, face_fac * 1./4. * w_nat[2] * theta * g * dzetady_e);
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetady_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_n  , face_fac * 3./4. * w_nat[1] * theta * g * dzetady_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_nn , face_fac * 3./4. * w_nat[2] * theta * g * dzetady_b);
        
        //face 1
        face_fac = 0.5;
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetady_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_n  , face_fac * 3./4. * w_nat[1] * theta * g * dzetady_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_nn , face_fac * 3./4. * w_nat[2] * theta * g * dzetady_b);
        add_value(values, col_wb , face_fac * 1./4. * w_nat[0] * theta * g * dzetady_w - face_fac * 1./4. * theta * g * htheta_w);
        add_value(values, col_wn , face_fac * 1./4. * w_nat[1] * theta * g * dzetady_w + face_fac * 1./4. * theta * g * htheta_w);
        add_value(values, col_wnn, face_fac * 1./4. * w_nat[2] * theta * g * dzetady_w);

        // Contribution Delta q
        face_fac = dx_dxi_0;
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_n  + 1, 0.0);
        add_value(values, col_nn + 1, 0.0);

        // Contribution Delta r
        // face 0
        face_fac = 0.5 * dx_dxi_0;
        add_value(values, col_eb  + 2, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_en  + 2, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_enn + 2, face_fac * 1./4. * dtinv * w_nat[2]);
        add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_n   + 2, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_nn  + 2, face_fac * 3./4. * dtinv * w_nat[2]);
        
        //face 1
        face_fac = 0.5 * dx_dxi_1;
        add_value(values, col_b   + 2, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_n   + 2, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_nn  + 2, face_fac * 3./4. * dtinv * w_nat[2]);
        add_value(values, col_wb  + 2, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_wn  + 2, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_wnn + 2, face_fac * 1./4. * dtinv * w_nat[2]);
        //
        double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_e);
        double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_w);
        double rtheta_0 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_e);
        double rtheta_1 = 0.25 * (3.0 * rtheta_b + 1.0 * rtheta_w);
        double dzetady_0 = 0.25 * (3.0 * dzetady_b + 1.0 * dzetady_e);
        double dzetady_1 = 0.25 * (3.0 * dzetady_b + 1.0 * dzetady_w);
        double drdt_e = dtinv * (rp_e - rn_e);
        double drdt_b = dtinv * (rp_b - rn_b);
        double drdt_w = dtinv * (rp_w - rn_w);
        double drdt_0 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_e);
        double drdt_1 = 0.25 * ( 3.0 * drdt_b + 1.0 * drdt_w);

        rhs[row + 2] = - ( 0.5 * ( dx_dxi_0 * drdt_0 + g * htheta_0 * dzetady_0 ) 
                         + 0.5 * ( dx_dxi_1 * drdt_1 + g * htheta_1 * dzetady_1 ) );
        //
        if (do_convection)
        {
            // South boundary convection
            
            double dh_deta_0 = 0.25 * (3.0 * (htheta[p_4] - htheta[p_3]) +  1.0 * (htheta[p_7] - htheta[p_6]));
            double dr_deta_0 = 0.25 * (3.0 * (rtheta[p_4] - rtheta[p_3]) +  1.0 * (rtheta[p_7] - rtheta[p_6]));
            double dh_deta_1 = 0.25 * (3.0 * (htheta[p_4] - htheta[p_3]) +  1.0 * (htheta[p_1] - htheta[p_0]));
            double dr_deta_1 = 0.25 * (3.0 * (rtheta[p_4] - rtheta[p_3]) +  1.0 * (rtheta[p_1] - rtheta[p_0]));

            double dh_dxi_0 = 0.5 * ((htheta[p_7] - htheta[p_4]) + (htheta[p_6] - htheta[p_3]));
            double dr_dxi_0 = 0.5 * ((rtheta[p_7] - rtheta[p_4]) + (rtheta[p_6] - rtheta[p_3]));
            double dh_dxi_1 = 0.5 * ((htheta[p_4] - htheta[p_1]) + (htheta[p_3] - htheta[p_0]));
            double dr_dxi_1 = 0.5 * ((rtheta[p_4] - rtheta[p_1]) + (rtheta[p_3] - rtheta[p_0]));

            double dx_deta_0 = dcdy_scvf_n(y[p_4], y[p_3], y[p_7], y[p_6]);
            double dx_deta_1 = dcdy_scvf_n(y[p_4], y[p_3], y[p_1], y[p_0]);

            double aa_0 = - 2. * rtheta_0 / (htheta_0 * htheta_0) * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) + 2. * (rtheta_0 * rtheta_0) / (htheta_0 * htheta_0 * htheta_0) * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0);
            double aa_1 = - 2. * rtheta_1 / (htheta_1 * htheta_1) * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) + 2. * (rtheta_1 * rtheta_1) / (htheta_1 * htheta_1 * htheta_1) * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1);
            
            double bb_0 = 2. / htheta_0 * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) - 2. * rtheta_0 / (htheta_0 * htheta_0) * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0);
            double bb_1 = 2. / htheta_1 * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) - 2. * rtheta_1 / (htheta_1 * htheta_1) * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1);
            
            double cc_0 = - rtheta_0 * rtheta_0 / (htheta_0 * htheta_0);
            double cc_1 = - rtheta_1 * rtheta_1 / (htheta_1 * htheta_1);
            
            double dd_0 = 2. * rtheta_0 / htheta_0;
            double dd_1 = 2. * rtheta_1 / htheta_1;
            
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * theta / dy_deta_0;
            add_value(values, col_eb , face_fac * 1./4. * (aa_0 * w_nat[0] - dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_en , face_fac * 1./4. * (aa_0 * w_nat[1] + dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_enn, face_fac * 1./4. * (aa_0 * w_nat[2]));
            add_value(values, col_b  , face_fac * 3./4. * (aa_0 * w_nat[0] - dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_n  , face_fac * 3./4. * (aa_0 * w_nat[1] + dx_dxi_0 * cc_0));  // not yet implemeted term with dx_deta
            add_value(values, col_nn , face_fac * 3./4. * (aa_0 * w_nat[2]));

            //face 1
            face_fac = 0.5 * theta / dy_deta_1;
            add_value(values, col_b  , face_fac * 3./4. * (aa_1 * w_nat[0] - dx_dxi_1 * cc_1));  // not yet implemeted term with dx_deta
            add_value(values, col_n  , face_fac * 3./4. * (aa_1 * w_nat[1] + dx_dxi_1 * cc_1));  // not yet implemeted term with dx_deta
            add_value(values, col_nn , face_fac * 3./4. * (aa_1 * w_nat[2]));
            add_value(values, col_wb , face_fac * 1./4. * (aa_1 * w_nat[0] - dx_dxi_1 * cc_1));  // not yet implemeted term with dx_deta
            add_value(values, col_wn , face_fac * 1./4. * (aa_1 * w_nat[1] + dx_dxi_1 * cc_1));  // not yet implemeted term with dx_deta
            add_value(values, col_wnn, face_fac * 1./4. * (aa_1 * w_nat[2]));
            
            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_n  + 1, 0.0);
            add_value(values, col_nn + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * theta / dy_deta_0;
            add_value(values, col_eb  + 2, face_fac * 1./4. * (bb_0 * w_nat[0] - dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_en  + 2, face_fac * 1./4. * (bb_0 * w_nat[1] + dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_enn + 2, face_fac * 1./4. * (bb_0 * w_nat[2]));
            add_value(values, col_b   + 2, face_fac * 3./4. * (bb_0 * w_nat[0] - dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_n   + 2, face_fac * 3./4. * (bb_0 * w_nat[1] + dx_dxi_0 * dd_0));  // not yet implemeted term with dx_deta
            add_value(values, col_nn  + 2, face_fac * 3./4. * (bb_0 * w_nat[2]));

            // face 1
            face_fac = 0.5 * theta / dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3./4. * (bb_1 * w_nat[0] - dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_n   + 2, face_fac * 3./4. * (bb_1 * w_nat[1] + dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_nn  + 2, face_fac * 3./4. * (bb_1 * w_nat[2]));
            add_value(values, col_wb  + 2, face_fac * 1./4. * (bb_1 * w_nat[0] - dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_wn  + 2, face_fac * 1./4. * (bb_1 * w_nat[1] + dx_dxi_1 * dd_1));  // not yet implemeted term with dx_deta
            add_value(values, col_wnn + 2, face_fac * 1./4. * (bb_1 * w_nat[2]));
                                                
           rhs[row + 2] += - (
                0.5 * dd_0 * (-dx_deta_0 * dr_dxi_0 + dx_dxi_0 * dr_deta_0) / dy_deta_0
              + 0.5 * cc_0 * (-dx_deta_0 * dh_dxi_0 + dx_dxi_0 * dh_deta_0) / dy_deta_0
              + 0.5 * dd_1 * (-dx_deta_1 * dr_dxi_1 + dx_dxi_1 * dr_deta_1) / dy_deta_1
              + 0.5 * cc_1 * (-dx_deta_1 * dh_dxi_1 + dx_dxi_1 * dh_deta_1) / dy_deta_1
                );
        }
        if (do_bed_shear_stress)
        {
            // South boundary bed shear stress

            double cf_w = cf;
            double cf_b = cf;
            double cf_e = cf;
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dy_deta_0;
            add_value(values, col_eb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_en , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_enn, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_n  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nn , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_n  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nn , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_wb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wn , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wnn, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_21(htheta_w, qtheta_w, rtheta_w, cf_w)) );    

            // Contribution Delta q
            add_value(values, col_b  + 1, 0.0);
            add_value(values, col_n  + 1, 0.0);
            add_value(values, col_nn + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dy_deta_0;
            add_value(values, col_eb  + 2, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_en  + 2, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_enn + 2, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_e, qtheta_e, rtheta_e, cf_e)) );    
            add_value(values, col_b   + 2, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_n   + 2, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nn  + 2, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dy_deta_1;
            add_value(values, col_b   + 2, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_n   + 2, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nn  + 2, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_wb  + 2, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wn  + 2, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );    
            add_value(values, col_wnn + 2, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_23(htheta_w, qtheta_w, rtheta_w, cf_w)) );

            // right hand side
            double qtheta_0 = 0.0;  // no tangential discharge q
            double qtheta_1 = 0.0;
            double abs_qtheta_0 = abs_vecq(qtheta_0, rtheta_0, 1.0);
            double abs_qtheta_1 = abs_vecq(qtheta_1, rtheta_1, 1.0);

            double cf_0 = 0.25 * (3.0 * cf_b + 1.0 * cf_e);
            double cf_1 = 0.25 * (3.0 * cf_b + 1.0 * cf_w);

            rhs[row + 2] += - (
                  0.5 * dy_deta_0 * cf_0 * rtheta_0 * abs_qtheta_0 / (htheta_0 * htheta_0)
                + 0.5 * dy_deta_1 * cf_1 * rtheta_1 * abs_qtheta_1 / (htheta_1 * htheta_1)

            );
        }
        do_viscosity = false;  // south boundary
        if (do_viscosity)
        {
            // South boundary viscosity
        }
        //
        // continuity part -c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        // face 0
        face_fac = 0.5 * dy_deta_0;
        add_value(values, col_eb , face_fac * 1./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_en , face_fac * 1./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_enn, face_fac * 1./4. * -con_fac * dtinv * w_nat[2]);
        add_value(values, col_b  , face_fac * 3./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_n  , face_fac * 3./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_nn , face_fac * 3./4. * -con_fac * dtinv * w_nat[2]);

        // face 1
        face_fac = 0.5 * dy_deta_1;
        add_value(values, col_b  , face_fac * 3./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_n  , face_fac * 3./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_nn , face_fac * 3./4. * -con_fac * dtinv * w_nat[2]);
        add_value(values, col_wb , face_fac * 1./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_wn , face_fac * 1./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_wnn, face_fac * 1./4. * -con_fac * dtinv * w_nat[2]);

        // Contribution Delta q
        face_fac = dx_dxi_0;
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_n  + 1, 0.0);
        add_value(values, col_nn + 1, 0.0);

        // Contribution Delta r
        // face 0
        face_fac = 0.5;
        add_value(values, col_eb  + 2, face_fac * 1./4. * -con_fac * -theta);
        add_value(values, col_en  + 2, face_fac * 1./4. * -con_fac *  theta);
        add_value(values, col_enn + 2, face_fac * 1./4. * 0.0);
        add_value(values, col_b   + 2, face_fac * 3./4. * -con_fac * -theta);
        add_value(values, col_n   + 2, face_fac * 3./4. * -con_fac *  theta);
        add_value(values, col_nn  + 2, face_fac * 3./4. * 0.0);

        // face 1
        face_fac = 0.5;
        add_value(values, col_b   + 2, face_fac * 3./4. * -con_fac * -theta);
        add_value(values, col_n   + 2, face_fac * 3./4. * -con_fac *  theta);
        add_value(values, col_nn  + 2, face_fac * 3./4. * 0.0);
        add_value(values, col_wb  + 2, face_fac * 1./4. * -con_fac * -theta);
        add_value(values, col_wn  + 2, face_fac * 1./4. * -con_fac *  theta);
        add_value(values, col_wnn + 2, face_fac * 1./4. * 0.0);

        //
        double dhdt_e = dtinv * (hp_e - hn_e);
        double dhdt_b = dtinv * (hp_b - hn_b);
        double dhdt_w = dtinv * (hp_w - hn_w);
        double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_e);
        double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_w);
        double drdy_0 = 0.25 * (3.0 * drdy_b + 1.0 * drdy_e);
        double drdy_1 = 0.25 * (3.0 * drdy_b + 1.0 * drdy_w);

        rhs[row + 2] += - ( - 0.5 * con_fac * (dy_deta_0 * dhdt_0 + drdy_0)
                            - 0.5 * con_fac * (dy_deta_1 * dhdt_1 + drdy_1) );
    }
    return 0;
}
int boundary_west(double* values, size_t row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    size_t nx, size_t ny,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, double cf,
    std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_WEST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    std::fill_n(values + c_eq, 3 * 27, 0.0);  // set all coefficients for one row of Delta c-, Delta q- and Delta r-equation to zero

    size_t p_1 = c_eq/(3*27);  // node number of boundary point, ie west point of molucule
    size_t p_0 = p_1 - 1;
    size_t p_2 = p_1 + 1;
    size_t p_3 = p_1 + ny - 1;
    size_t p_4 = p_1 + ny;
    size_t p_5 = p_1 + ny + 1;
    size_t p_6 = p_1 + 2 * ny - 1;
    size_t p_7 = p_1 + 2 * ny;
    size_t p_8 = p_1 + 2 * ny + 1;

    double htheta_b = w_nat[0] * htheta[p_1] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_7];
    double zb_b = w_nat[0] * zb[p_1] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_7];

    double h_given = bc[BC_WEST] - zb_b;
    double h_infty = -zb_b;
    double c_wave = std::sqrt(g * htheta_b);
    double con_fac = c_wave;

    // wwest
    double dx_dxi_s = x[p_3] - x[p_0];
    double dx_dxi_b = x[p_4] - x[p_1];
    double dx_dxi_n = x[p_5] - x[p_2];
    double dy_deta_0 = 0.5 * (y[p_4] - y[p_3]) + 0.5 * (y[p_1] - y[p_0]);
    double dy_deta_1 = 0.5 * (y[p_5] - y[p_4]) + 0.5 * (y[p_2] - y[p_1]);
    double dx_dxi_0 = 0.25 * (3. * (x[p_4] - x[p_1]) + 1. * (x[p_3] - x[p_0]));
    double dx_dxi_1 = 0.25 * (3. * (x[p_4] - x[p_1]) + 1. * (x[p_5] - x[p_2]));

    if (bc_type[BC_WEST] == "no_slip" || bc_type[BC_WEST] == "free_slip")
    {
        if (bc_type[BC_WEST] == "no_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 1 * 3; // point of boundary, ie west point of molecule
            size_t col_e  = c_eq + 4 * 3;
            size_t col_ee = c_eq + 7 * 3;
            values[col_b ] = -1.0 * theta;
            values[col_e ] =  1.0 * theta;
            values[col_ee] = 0.0;
            rhs[row] = -(htheta[p_4] + zb[p_4] - htheta[p_1] - zb[p_1]);

            // Contribution Delta q
            col_b  = q_eq + 1 * 3;
            col_e  = q_eq + 4 * 3;
            col_ee = q_eq + 7 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_e  + 1] = 0.5 * theta;
            values[col_ee + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 1 * 3;
            col_e  = r_eq + 4 * 3;
            col_ee = r_eq + 7 * 3;
            values[col_b  + 2]  = 0.5 * theta;
            values[col_e  + 2]  = 0.5 * theta;
            values[col_ee + 2]  = 0.0;
            rhs[row + 2] = 0.0;
        }
        else if (bc_type[BC_WEST] == "free_slip")
        {
            // Contribution Delta h
            size_t col_b  = c_eq + 1 * 3; // point of boundary, ie west point of molecule
            size_t col_e  = c_eq + 4 * 3;
            size_t col_ee = c_eq + 7 * 3;
            values[col_b ] = -1.0 * theta;
            values[col_e ] =  1.0 * theta;
            values[col_ee] = 0.0;
            rhs[row] = -(htheta[p_4] + zb[p_4] - htheta[p_1] - zb[p_1]);

            // Contribution Delta q
            col_b  = q_eq + 1 * 3;
            col_e  = q_eq + 4 * 3;
            col_ee = q_eq + 7 * 3;
            values[col_b  + 1] = 0.5 * theta;
            values[col_e  + 1] = 0.5 * theta;
            values[col_ee + 1] = 0.0;
            rhs[row + 1] = 0.0;

            // Contribution Delta r
            col_b  = r_eq + 1 * 3;
            col_e  = r_eq + 4 * 3;
            col_ee = r_eq + 7 * 3;
            values[col_b  + 2]  = -1.0 * theta;
            values[col_e  + 2]  =  1.0 * theta;
            values[col_ee + 2]  = 0.0;
            rhs[row + 2] = -(rtheta[p_4] - rtheta[p_1]);
        }
    }
    else
    {
        if (bc_type[BC_WEST] == "mooiman")
        {
            // essential boundary condition
            // first row
            size_t col_sb  = c_eq + 0 * 3;
            size_t col_se  = c_eq + 3 * 3;
            size_t col_see = c_eq + 6 * 3;
            size_t col_b  = c_eq + 1 * 3;
            size_t col_e  = c_eq + 4 * 3;
            size_t col_ee = c_eq + 7 * 3;
            size_t col_nb  = c_eq + 2 * 3;
            size_t col_ne  = c_eq + 5 * 3;
            size_t col_nee = c_eq + 8 * 3;
            //
            // face 0
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb , face_fac * 1./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_se , face_fac * 1./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_see, face_fac * 1./4. * theta * w_ess[2] * c_wave);
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_e  , face_fac * 3./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_ee , face_fac * 3./4. * theta * w_ess[2] * c_wave);
            
            //face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_e  , face_fac * 3./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_ee , face_fac * 3./4. * theta * w_ess[2] * c_wave);
            add_value(values, col_nb , face_fac * 1./4. * theta * w_ess[0] * c_wave);
            add_value(values, col_ne , face_fac * 1./4. * theta * w_ess[1] * c_wave);
            add_value(values, col_nee, face_fac * 1./4. * theta * w_ess[2] * c_wave);
            
            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb  + 1, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_se  + 1, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_see + 1, face_fac * 1./4. * theta * w_ess[2]);
            add_value(values, col_b   + 1, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_e   + 1, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ee  + 1, face_fac * 3./4. * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * theta * w_ess[0]);
            add_value(values, col_e   + 1, face_fac * 3./4. * theta * w_ess[1]);
            add_value(values, col_ee  + 1, face_fac * 3./4. * theta * w_ess[2]);
            add_value(values, col_nb  + 1, face_fac * 1./4. * theta * w_ess[0]);
            add_value(values, col_ne  + 1, face_fac * 1./4. * theta * w_ess[1]);
            add_value(values, col_nee + 1, face_fac * 1./4. * theta * w_ess[2]);
            
            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_e  + 2, 0.0);
            add_value(values, col_ee + 2, 0.0);
            //
            // Right hand side
            face_fac = 1.0;
            double htheta_b   = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_2]);
            double htheta_ip1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_5]);
            double htheta_ip2 = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_8]);
            double ht_0 = face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2);
            htheta_b   = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_0]);
            htheta_ip1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_3]);
            htheta_ip2 = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_7]);            
            double ht_1 = face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2);

            double qtheta_b   = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_2]);
            double qtheta_ip1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_5]);
            double qtheta_ip2 = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_8]);
            double qt_0 = face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2);
            qtheta_b   = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_0]);
            qtheta_ip1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_3]);
            qtheta_ip2 = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_6]);
            double qt_1 = face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2);

            rhs[row] = -( 0.5 * dx_dxi_0 * ( qt_0 + c_wave * (ht_0 - h_given) ) +
                          0.5 * dx_dxi_1 * ( qt_1 + c_wave * (ht_1 - h_given) ) )
                ;
            if (bc_vars[BC_WEST] == "zeta")
            {
                rhs[row] += 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * c_wave * bc[BC_WEST];
            }
            if (bc_vars[BC_WEST] == "q")
            {
                rhs[row] += 0.5 * (dx_dxi_0 + dx_dxi_1) * 2. * bc[BC_WEST];
            }
        }
        else if (bc_type[BC_WEST] == "borsboom")
        {
            //----------------------------------------------------------------------
            // Essential boundary condition
            //----------------------------------------------------------------------
            double hp_n = w_ess[0] * hp[p_2] + w_ess[1] * hp[p_5] + w_ess[2] * hp[p_8];
            double hp_b = w_ess[0] * hp[p_1] + w_ess[1] * hp[p_4] + w_ess[2] * hp[p_7];
            double hp_s = w_ess[0] * hp[p_0] + w_ess[1] * hp[p_3] + w_ess[2] * hp[p_6];
        
            double hn_n = w_ess[0] * hn[p_2] + w_ess[1] * hn[p_5] + w_ess[2] * hn[p_8];
            double hn_b = w_ess[0] * hn[p_1] + w_ess[1] * hn[p_4] + w_ess[2] * hn[p_7];
            double hn_s = w_ess[0] * hn[p_0] + w_ess[1] * hn[p_3] + w_ess[2] * hn[p_6];

            double htheta_n = w_ess[0] * htheta[p_2] + w_ess[1] * htheta[p_5] + w_ess[2] * htheta[p_8];  // north of boundary location
            double htheta_b = w_ess[0] * htheta[p_1] + w_ess[1] * htheta[p_4] + w_ess[2] * htheta[p_7];
            double htheta_s = w_ess[0] * htheta[p_0] + w_ess[1] * htheta[p_3] + w_ess[2] * htheta[p_6];  // south of boundary location

            double qp_n = w_ess[0] * qp[p_2] + w_ess[1] * qp[p_5] + w_ess[2] * qp[p_8];
            double qp_b = w_ess[0] * qp[p_1] + w_ess[1] * qp[p_4] + w_ess[2] * qp[p_7];
            double qp_s = w_ess[0] * qp[p_0] + w_ess[1] * qp[p_3] + w_ess[2] * qp[p_6];
        
            double qn_n = w_ess[0] * qn[p_2] + w_ess[1] * qn[p_5] + w_ess[2] * qn[p_8];
            double qn_b = w_ess[0] * qn[p_1] + w_ess[1] * qn[p_4] + w_ess[2] * qn[p_7];
            double qn_s = w_ess[0] * qn[p_0] + w_ess[1] * qn[p_3] + w_ess[2] * qn[p_6];
        
            double qtheta_n = w_ess[0] * qtheta[p_2] + w_ess[1] * qtheta[p_5] + w_ess[2] * qtheta[p_8];  // north of boundary location
            double qtheta_b = w_ess[0] * qtheta[p_1] + w_ess[1] * qtheta[p_4] + w_ess[2] * qtheta[p_7];
            double qtheta_s = w_ess[0] * qtheta[p_0] + w_ess[1] * qtheta[p_3] + w_ess[2] * qtheta[p_6];  // south of boundary location
            //------------------------------------------------------------------
            // first row
            //
            size_t col_sb  = c_eq + 0 * 3;
            size_t col_se  = c_eq + 3 * 3;
            size_t col_see = c_eq + 6 * 3;
            size_t col_b   = c_eq + 1 * 3;
            size_t col_e   = c_eq + 4 * 3;
            size_t col_ee  = c_eq + 7 * 3;
            size_t col_nb  = c_eq + 2 * 3;
            size_t col_ne  = c_eq + 5 * 3;
            size_t col_nee = c_eq + 8 * 3;
            //
            if (do_convection) { con_fac = c_wave - qp_b / hp_b; }
            //
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb , face_fac * 1./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_se , face_fac * 1./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_see, face_fac * 1./4. * dtinv * con_fac * w_ess[2]);
            add_value(values, col_b  , face_fac * 3./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_e  , face_fac * 3./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_ee , face_fac * 3./4. * dtinv * con_fac * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_e  , face_fac * 3./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_ee , face_fac * 3./4. * dtinv * con_fac * w_ess[2]);
            add_value(values, col_nb , face_fac * 1./4. * dtinv * con_fac * w_ess[0]);
            add_value(values, col_ne , face_fac * 1./4. * dtinv * con_fac * w_ess[1]);
            add_value(values, col_nee, face_fac * 1./4. * dtinv * con_fac * w_ess[2]);

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb  + 1, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_se  + 1, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_see + 1, face_fac * 1./4. * dtinv * w_ess[2]);
            add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_e   + 1, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ee  + 1, face_fac * 3./4. * dtinv * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_ess[0]);
            add_value(values, col_e   + 1, face_fac * 3./4. * dtinv * w_ess[1]);
            add_value(values, col_ee  + 1, face_fac * 3./4. * dtinv * w_ess[2]);
            add_value(values, col_nb  + 1, face_fac * 1./4. * dtinv * w_ess[0]);
            add_value(values, col_ne  + 1, face_fac * 1./4. * dtinv * w_ess[1]);
            add_value(values, col_nee + 1, face_fac * 1./4. * dtinv * w_ess[2]);

            // Contribution for Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_e  + 2, 0.0);
            add_value(values, col_ee + 2, 0.0);
            //
            double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_s);
            double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_n);
            double dhdt_s = dtinv * (hp_s - hn_s);
            double dhdt_b = dtinv * (hp_b - hn_b);
            double dhdt_n = dtinv * (hp_n - hn_n);
            double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_s);
            double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_n);
            double dqdt_s = dtinv * (qp_s - qn_s);
            double dqdt_b = dtinv * (qp_b - qn_b);
            double dqdt_n = dtinv * (qp_n - qn_n);
            double dqdt_0 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_s);
            double dqdt_1 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_n);

            rhs[row] = - ( 0.5 * dx_dxi_0 * (dqdt_0 + con_fac * dhdt_0) 
                         + 0.5 * dx_dxi_1 * (dqdt_1 + con_fac * dhdt_1) 
                         );

            double corr_term = 0.0;
            if (bc_vars[BC_WEST] == "zeta")
            {
                // face 0
                face_fac = 0.5 * dx_dxi_0;
                add_value(values, col_sb , face_fac * 1./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_se , face_fac * 1./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_see, face_fac * 1./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b  , face_fac * 3./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_e  , face_fac * 3./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ee , face_fac * 3./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dx_dxi_1;
                add_value(values, col_b  , face_fac * 3./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_e  , face_fac * 3./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ee , face_fac * 3./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_nb , face_fac * 1./4. * (-dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_ne , face_fac * 1./4. * (-dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nee, face_fac * 1./4. * (-dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]));

                double zb_s = w_ess[0] * zb[p_0] + w_ess[1] * zb[p_3] + w_ess[2] * zb[p_6];
                double zb_b = w_ess[0] * zb[p_1] + w_ess[1] * zb[p_4] + w_ess[2] * zb[p_7];
                double zb_n = w_ess[0] * zb[p_2] + w_ess[1] * zb[p_5] + w_ess[2] * zb[p_8];
                double zb_0 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_s);
                double zb_1 = 0.25 * ( 3.0 * zb_b + 1.0 * zb_n);

                corr_term = ( 0.5 * dx_dxi_0 * ( dhdt_0 + ( eps_bc_corr * ((bc[BC_WEST] - zb_0) - htheta_0) ))
                            + 0.5 * dx_dxi_1 * ( dhdt_1 + ( eps_bc_corr * ((bc[BC_WEST] - zb_1) - htheta_1) )) 
                            );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_WEST] == "q")
            {
                // face 0
                face_fac = 0.5 * dx_dxi_0;
                add_value(values, col_sb  + 1, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_se  + 1, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_see + 1, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_b   + 1, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_e   + 1, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ee  + 1, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                // face 1
                face_fac = 0.5 * dx_dxi_1;
                add_value(values, col_b   + 1, face_fac * 3./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_e   + 1, face_fac * 3./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_ee  + 1, face_fac * 3./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));
                add_value(values, col_nb  + 1, face_fac * 1./4. * (dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]));
                add_value(values, col_ne  + 1, face_fac * 1./4. * (dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]));
                add_value(values, col_nee + 1, face_fac * 1./4. * (dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]));

                double qtheta_0 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_s);
                double qtheta_1 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_n);

                corr_term = ( 0.5 * dx_dxi_0 * (- dqdt_0 + eps_bc_corr * (bc[BC_WEST] - qtheta_0)) +
                              0.5 * dx_dxi_1 * (- dqdt_1 + eps_bc_corr * (bc[BC_WEST] - qtheta_1))
                            );
                rhs[row] += corr_term;
            }
        }
        //----------------------------------------------------------------------
        // natural boundary condition
        //----------------------------------------------------------------------
        double hp_s = w_nat[0] * hp[p_0] + w_nat[1] * hp[p_3] + w_nat[2] * hp[p_6];
        double hp_b = w_nat[0] * hp[p_1] + w_nat[1] * hp[p_4] + w_nat[2] * hp[p_7];
        double hp_n = w_nat[0] * hp[p_2] + w_nat[1] * hp[p_5] + w_nat[2] * hp[p_8];
        
        double hn_s = w_nat[0] * hn[p_0] + w_nat[1] * hn[p_3] + w_nat[2] * hn[p_6];
        double hn_b = w_nat[0] * hn[p_1] + w_nat[1] * hn[p_4] + w_nat[2] * hn[p_7];
        double hn_n = w_nat[0] * hn[p_2] + w_nat[1] * hn[p_5] + w_nat[2] * hn[p_8];

        double qp_s = w_nat[0] * qp[p_0] + w_nat[1] * qp[p_3] + w_nat[2] * qp[p_6];
        double qp_b = w_nat[0] * qp[p_1] + w_nat[1] * qp[p_4] + w_nat[2] * qp[p_7];
        double qp_n = w_nat[0] * qp[p_2] + w_nat[1] * qp[p_5] + w_nat[2] * qp[p_8];
        
        double qn_s = w_nat[0] * qn[p_0] + w_nat[1] * qn[p_3] + w_nat[2] * qn[p_6];
        double qn_b = w_nat[0] * qn[p_1] + w_nat[1] * qn[p_4] + w_nat[2] * qn[p_7];
        double qn_n = w_nat[0] * qn[p_2] + w_nat[1] * qn[p_5] + w_nat[2] * qn[p_8];
        
        double htheta_s = w_nat[0] * htheta[p_0] + w_nat[1] * htheta[p_3] + w_nat[2] * htheta[p_6];  // south of boundary location
        double htheta_b = w_nat[0] * htheta[p_1] + w_nat[1] * htheta[p_4] + w_nat[2] * htheta[p_7];
        double htheta_n = w_nat[0] * htheta[p_2] + w_nat[1] * htheta[p_5] + w_nat[2] * htheta[p_8];  // north of boundary location

        double qtheta_s = w_nat[0] * qtheta[p_0] + w_nat[1] * qtheta[p_3] + w_nat[2] * qtheta[p_6];  // south of boundary location
        double qtheta_b = w_nat[0] * qtheta[p_1] + w_nat[1] * qtheta[p_4] + w_nat[2] * qtheta[p_7];
        double qtheta_n = w_nat[0] * qtheta[p_2] + w_nat[1] * qtheta[p_5] + w_nat[2] * qtheta[p_8];  // north of boundary location

        double rtheta_s = w_nat[0] * rtheta[p_0] + w_nat[1] * rtheta[p_3] + w_nat[2] * rtheta[p_6];  // south of boundary location
        double rtheta_b = w_nat[0] * rtheta[p_1] + w_nat[1] * rtheta[p_4] + w_nat[2] * rtheta[p_7];
        double rtheta_n = w_nat[0] * rtheta[p_2] + w_nat[1] * rtheta[p_5] + w_nat[2] * rtheta[p_5];  // north of boundary location

        double dqdx_s = qtheta[p_3] - qtheta[p_0];
        double dqdx_b = qtheta[p_4] - qtheta[p_1];
        double dqdx_n = qtheta[p_5] - qtheta[p_2];

        double dzetadx_s = htheta[p_3] + zb[p_3] - htheta[p_0] - zb[p_0];
        double dzetadx_b = htheta[p_4] + zb[p_4] - htheta[p_1] - zb[p_1];
        double dzetadx_n = htheta[p_5] + zb[p_5] - htheta[p_2] - zb[p_2];
        // ---------------------------------------------------------------------
        if (do_convection) { con_fac = c_wave + qp_b / hp_b; }
        // ---------------------------------------------------------------------
        // second row
        // momentum part dq/dt + gh d(zeta)/dx 
        //
        size_t col_sb  = q_eq + 0 * 3;
        size_t col_se  = q_eq + 3 * 3;
        size_t col_see = q_eq + 6 * 3;
        size_t col_b   = q_eq + 1 * 3;
        size_t col_e   = q_eq + 4 * 3;
        size_t col_ee  = q_eq + 7 * 3;
        size_t col_nb  = q_eq + 2 * 3;
        size_t col_ne  = q_eq + 5 * 3;
        size_t col_nee = q_eq + 8 * 3;
        //
        // Contribution Delta h
        // face 0
        double face_fac = 0.5;
        add_value(values, col_sb , face_fac * 1./4. * w_nat[0] * theta * g * dzetadx_s - face_fac * 1./4. * theta * g * htheta_s);
        add_value(values, col_se , face_fac * 1./4. * w_nat[1] * theta * g * dzetadx_s + face_fac * 1./4. * theta * g * htheta_s);
        add_value(values, col_see, face_fac * 1./4. * w_nat[2] * theta * g * dzetadx_s);
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetadx_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_e  , face_fac * 3./4. * w_nat[1] * theta * g * dzetadx_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ee , face_fac * 3./4. * w_nat[2] * theta * g * dzetadx_b);

        //face 1
        face_fac = 0.5;
        add_value(values, col_b  , face_fac * 3./4. * w_nat[0] * theta * g * dzetadx_b - face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_e  , face_fac * 3./4. * w_nat[1] * theta * g * dzetadx_b + face_fac * 3./4. * theta * g * htheta_b);
        add_value(values, col_ee , face_fac * 3./4. * w_nat[2] * theta * g * dzetadx_b);
        add_value(values, col_nb , face_fac * 1./4. * w_nat[0] * theta * g * dzetadx_n - face_fac * 1./4. * theta * g * htheta_n);
        add_value(values, col_ne , face_fac * 1./4. * w_nat[1] * theta * g * dzetadx_n + face_fac * 1./4. * theta * g * htheta_n);
        add_value(values, col_nee, face_fac * 1./4. * w_nat[2] * theta * g * dzetadx_n);

        // Contribution Delta q
        // face 0
        face_fac = 0.5 * dy_deta_0;
        add_value(values, col_sb  + 1, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_se  + 1, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_see + 1, face_fac * 1./4. * dtinv * w_nat[2]);
        add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_e   + 1, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ee  + 1, face_fac * 3./4. * dtinv * w_nat[2]);

        //face 1
        face_fac = 0.5 * dy_deta_1;
        add_value(values, col_b   + 1, face_fac * 3./4. * dtinv * w_nat[0]);
        add_value(values, col_e   + 1, face_fac * 3./4. * dtinv * w_nat[1]);
        add_value(values, col_ee  + 1, face_fac * 3./4. * dtinv * w_nat[2]);
        add_value(values, col_nb  + 1, face_fac * 1./4. * dtinv * w_nat[0]);
        add_value(values, col_ne  + 1, face_fac * 1./4. * dtinv * w_nat[1]);
        add_value(values, col_nee + 1, face_fac * 1./4. * dtinv * w_nat[2]); 

        // Contribution Delta r
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_e  + 2, 0.0);
        add_value(values, col_ee + 2, 0.0);

        //
        double htheta_0 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_s);
        double htheta_1 = 0.25 * (3.0 * htheta_b + 1.0 * htheta_n);
        double qtheta_0 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_s);
        double qtheta_1 = 0.25 * (3.0 * qtheta_b + 1.0 * qtheta_n);
        double dzetadx_0 = 0.25 * (3.0 * dzetadx_b + 1.0 * dzetadx_s);
        double dzetadx_1 = 0.25 * (3.0 * dzetadx_b + 1.0 * dzetadx_n);
        double dqdt_s = dtinv * (qp_s - qn_s);
        double dqdt_b = dtinv * (qp_b - qn_b);
        double dqdt_n = dtinv * (qp_n - qn_n);
        double dqdt_0 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_s);
        double dqdt_1 = 0.25 * ( 3.0 * dqdt_b + 1.0 * dqdt_n);

        rhs[row + 1] = - ( 0.5 * ( dy_deta_0 * dqdt_0 + g * htheta_0 * dzetadx_0 ) 
                         + 0.5 * ( dy_deta_1 * dqdt_1 + g * htheta_1 * dzetadx_1 ) );
        //
        if (do_convection)
        {
            // West boundary convection
            
            double dh_dxi_0 = 0.25 * (3.0 * (htheta[p_4] - htheta[p_1]) +  1.0 * (htheta[p_3] - htheta[p_0]));
            double dq_dxi_0 = 0.25 * (3.0 * (qtheta[p_4] - qtheta[p_1]) +  1.0 * (qtheta[p_3] - qtheta[p_0]));
            double dh_dxi_1 = 0.25 * (3.0 * (htheta[p_4] - htheta[p_1]) +  1.0 * (htheta[p_5] - htheta[p_2]));
            double dq_dxi_1 = 0.25 * (3.0 * (qtheta[p_4] - qtheta[p_1]) +  1.0 * (qtheta[p_5] - qtheta[p_2]));

            double dh_deta_0 = 0.5 * ((htheta[p_4] - htheta[p_3]) + (htheta[p_1] - htheta[p_0]));
            double dq_deta_0 = 0.5 * ((qtheta[p_4] - qtheta[p_3]) + (qtheta[p_1] - qtheta[p_0]));
            double dh_deta_1 = 0.5 * ((htheta[p_5] - htheta[p_4]) + (htheta[p_2] - htheta[p_1]));
            double dq_deta_1 = 0.5 * ((qtheta[p_5] - qtheta[p_4]) + (qtheta[p_2] - qtheta[p_1]));

            double dy_dxi_0 = dcdx_scvf_n(y[p_4], y[p_1], y[p_3], y[p_0]);
            double dy_dxi_1 = dcdx_scvf_n(y[p_4], y[p_1], y[p_5], y[p_2]);

            double aa_0 = - 2. * qtheta_0 / (htheta_0 * htheta_0) * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) + 2. * (qtheta_0 * qtheta_0) / (htheta_0 * htheta_0 * htheta_0) * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0);
            double aa_1 = - 2. * qtheta_1 / (htheta_1 * htheta_1) * (dy_deta_1 * dq_dxi_1 - dy_dxi_1 * dq_deta_1) + 2. * (qtheta_1 * qtheta_1) / (htheta_1 * htheta_1 * htheta_1) * (dy_deta_1 * dh_dxi_1 - dy_dxi_1 * dh_deta_1);

            double bb_0 = 2. / htheta_0 * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) - 2. * qtheta_0 / (htheta_0 * htheta_0) * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0);
            double bb_1 = 2. / htheta_1 * (dy_deta_1 * dq_dxi_1 - dy_dxi_1 * dq_deta_1) - 2. * qtheta_1 / (htheta_1 * htheta_1) * (dy_deta_1 * dh_dxi_1 - dy_dxi_1 * dh_deta_1);

            double cc_0 = - qtheta_0 * qtheta_0 / (htheta_0 * htheta_0);
            double cc_1 = - qtheta_1 * qtheta_1 / (htheta_1 * htheta_1);

            double dd_0 = 2. * qtheta_0 / htheta_0;
            double dd_1 = 2. * qtheta_1 / htheta_1;

            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * theta / dx_dxi_0;
            add_value(values, col_sb , face_fac * 1./4. * (aa_0 * w_nat[0] - dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_se , face_fac * 1./4. * (aa_0 * w_nat[1] + dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_see, face_fac * 1./4. * (aa_0 * w_nat[2]));
            add_value(values, col_b  , face_fac * 3./4. * (aa_0 * w_nat[0] - dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_e  , face_fac * 3./4. * (aa_0 * w_nat[1] + dy_deta_0 * cc_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_ee , face_fac * 3./4. * (aa_0 * w_nat[2]));

            //face 1
            face_fac = 0.5 * theta / dx_dxi_1;
            add_value(values, col_b  , face_fac * 3./4. * (aa_1 * w_nat[0] - dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_e  , face_fac * 3./4. * (aa_1 * w_nat[1] + dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ee , face_fac * 3./4. * (aa_1 * w_nat[2]));
            add_value(values, col_nb , face_fac * 1./4. * (aa_1 * w_nat[0] - dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ne , face_fac * 1./4. * (aa_1 * w_nat[1] + dy_deta_1 * cc_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_nee, face_fac * 1./4. * (aa_1 * w_nat[2]));

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * theta / dx_dxi_0;
            add_value(values, col_sb  + 1, face_fac * 1./4. * (bb_0 * w_nat[0] - dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_se  + 1, face_fac * 1./4. * (bb_0 * w_nat[1] + dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_see + 1, face_fac * 1./4. * (bb_0 * w_nat[2]));
            add_value(values, col_b   + 1, face_fac * 3./4. * (bb_0 * w_nat[0] - dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_e   + 1, face_fac * 3./4. * (bb_0 * w_nat[1] + dy_deta_0 * dd_0));  // not yet implemeted term with dy_dxi
            add_value(values, col_ee  + 1, face_fac * 3./4. * (bb_0 * w_nat[2]));

            // face 1
            face_fac = 0.5 * theta / dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3./4. * (bb_1 * w_nat[0] - dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_e   + 1, face_fac * 3./4. * (bb_1 * w_nat[1] + dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ee  + 1, face_fac * 3./4. * (bb_1 * w_nat[2]));
            add_value(values, col_nb  + 1, face_fac * 1./4. * (bb_1 * w_nat[0] - dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_ne  + 1, face_fac * 1./4. * (bb_1 * w_nat[1] + dy_deta_1 * dd_1));  // not yet implemeted term with dy_dxi
            add_value(values, col_nee + 1, face_fac * 1./4. * (bb_1 * w_nat[2]));

            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_e  + 2, 0.0);
            add_value(values, col_ee + 2, 0.0);

            rhs[row + 1] += - (
                0.5 * dd_0 * (dy_deta_0 * dq_dxi_0 - dy_dxi_0 * dq_deta_0) / dx_dxi_0
              + 0.5 * cc_0 * (dy_deta_0 * dh_dxi_0 - dy_dxi_0 * dh_deta_0) / dx_dxi_0
              + 0.5 * dd_1 * (dy_deta_1 * dq_dxi_1 - dy_dxi_0 * dq_deta_1) / dx_dxi_1
              + 0.5 * cc_1 * (dy_deta_1 * dh_dxi_1 - dy_dxi_0 * dh_deta_1) / dx_dxi_1
                );
        }
        if (do_bed_shear_stress)
        {
            // West boundary bed shear stress 

            double cf_s = cf;
            double cf_b = cf;
            double cf_n = cf;
            // Contribution Delta h
            // face 0
            double face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_se , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_see, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_e  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ee , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b  , face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_e  , face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ee , face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nb , face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_ne , face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nee, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_11(htheta_n, qtheta_n, rtheta_n, cf_n)) );    

            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dx_dxi_0;
            add_value(values, col_sb  + 1, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_se  + 1, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_see + 1, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_s, qtheta_s, rtheta_s, cf_s)) );    
            add_value(values, col_b   + 1, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_e   + 1, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ee  + 1, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
                
            // face 1
            face_fac = 0.5 * dx_dxi_1;
            add_value(values, col_b   + 1, face_fac * 3.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_e   + 1, face_fac * 3.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_ee  + 1, face_fac * 3.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_b, qtheta_b, rtheta_b, cf_b)) );    
            add_value(values, col_nb  + 1, face_fac * 1.0/4.0 * (theta * w_nat[0] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_ne  + 1, face_fac * 1.0/4.0 * (theta * w_nat[1] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );    
            add_value(values, col_nee + 1, face_fac * 1.0/4.0 * (theta * w_nat[2] * bed_shear_stress_J_12(htheta_n, qtheta_n, rtheta_n, cf_n)) );

            // Contribution Delta r
            add_value(values, col_b  + 2, 0.0);
            add_value(values, col_e  + 2, 0.0);
            add_value(values, col_ee + 2, 0.0);

            // right hand side
            double rtheta_0 = 0.0;  // no tangential discharge r
            double rtheta_1 = 0.0;
            double abs_qtheta_0 = abs_vecq(qtheta_0, rtheta_0, 1.0);
            double abs_qtheta_1 = abs_vecq(qtheta_1, rtheta_1, 1.0);

            double cf_0 = 0.25 * (3.0 * cf_b + 1.0 * cf_s);
            double cf_1 = 0.25 * (3.0 * cf_b + 1.0 * cf_n);

            rhs[row + 1] += - (
                  0.5 * dx_dxi_0 * cf_0 * qtheta_0 * abs_qtheta_0 / (htheta_0 * htheta_0)
                + 0.5 * dx_dxi_1 * cf_1 * qtheta_1 * abs_qtheta_1 / (htheta_1 * htheta_1)
                );
        }
        do_viscosity = false;  // west boundary
        if (do_viscosity)
        {
            // West boundary viscosity
        }
        //
        // continuity part -c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        // face 0
        face_fac = 0.5 * dx_dxi_0;
        add_value(values, col_sb , face_fac * 1./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_se , face_fac * 1./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_see, face_fac * 1./4. * -con_fac * dtinv * w_nat[2]);
        add_value(values, col_b  , face_fac * 3./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_e  , face_fac * 3./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_ee , face_fac * 3./4. * -con_fac * dtinv * w_nat[2]);

        // face 1
        face_fac = 0.5 * dx_dxi_1;
        add_value(values, col_b  , face_fac * 3./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_e  , face_fac * 3./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_ee , face_fac * 3./4. * -con_fac * dtinv * w_nat[2]);
        add_value(values, col_nb , face_fac * 1./4. * -con_fac * dtinv * w_nat[0]);
        add_value(values, col_ne , face_fac * 1./4. * -con_fac * dtinv * w_nat[1]);
        add_value(values, col_nee, face_fac * 1./4. * -con_fac * dtinv * w_nat[2]);

        // Contribution Delta q
        // face 0
        face_fac = 0.5;
        add_value(values, col_sb  + 1, face_fac * 1./4. * -con_fac * -theta);
        add_value(values, col_se  + 1, face_fac * 1./4. * -con_fac *  theta);
        add_value(values, col_see + 1, face_fac * 1./4. * 0.0);
        add_value(values, col_b   + 1, face_fac * 3./4. * -con_fac * -theta);
        add_value(values, col_e   + 1, face_fac * 3./4. * -con_fac *  theta);
        add_value(values, col_ee  + 1, face_fac * 3./4. * 0.0);

        // face 1
        face_fac = 0.5;
        add_value(values, col_b   + 1, face_fac * 3./4. * -con_fac * -theta);
        add_value(values, col_e   + 1, face_fac * 3./4. * -con_fac *  theta);
        add_value(values, col_ee  + 1, face_fac * 3./4. * 0.0);
        add_value(values, col_nb  + 1, face_fac * 1./4. * -con_fac * -theta);
        add_value(values, col_ne  + 1, face_fac * 1./4. * -con_fac *  theta);
        add_value(values, col_nee + 1, face_fac * 1./4. * 0.0);

        // Contribution Delta r
        add_value(values, col_b  + 2, 0.0);
        add_value(values, col_e  + 2, 0.0);
        add_value(values, col_ee + 2, 0.0);
               
        double dhdt_s = dtinv * (hp_s - hn_s);
        double dhdt_b = dtinv * (hp_b - hn_b);
        double dhdt_n = dtinv * (hp_n - hn_n);
        double dhdt_0 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_s);
        double dhdt_1 = 0.25 * (3.0 * dhdt_b + 1.0 * dhdt_n);
        double dqdx_0 = 0.25 * (3.0 * dqdx_b + 1.0 * dqdx_s);
        double dqdx_1 = 0.25 * (3.0 * dqdx_b + 1.0 * dqdx_n);

        rhs[row + 1] += - ( - 0.5 * con_fac * (dx_dxi_0 * dhdt_0 + dqdx_0) 
                            - 0.5 * con_fac * (dx_dxi_1 * dhdt_1 + dqdx_1) );
        //------------------------------------------------------------------
        // third row 
        //------------------------------------------------------------------
        col_sb  = r_eq + 0 * 3;
        col_se  = r_eq + 3 * 3;
        col_see = r_eq + 6 * 3;
        col_b   = r_eq + 1 * 3;
        col_e   = r_eq + 4 * 3;
        col_ee  = r_eq + 7 * 3;
        col_nb  = r_eq + 2 * 3;
        col_ne  = r_eq + 5 * 3;
        col_nee = r_eq + 8 * 3;
        //
        // Contribution Delta h
        add_value(values, col_b , 0.0);
        add_value(values, col_e , 0.0);
        add_value(values, col_ee, 0.0);
        // Contribution Delta q
        add_value(values, col_b  + 1, 0.0);
        add_value(values, col_e  + 1, 0.0);
        add_value(values, col_ee + 1, 0.0);
        // Contribution Delta r
        add_value(values, col_b  + 2, 1.0);
        add_value(values, col_e  + 2, 0.0);
        add_value(values, col_ee + 2, 0.0);
        //
        rhs[row + 2] = 0.0;
    }
    return 0;
}
void  molecule(std::vector<size_t>& p, size_t p_sw, size_t ny)
{
        p[0] = p_sw;
        p[1] = p_sw + 1;
        p[2] = p_sw + 2;
        p[3] = p_sw + ny;
        p[4] = p_sw + ny + 1;
        p[5] = p_sw + ny + 2;
        p[6] = p_sw + 2 * ny;
        p[7] = p_sw + 2 * ny + 1;
        p[8] = p_sw + 2 * ny + 2;
}
inline void add_value(double * values, size_t col, double data)
{ 
    values[col] += data; 
}
inline size_t ma_index(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}

