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
int boundary_north(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & dyinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_NORTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    memset(&values[c_eq], 0, 3 * 27 * sizeof(double));  // set all coefficients for one row of c-, q- and r-equation to zero

    int p_5 = c_eq/(3*27);  // node number of boundary point, ie north point of molecule
    int p_4 = p_5 - 1;
    int p_3 = p_5 - 2;
    int p_0 = p_5 - ny - 2;
    int p_1 = p_5 - ny - 1;
    int p_2 = p_5 - ny;
    int p_6 = p_5 + ny - 2;
    int p_7 = p_5 + ny - 1;
    int p_8 = p_5 + ny;

    double hn_b   = hn[p_5];
    double hn_jm1 = hn[p_4];
    double hn_jm2 = hn[p_3];
    double hp_b   = hp[p_5];
    double hp_jm1 = hp[p_4];
    double hp_jm2 = hp[p_3];
    double htheta_b   = htheta[p_5];
    double htheta_jm1 = htheta[p_4];
    double htheta_jm2 = htheta[p_3];

    double hn_jm12 = w_nat[0] * hn_b + w_nat[1] * hn_jm1 + w_nat[2] * hn_jm2;
    double hp_jm12 = w_nat[0] * hp_b + w_nat[1] * hp_jm1 + w_nat[2] * hp_jm2;
    double htheta_jm12 = w_nat[0] * htheta_b + w_nat[1] * htheta_jm1 + w_nat[2] * htheta_jm2;

    double qn_b   = qn[p_5];
    double qn_jm1 = qn[p_4];
    double qn_jm2 = qn[p_3];
    double qp_b   = qp[p_5];
    double qp_jm1 = qp[p_4];
    double qp_jm2 = qp[p_3];
    double qtheta_b   = qtheta[p_5];
    double qtheta_jm1 = qtheta[p_4];
    double qtheta_jm2 = qtheta[p_3];

    double qn_jm12 = w_nat[0] * qn_b + w_nat[1] * qn_jm1 + w_nat[2] * qn_jm2;
    double qp_jm12 = w_nat[0] * qp_b + w_nat[1] * qp_jm1 + w_nat[2] * qp_jm2;
    double qtheta_jm12 = w_nat[0] * qtheta_b + w_nat[1] * qtheta_jm1 + w_nat[2] * qtheta_jm2;

    double rn_b   = rn[p_5];
    double rn_jm1 = rn[p_4];
    double rn_jm2 = rn[p_3];
    double rp_b   = rp[p_5];
    double rp_jm1 = rp[p_4];
    double rp_jm2 = rp[p_3];
    double rtheta_b   = rtheta[p_5];
    double rtheta_jm1 = rtheta[p_4];
    double rtheta_jm2 = rtheta[p_3];

    double rn_jm12 = w_nat[0] * rn_b + w_nat[1] * rn_jm1 + w_nat[2] * rn_jm2;
    double rp_jm12 = w_nat[0] * rp_b + w_nat[1] * rp_jm1 + w_nat[2] * rp_jm2;
    double rtheta_jm12 = w_nat[0] * rtheta_b + w_nat[1] * rtheta_jm1 + w_nat[2] * rtheta_jm2;

    double zb_jm12 = w_nat[0] * zb[p_5] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_3];

    double sign = 1.0;
    double h_given = bc[BC_NORTH] - zb_jm12;
    double h_infty = -zb_jm12;
    double c_wave = std::sqrt(g * htheta_jm12);
    double con_fac = c_wave;

    // nnorth
    if (bc_type[BC_NORTH] == "dirichlet")
    {
        // Contribution Delta h
        int col_b  = c_eq + 5 * 3;
        int col_s  = c_eq + 4 * 3;
        int col_ss = c_eq + 3 * 3;
        //
        values[col_b ] =  1.0 * theta;
        values[col_s ] = -1.0 * theta;
        values[col_ss] = 0.0;
        rhs[row] = -( htheta_b + zb[p_5] - htheta_jm1 - zb[p_4] );
        // Contribution Delta q
        col_b  = q_eq + 5 * 3;
        col_s  = col_b - 3;
        col_ss = col_s - 3;
        values[col_b  + 1] = 1.0 * theta;
        values[col_s  + 1] =-1.0 * theta;
        values[col_ss + 1] = 0.0;
        rhs[row + 1] = -( qtheta_b - qtheta_jm1 );
        // Contribution Delta r0.5
        col_b  = r_eq + 5 * 3;
        col_s  = col_b - 3;
        col_ss = col_s - 3;
        values[col_b  + 2]  = 0.5 * theta;
        values[col_s  + 2]  = 0.5 * theta;
        values[col_ss + 2]  = 0.0;
        rhs[row + 2] = 0.0;
    }
    else
    {
        if (bc_type[BC_NORTH] == "mooiman")
        {
            // essential boundary condition
            // ----------------------------
            // first row
            int col_wb  = c_eq + 2 * 3;
            int col_ws  = c_eq + 1 * 3;
            int col_wss = c_eq + 0 * 3;
            int col_b   = c_eq + 5 * 3;
            int col_s   = c_eq + 4 * 3;
            int col_ss  = c_eq + 3 * 3;
            int col_eb  = c_eq + 8 * 3;
            int col_es  = c_eq + 7 * 3;
            int col_ess = c_eq + 6 * 3;
            //
            // face_0 
            double face_fac = 0.5 * dx * 0.25;
            set_value(values, col_wb , face_fac * 1.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_ws , face_fac * 1.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_wss, face_fac * 1.0 * theta * w_ess[2] * -c_wave);

            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_s , face_fac * 3.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_ss, face_fac * 3.0 * theta * w_ess[2] * -c_wave);

            // face 1
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_s , face_fac * 3.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_ss, face_fac * 3.0 * theta * w_ess[2] * -c_wave);

            set_value(values, col_eb , face_fac * 1.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_es , face_fac * 1.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_ess, face_fac * 1.0 * theta * w_ess[2] * -c_wave);

            // Contribution Delta q
            set_value(values, col_b  + 1, 0.0);
            set_value(values, col_s  + 1, 0.0);
            set_value(values, col_ss + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_wb  + 2, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_ws  + 2, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_wss + 2, face_fac * 1.0 * theta * w_ess[2]);

            set_value(values, col_b  + 2, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_s  + 2, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ss + 2, face_fac * 3.0 * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_b  + 2, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_s  + 2, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ss + 2, face_fac * 3.0 * theta * w_ess[2]);

            set_value(values, col_eb  + 2, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_es  + 2, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_ess + 2, face_fac * 1.0 * theta * w_ess[2]);
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

            rhs[row] = -dx * ( rt_jm12 - c_wave * (ht_jm12 - h_given) );
            
            if (bc_vars[BC_NORTH] == "zeta")
            {
                rhs[row] += - dx * 2. * c_wave * bc[BC_NORTH];
            }
            if (bc_vars[BC_NORTH] == "q")
            {
                rhs[row] += - dx * 2. * bc[BC_NORTH];
            }
        }
        else if (bc_type[BC_NORTH] == "borsboom")
        {
            //
            // Essential boundary condition
            // first row
            // ----------------------------
            int col_b  = c_eq + 5 * 3;
            int col_s  = c_eq + 4 * 3;
            int col_ss = c_eq + 3 * 3;
            //
            if (do_convection) { con_fac = c_wave - rp_jm12 / hp_jm12; }
            //
            set_value(values, col_b , dtinv * -con_fac * w_ess[0]);
            set_value(values, col_s , dtinv * -con_fac * w_ess[1]);
            set_value(values, col_ss, dtinv * -con_fac * w_ess[2]);
            //
            set_value(values, col_b  + 1, 0.0);
            set_value(values, col_s  + 1, 0.0);
            set_value(values, col_ss + 1, 0.0);
            //
            set_value(values, col_b  + 2, dtinv * w_ess[0]);
            set_value(values, col_s  + 2, dtinv * w_ess[1]);
            set_value(values, col_ss + 2, dtinv * w_ess[2]);
            //
            double dhdt = dtinv * (hp_b - hn_b) * w_ess[0]
                + dtinv * (hp_jm1 - hn_jm1) * w_ess[1]
                + dtinv * (hp_jm2 - hn_jm2) * w_ess[2];
            double drdt = dtinv * (rp_b - rn_b) * w_ess[0]
                + dtinv * (rp_jm1 - rn_jm1) * w_ess[1]
                + dtinv * (rp_jm2 - rn_jm2) * w_ess[2];
            rhs[row] = -(drdt - con_fac * dhdt);

            double corr_term = 0.0;
            if (bc_vars[BC_NORTH] == "zeta")
            {
                if (stationary) { sign = -1.0; }
                set_value(values, col_b , dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_s , dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ss, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double htheta_jm12 = w_ess[0] * htheta_b + w_ess[1] * htheta_jm1 + w_ess[2] * htheta_jm2;
                corr_term = - dhdt - ( eps_bc_corr * ((bc[BC_NORTH] - zb_jm12) - htheta_jm12) );
                rhs[row] += corr_term;
                sign = 1.0;
            }
            if (bc_vars[BC_NORTH] == "q")
            {
                if (stationary) { sign = -1.0; }
                set_value(values, col_b  + 2, dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_s  + 2, dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ss + 2, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double rtheta_jm12 = w_ess[0] * rtheta_b + w_ess[1] * rtheta_jm1 + w_ess[2] * rtheta_jm2;
                corr_term = - drdt + sign * eps_bc_corr * (bc[BC_NORTH]- qtheta_jm12);
                rhs[row] += corr_term;
                sign = 1.0;
            }
        }
        //------------------------------------------------------------------
        // natural boundary condition
        // -----------------------------------------------------------------
        double dhdt = dtinv * (hp_b - hn_b) * w_nat[0]
            + dtinv * (hp_jm1 - hn_jm1) * w_nat[1]
            + dtinv * (hp_jm2 - hn_jm2) * w_nat[2];
        double drdy = dyinv * (rtheta_b - rtheta_jm1);
        double drdt = dtinv * (rp_b - rn_b) * w_nat[0]
            + dtinv * (rp_jm1 - rn_jm1) * w_nat[1]
            + dtinv * (rp_jm2 - rn_jm2) * w_nat[2];
        double dzetady = dyinv * (htheta_b + zb[p_5] - htheta_jm1 - zb[p_4]);
        //
        // second equation
        // q-momentum (tangential equation, q == 0)
        //
        int col_b  = q_eq + 5 * 3;
        int col_s  = q_eq + 4 * 3;
        int col_ss = q_eq + 3 * 3;
        //
        set_value(values, col_b , 0.0);
        set_value(values, col_s , 0.0);
        set_value(values, col_ss, 0.0);
        //
        set_value(values, col_b  + 1, 1.0);
        set_value(values, col_s  + 1, 0.0);
        set_value(values, col_ss + 1, 0.0);
        //
        set_value(values, col_b  + 2, 0.0);
        set_value(values, col_s  + 2, 0.0);
        set_value(values, col_ss + 2, 0.0);
        //
        rhs[row + 1] = 0.0;
        //
        // third equation
        // momentum part dr/dt + gh d(zeta)/dy
        //
        col_b  = r_eq + 5 * 3;
        col_s  = r_eq + 4 * 3;
        col_ss = r_eq + 3 * 3;
        //
        // Contribution Delta h
        set_value(values, col_b , w_nat[0] * theta * g * dzetady + dyinv * theta * g * htheta_jm12);
        set_value(values, col_s , w_nat[1] * theta * g * dzetady - dyinv * theta * g * htheta_jm12);
        set_value(values, col_ss, w_nat[2] * theta * g * dzetady);
        // Contribution Delta q
        set_value(values, col_b  + 1, 0.0);
        set_value(values, col_s  + 1, 0.0);
        set_value(values, col_ss + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, dtinv * w_nat[0]);
        set_value(values, col_s  + 2, dtinv * w_nat[1]);
        set_value(values, col_ss + 2, dtinv * w_nat[2]);
                        //
        rhs[row + 2] = -( drdt + g * htheta_jm12 * dzetady );
        //
        if (do_convection)
        {
        }
        if (do_bed_shear_stress)
        {
        }
        if (do_viscosity)
        {
        }
        //
        // continuity part +c_wave * (dhdt + dq/dx + dr/dy)
        //
        if (do_convection) { con_fac = c_wave + rp_jm12 / hp_jm12; }
        // Contribution Delta h
        set_value(values, col_b , con_fac * dtinv * w_nat[0]);
        set_value(values, col_s , con_fac * dtinv * w_nat[1]);
        set_value(values, col_ss, con_fac * dtinv * w_nat[2]);
        // Contribution Delta q
        set_value(values, col_b  + 1, 0.0);
        set_value(values, col_s  + 1, 0.0);
        set_value(values, col_ss + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, con_fac *  dyinv * theta);
        set_value(values, col_s  + 2, con_fac * -dyinv * theta);
        set_value(values, col_ss + 2, 0.0);
        //
        rhs[row + 2] += - ( con_fac * (dhdt + drdy) );
    }
    return 0;
}
//==============================================================================
int boundary_east(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & dxinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_EAST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    memset(&values[c_eq], 0, 3 * 27 * sizeof(double));  // set all coefficients for one row of c-, q- and r-equation to zero

    int p_7 = c_eq/(3*27);  // node number of boundary point, ie east point of molecule
    int p_8 = p_7 + 1;
    int p_6 = p_7 - 1;
    int p_5 = p_7 - ny + 1;
    int p_4 = p_7 - ny;
    int p_3 = p_7 - ny - 1;
    int p_2 = p_7 - 2 * ny + 1;
    int p_1 = p_7 - 2 * ny;
    int p_0 = p_7 - 2 * ny - 1;

    double hn_b   = hn[p_7];
    double hn_im1 = hn[p_4];
    double hn_im2 = hn[p_1];
    double hp_b   = hp[p_7];
    double hp_im1 = hp[p_4];
    double hp_im2 = hp[p_1];
    double htheta_b   = htheta[p_7];
    double htheta_im1 = htheta[p_4];
    double htheta_im2 = htheta[p_1];

    double hn_im12 = w_nat[0] * hn_b + w_nat[1] * hn_im1 + w_nat[2] * hn_im2;
    double hp_im12 = w_nat[0] * hp_b + w_nat[1] * hp_im1 + w_nat[2] * hp_im2;
    double htheta_im12 = w_nat[0] * htheta_b + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;

    double qn_b   = qn[p_7];
    double qn_im1 = qn[p_4];
    double qn_im2 = qn[p_1];
    double qp_b   = qp[p_7];
    double qp_im1 = qp[p_4];
    double qp_im2 = qp[p_1];
    double qtheta_b   = qtheta[p_7];
    double qtheta_im1 = qtheta[p_4];
    double qtheta_im2 = qtheta[p_1];

    double qn_im12 = w_nat[0] * qn_b + w_nat[1] * qn_im1 + w_nat[2] * qn_im2;
    double qp_im12 = w_nat[0] * qp_b + w_nat[1] * qp_im1 + w_nat[2] * qp_im2;
    double qtheta_im12 = w_nat[0] * qtheta_b + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;

    double rn_b   = rn[p_7];
    double rn_im1 = rn[p_4];
    double rn_im2 = rn[p_1];
    double rp_b   = rp[p_7];
    double rp_im1 = rp[p_4];
    double rp_im2 = rp[p_1];
    double rtheta_b   = rtheta[p_7];
    double rtheta_im1 = rtheta[p_4];
    double rtheta_im2 = rtheta[p_1];

    double rn_im12 = w_nat[0] * rn_b + w_nat[1] * rn_im1 + w_nat[2] * rn_im2;
    double rp_im12 = w_nat[0] * rp_b + w_nat[1] * rp_im1 + w_nat[2] * rp_im2;
    double rtheta_im12 = w_nat[0] * rtheta_b + w_nat[1] * rtheta_im1 + w_nat[2] * rtheta_im2;

    double zb_im12 = w_nat[0] * zb[p_7] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_1];

    double sign = 1.0;
    double h_given = bc[BC_EAST] - zb_im12;
    double h_infty = -zb_im12;
    double c_wave = std::sqrt(g * htheta_im12);
    double con_fac = c_wave;

    // eeast
    if (bc_type[BC_EAST] == "dirichlet")
    {
        // Contribution Delta h
        int col_b  = c_eq + 7 * 3;
        int col_w  = c_eq + 4 * 3;
        int col_ww = c_eq + 1 * 3;
        //
        values[col_b ] =  1.0 * theta;
        values[col_w ] = -1.0 * theta;
        values[col_ww] = 0.0;
        rhs[row] = -( htheta_b + zb[p_7] - htheta_im1 - zb[p_4]);
        // Contribution Delta q
        col_b  = q_eq + 7 * 3;
        col_w  = col_b - 9;
        col_ww = col_w - 9;
        values[col_b  + 1] = 0.5 * theta;
        values[col_w  + 1] = 0.5 * theta;
        values[col_ww + 1] = 0.0;
        rhs[row + 1] = 0.0;
        // Contribution Delta r0.5
        col_b  = r_eq + 7 * 3;
        col_w  = col_b - 9;
        col_ww = col_w - 9;
        values[col_b  + 2]  =  1.0 * theta;
        values[col_w  + 2]  = -1.0 * theta;
        values[col_ww + 2]  = 0.0;
        rhs[row + 2] = -( rtheta_b - rtheta_im1 );
    }
    else
    {
        if (bc_type[BC_EAST] == "mooiman")
        {
            // essential boundary condition
            // ----------------------------
            // first row
            int col_nb  = c_eq + 8 * 3;
            int col_nw  = c_eq + 5 * 3;
            int col_nww = c_eq + 2 * 3;
            int col_b  = c_eq + 7 * 3;
            int col_w  = c_eq + 4 * 3;
            int col_ww = c_eq + 1 * 3;
            int col_sb  = c_eq + 6 * 3;
            int col_sw  = c_eq + 3 * 3;
            int col_sww = c_eq + 0 * 3;
            //
            // face_0 
            double face_fac = 0.5 * dy * 0.25;
            set_value(values, col_nb , face_fac * 1.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_nw , face_fac * 1.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_nww, face_fac * 1.0 * theta * w_ess[2] * -c_wave);
            
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_w , face_fac * 3.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_ww, face_fac * 3.0 * theta * w_ess[2] * -c_wave);
            
             // face 1
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_w , face_fac * 3.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_ww, face_fac * 3.0 * theta * w_ess[2] * -c_wave);
           
            set_value(values, col_sb , face_fac * 1.0 * theta * w_ess[0] * -c_wave);
            set_value(values, col_sw , face_fac * 1.0 * theta * w_ess[1] * -c_wave);
            set_value(values, col_sww, face_fac * 1.0 * theta * w_ess[2] * -c_wave);
            
            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_nb  + 1, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_nw  + 1, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_nww + 1, face_fac * 1.0 * theta * w_ess[2]);

            set_value(values, col_b  + 1, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_w  + 1, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ww + 1, face_fac * 3.0 * theta * w_ess[2]);
            
            // face 1
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_b  + 1, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_w  + 1, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ww + 1, face_fac * 3.0 * theta * w_ess[2]);
            
            set_value(values, col_sb  + 1, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_sw  + 1, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_sww + 1, face_fac * 1.0 * theta * w_ess[2]);
            
            // Contribution Delta r
            set_value(values, col_b  + 2, 0.0);
            set_value(values, col_w  + 2, 0.0);
            set_value(values, col_ww + 2, 0.0);
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
            
            rhs[row] = -dy * (qt_im12 - c_wave * (ht_im12 - h_given));
            if (bc_vars[BC_EAST] == "zeta")
            {
                rhs[row] += -dy * 2. * c_wave * bc[BC_EAST];
            }
            if (bc_vars[BC_EAST] == "q")
            {
                rhs[row] += -dy * 2. * bc[BC_EAST];
            }
        }
        else if (bc_type[BC_EAST] == "borsboom")
        {
            //
            // Essential boundary condition
            // first row
            // ----------------------------
            int col_b  = c_eq + 7 * 3;
            int col_w  = c_eq + 4 * 3;
            int col_ww = c_eq + 1 * 3;
            //
            if (do_convection) { con_fac = c_wave - qp_im12 / hp_im12; }
            //
            // face_0 [(i-1/2, j-1/2), (i-1/2,j)]
            // face_1 [(i-1/2, j), (i-1/2,j+1/2)]
            set_value(values, col_b , dtinv * -con_fac * w_ess[0]);
            set_value(values, col_w , dtinv * -con_fac * w_ess[1]);
            set_value(values, col_ww, dtinv * -con_fac * w_ess[2]);
            // flow x
            set_value(values, col_b  + 1, dtinv * w_ess[0]);
            set_value(values, col_w  + 1, dtinv * w_ess[1]);
            set_value(values, col_ww + 1, dtinv * w_ess[2]);
            // flow y
            set_value(values, col_b  + 2, 0.0);
            set_value(values, col_w  + 2, 0.0);
            set_value(values, col_ww + 2, 0.0);
            //
            double dhdt = dtinv * (hp_b - hn_b) * w_ess[0]
                + dtinv * (hp_im1 - hn_im1) * w_ess[1]
                + dtinv * (hp_im2 - hn_im2) * w_ess[2];
            double dqdt = dtinv * (qp_b - qn_b) * w_ess[0]
                + dtinv * (qp_im1 - qn_im1) * w_ess[1]
                + dtinv * (qp_im2 - qn_im2) * w_ess[2];
            rhs[row] = -(dqdt - con_fac * dhdt);

            double corr_term = 0.0;
            if (bc_vars[BC_EAST] == "zeta")
            {
                if (stationary) { sign = -1.0; }
                set_value(values, col_b , dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_w , dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ww, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double htheta_im12 = w_ess[0] * htheta_b + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2;
                corr_term = - dhdt - ( eps_bc_corr * ( (bc[BC_EAST] - zb_im12) - htheta_im12) );
                rhs[row] += corr_term;
                sign = 1.0;
            }
            if (bc_vars[BC_EAST] == "q")
            {
                if (stationary) { sign = -1.0; }
                set_value(values, col_b  + 1, dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_w  + 1, dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ww + 1, dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2]);

                double qtheta_im12 = w_ess[0] * qtheta_b + w_ess[1] * qtheta_im1 + w_ess[2] * qtheta_im2;
                corr_term = - dqdt + sign * eps_bc_corr * (bc[BC_EAST] - qtheta_im12);
                rhs[row] += corr_term;
                sign = 1.0;
            }
        }
        //
        // natural boundary condition
        //  second row
        // --------------------------
        //
        double dhdt = dtinv * (hp_b - hn_b) * w_nat[0]
                + dtinv * (hp_im1 - hn_im1) * w_nat[1]
                + dtinv * (hp_im2 - hn_im2) * w_nat[2];
        double dqdx = dxinv * (qtheta_b - qtheta_im1);
        double dqdt = dtinv * (qp_b - qn_b) * w_nat[0]
                + dtinv * (qp_im1 - qn_im1) * w_nat[1]
                + dtinv * (qp_im2 - qn_im2) * w_nat[2];
        double dzetadx = dxinv * (htheta_b + zb[p_7] - htheta_im1 - zb[p_4]);
        //
        // second equation
        // momentum part dq/dt + gh d(zeta)/dx
        //
        int col_b  = q_eq + 7 * 3;
        int col_w  = q_eq + 4 * 3;
        int col_ww = q_eq + 1 * 3;
        //
        set_value(values, col_b , w_nat[0] * theta * g * dzetadx + dxinv * theta * g * htheta_im12);
        set_value(values, col_w , w_nat[1] * theta * g * dzetadx - dxinv * theta * g * htheta_im12);
        set_value(values, col_ww, w_nat[2] * theta * g * dzetadx);
        //
        set_value(values, col_b  + 1, dtinv * w_nat[0]);
        set_value(values, col_w  + 1, dtinv * w_nat[1]);
        set_value(values, col_ww + 1, dtinv * w_nat[2]);
        // flow y
        set_value(values, col_b  + 2, 0.0);
        set_value(values, col_w  + 2, 0.0);
        set_value(values, col_ww + 2, 0.0);
        //
        rhs[row + 1] = -(dqdt + g * htheta_im12 * dzetadx);
        //
        if (do_convection)
        {
            double nxi = -1.0;
            double neta = 0.0;
            //double h_jm14 = h_0 * + h_w * + h_ww + h_jm1 * + h_im1j * + h_im2j *  + h_im2j * ;
            double h_jm14 = 0.0;
            double q_jm14 = 0.0;
            double r_jm14 = 0.0;
            double h_jp14 = 0.0;
            double q_jp14 = 0.0;
            double r_jp14 = 0.0;
            double qp_1 = 0.5 * dy * (convection_J_10(h_jm14, q_jm14, r_jm14, nxi, neta));
            double qp_2 = 0.5 * dy * (convection_J_10(h_jp14, q_jp14, r_jp14, nxi, neta));
            rhs[row + 1] += qp_1 + qp_2;
        }
        if (do_bed_shear_stress)
        {
        }
        if (do_viscosity)
        {
        }
        //
        // continuity part +c_wave * (dhdt + dq/dx + dr/dy)
        //
        if (do_convection) { con_fac = c_wave + qp_im12 / hp_im12; }
        // Contribution Delta h
        set_value(values, col_b , con_fac * dtinv * w_nat[0]);
        set_value(values, col_w , con_fac * dtinv * w_nat[1]);
        set_value(values, col_ww, con_fac * dtinv * w_nat[2]);
        // Contribution Delta q
        set_value(values, col_b  + 1, con_fac * dxinv * theta);
        set_value(values, col_w  + 1, con_fac * -dxinv * theta);
        set_value(values, col_ww + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, 0.0);
        set_value(values, col_w  + 2, 0.0);
        set_value(values, col_ww + 2, 0.0);
        //
        rhs[row + 1] += - ( con_fac * (dhdt + dqdx) );
        //
        // third equation
        // r-momentum (tangential equation, r == 0)
        //
        col_b  = r_eq + 7 * 3;
        col_w  = r_eq + 4 * 3;
        col_ww = r_eq + 1 * 3;
        // 
        // Contribution Delta h
        set_value(values, col_b , 0.0);
        set_value(values, col_w , 0.0);
        set_value(values, col_ww, 0.0);
        // Contribution Delta q
        set_value(values, col_b  + 1, 0.0);
        set_value(values, col_w  + 1, 0.0);
        set_value(values, col_ww + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, 1.0);
        set_value(values, col_w  + 2, 0.0);
        set_value(values, col_ww + 2, 0.0);
        //
        rhs[row + 2] = 0.0;
    }
    return 0;
}
//==============================================================================
int boundary_south(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & dyinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_SOUTH, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    memset(&values[c_eq], 0, 3 * 27 * sizeof(double));  // set all coefficients for one row of c-, q- and r-equation to zero

    int p_3 = c_eq/(3*27);  // node number of boundary point, ie south point of molecule
    int p_4 = p_3 + 1;
    int p_5 = p_3 + 2;
    int p_0 = p_3 - ny;
    int p_1 = p_3 - ny + 1;
    int p_2 = p_3 - ny + 2;
    int p_6 = p_3 + ny;
    int p_7 = p_3 + ny + 1;
    int p_8 = p_3 + ny + 2;

    double hn_b   = hn[p_3];
    double hn_jp1 = hn[p_4];
    double hn_jp2 = hn[p_5];
    double hp_b   = hp[p_3];
    double hp_jp1 = hp[p_4];
    double hp_jp2 = hp[p_5];
    double htheta_b   = htheta[p_3];
    double htheta_jp1 = htheta[p_4];
    double htheta_jp2 = htheta[p_5];

    double hn_jp12 = w_nat[0] * hn_b + w_nat[1] * hn_jp1 + w_nat[2] * hn_jp2;
    double hp_jp12 = w_nat[0] * hp_b + w_nat[1] * hp_jp1 + w_nat[2] * hp_jp2;
    double htheta_jp12 = w_nat[0] * htheta_b + w_nat[1] * htheta_jp1 + w_nat[2] * htheta_jp2;

    double qn_b   = qn[p_3];
    double qn_jp1 = qn[p_4];
    double qn_jp2 = qn[p_5];
    double qp_b   = qp[p_3];
    double qp_jp1 = qp[p_4];
    double qp_jp2 = qp[p_5];
    double qtheta_b   = qtheta[p_3];
    double qtheta_jp1 = qtheta[p_4];
    double qtheta_jp2 = qtheta[p_5];

    double qn_jp12 = w_nat[0] * qn_b + w_nat[1] * qn_jp1 + w_nat[2] * qn_jp2;
    double qp_jp12 = w_nat[0] * qp_b + w_nat[1] * qp_jp1 + w_nat[2] * qp_jp2;
    double qtheta_jp12 = w_nat[0] * qtheta_b + w_nat[1] * qtheta_jp1 + w_nat[2] * qtheta_jp2;

    double rn_b   = rn[p_3];
    double rn_jp1 = rn[p_4];
    double rn_jp2 = rn[p_5];
    double rp_b   = rp[p_3];
    double rp_jp1 = rp[p_4];
    double rp_jp2 = rp[p_5];
    double rtheta_b   = rtheta[p_3];
    double rtheta_jp1 = rtheta[p_4];
    double rtheta_jp2 = rtheta[p_5];

    double rn_jp12 = w_nat[0] * rn_b + w_nat[1] * rn_jp1 + w_nat[2] * rn_jp2;
    double rp_jp12 = w_nat[0] * rp_b + w_nat[1] * rp_jp1 + w_nat[2] * rp_jp2;
    double rtheta_jp12 = w_nat[0] * rtheta_b + w_nat[1] * rtheta_jp1 + w_nat[2] * rtheta_jp2;

    double zb_jp12 = w_nat[0] * zb[p_3] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_5];
    double h_given = bc[BC_SOUTH] - zb_jp12;
    double h_infty = -zb_jp12;
    double c_wave = std::sqrt(g * htheta_jp12);
    double con_fac = c_wave;

    // ssouth
    if (bc_type[BC_SOUTH] == "dirichlet")
    {
        // Contribution Delta h
        int col_b  = c_eq + 3 * 3;
        int col_n  = c_eq + 4 * 3;
        int col_nn = c_eq + 5 * 3;
        //
        values[col_b ] = -1.0 * theta;
        values[col_n ] =  1.0 * theta;
        values[col_nn] = 0.0;
        rhs[row] = -( htheta_jp1 + zb[p_4] - htheta_b - zb[p_3] );
        // Contribution Delta q
        col_b  = q_eq + 3 * 3;
        col_n  = col_b + 3;
        col_nn = col_n + 3;
        values[col_b  + 1] = -1.0 * theta;
        values[col_n  + 1] =  1.0 * theta;
        values[col_nn + 1] = 0.0;
        rhs[row + 1] = -( qtheta_jp1 - qtheta_b );
        // Contribution Delta r0.5
        col_b  = r_eq + 3 * 3;
        col_n  = col_b + 3;
        col_nn = col_n + 3;
        values[col_b  + 2]  = 0.5 * theta;
        values[col_n  + 2]  = 0.5 * theta;
        values[col_nn + 2]  = 0.0;
        rhs[row + 2] = 0.0;
    }
    else
    {
        if (bc_type[BC_SOUTH] == "mooiman")
        {
            // Essential boundary condition
            // ----------------------------
            // first row
            int col_wb  = c_eq + 0* 3;
            int col_wn  = c_eq + 1 * 3;
            int col_wnn = c_eq + 2 * 3;
            int col_b  = c_eq + 3 * 3;
            int col_n  = c_eq + 4 * 3;
            int col_nn = c_eq + 5 * 3;
            int col_eb  = c_eq + 6 * 3;
            int col_en  = c_eq + 7 * 3;
            int col_enn = c_eq + 8 * 3;
            //
            // face_0 
            double face_fac = 0.5 * dx * 0.25;
            set_value(values, col_wb , face_fac * 1.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_wn , face_fac * 1.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_wnn, face_fac * 1.0 * theta * w_ess[2] * c_wave);

            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_n , face_fac * 3.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_nn, face_fac * 3.0 * theta * w_ess[2] * c_wave);

            // face 1
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_n , face_fac * 3.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_nn, face_fac * 3.0 * theta * w_ess[2] * c_wave);

            set_value(values, col_eb , face_fac * 1.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_en , face_fac * 1.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_enn, face_fac * 1.0 * theta * w_ess[2] * c_wave);

            // Contribution Delta q
            set_value(values, col_b  + 1, 0.0);
            set_value(values, col_n  + 1, 0.0);
            set_value(values, col_nn + 1, 0.0);

            // Contribution Delta r
            // face 0
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_wb  + 2, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_wn  + 2, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_wnn + 2, face_fac * 1.0 * theta * w_ess[2]);

            set_value(values, col_b  + 2, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_n  + 2, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_nn + 2, face_fac * 3.0 * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dx * 0.25;
            set_value(values, col_b  + 2, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_n  + 2, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_nn + 2, face_fac * 3.0 * theta * w_ess[2]);

            set_value(values, col_eb  + 2, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_en  + 2, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_enn + 2, face_fac * 1.0 * theta * w_ess[2]);
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

            rhs[row] = - dx * ( rt_jp12 + c_wave * (ht_jp12 - h_given) );

            if (bc_vars[BC_SOUTH] == "zeta")
            {
                rhs[row] += dx * 2. * c_wave * bc[BC_SOUTH];
            }
            if (bc_vars[BC_SOUTH] == "q")
            {
                rhs[row] += dx * 2. * bc[BC_SOUTH];
            }
        }
        else if (bc_type[BC_SOUTH] == "borsboom")
        {
            //
            // Essential boundary condition
            // first row
            // ----------------------------
            //
            int col_b  = c_eq + 3 * 3;
            int col_n  = c_eq + 4 * 3;
            int col_nn = c_eq + 5 * 3;
            //
            if (do_convection) { con_fac = c_wave + rp_jp12 / hp_jp12; }
            //
            set_value(values, col_b , dtinv * con_fac * w_ess[0]);
            set_value(values, col_n , dtinv * con_fac * w_ess[1]);
            set_value(values, col_nn, dtinv * con_fac * w_ess[2]);
            // Contribution Delta q
            set_value(values, col_b  + 1, 0.0);
            set_value(values, col_n  + 1, 0.0);
            set_value(values, col_nn + 1, 0.0);
            // Contribution Delta r
            set_value(values, col_b  + 2, dtinv * w_ess[0]);
            set_value(values, col_n  + 2, dtinv * w_ess[1]);
            set_value(values, col_nn + 2, dtinv * w_ess[2]);
            //
            double dhdt = dtinv * (hp_b - hn_b) * w_ess[0]
                + dtinv * (hp_jp1 - hn_jp1) * w_ess[1]
                + dtinv * (hp_jp2 - hn_jp2) * w_ess[2];
            double drdt = dtinv * (rp_b - rn_b) * w_ess[0]
                + dtinv * (rp_jp1 - rn_jp1) * w_ess[1]
                + dtinv * (rp_jp2 - rn_jp2) * w_ess[2];
            rhs[row] = -(drdt + con_fac * dhdt);

            double corr_term = 0.0;
            if (bc_vars[BC_SOUTH] == "zeta")
            {
                set_value(values, col_b , dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_n , dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_nn, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double htheta_jp12 = w_ess[0] * htheta_b + w_ess[1] * htheta_jp1 + w_ess[2] * htheta_jp2;
                corr_term = + dhdt + ( eps_bc_corr * ((bc[BC_SOUTH] - zb_jp12) - htheta_jp12) );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_SOUTH] == "q")
            {
                set_value(values, col_b  + 2, dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_n  + 2, dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_nn + 2, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double rtheta_jp12 = w_ess[0] * rtheta_b + w_ess[1] * rtheta_jp1 + w_ess[2] * rtheta_jp2;
                corr_term = - drdt - eps_bc_corr * (bc[BC_SOUTH] - rtheta_jp12);
                rhs[row] += corr_term;
            }
        }
        //
        // natural boundary condition
        // second row
        // --------------------------
        //
        if (do_convection) { con_fac = c_wave - rp_jp12 / hp_jp12; }
        double dhdt = dtinv * (hp_b - hn_b) * w_nat[0]
            + dtinv * (hp_jp1 - hn_jp1) * w_nat[1]
            + dtinv * (hp_jp2 - hn_jp2) * w_nat[2];
        double drdy = dyinv * (rtheta_jp1 - rtheta_b);
        double drdt = dtinv * (rp_b - rn_b) * w_nat[0]
            + dtinv * (rp_jp1 - rn_jp1) * w_nat[1]
            + dtinv * (rp_jp2 - rn_jp2) * w_nat[2];
        double dzetady = dyinv * (htheta_jp1 + zb[p_4] - htheta_b - zb[p_3]);
        //
        // second equation
        // q-momentum (tangential equation, q == 0)
        //
        int col_b  = q_eq + 3 * 3;
        int col_n  = q_eq + 4 * 3;
        int col_nn = q_eq + 5 * 3;
        //
        set_value(values, col_b, 0.0);
        set_value(values, col_n, 0.0);
        set_value(values, col_nn, 0.0);
        //
        set_value(values, col_b  + 1, 1.0);
        set_value(values, col_n  + 1, 0.0);
        set_value(values, col_nn + 1, 0.0);
        //
        set_value(values, col_b  + 2, 0.0);
        set_value(values, col_n  + 2, 0.0);
        set_value(values, col_nn + 2, 0.0);
        //
        rhs[row + 1] = 0.0;
        //
        // third equation
        // momentum part dr/dt + gh d(zeta)/dy
        //
        col_b  = r_eq + 3 * 3;
        col_n  = r_eq + 4 * 3;
        col_nn = r_eq + 5 * 3;
        //
        // Contribution Delta h
        set_value(values, col_b , w_nat[0] * theta * g * dzetady - dyinv * theta * g * htheta_jp12);
        set_value(values, col_n , w_nat[1] * theta * g * dzetady + dyinv * theta * g * htheta_jp12);
        set_value(values, col_nn, w_nat[2] * theta * g * dzetady);
        // Contribution Delta q
        set_value(values, col_b  + 1, 0.0);
        set_value(values, col_n  + 1, 0.0);
        set_value(values, col_nn + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, dtinv * w_nat[0]);
        set_value(values, col_n  + 2, dtinv * w_nat[1]);
        set_value(values, col_nn + 2, dtinv * w_nat[2]);
        //
        rhs[row + 2] = -( drdt + g * htheta_jp12 * dzetady );
        //
        if (do_convection)
        {
        }
        if (do_bed_shear_stress)
        {
        }
        if (do_viscosity)
        {
        }
        //
        // continuity part -c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        set_value(values, col_b , -con_fac * dtinv * w_nat[0]);
        set_value(values, col_n , -con_fac * dtinv * w_nat[1]);
        set_value(values, col_nn, -con_fac * dtinv * w_nat[2]);
        // Contribution Delta q
        set_value(values, col_b  + 1, 0.0);
        set_value(values, col_n  + 1, 0.0);
        set_value(values, col_nn + 1, 0.0);
        // Contribution Delta r
        set_value(values, col_b  + 2, -con_fac * -dyinv * theta);
        set_value(values, col_n  + 2, -con_fac * dyinv * theta);
        set_value(values, col_nn + 2, 0.0);
        //
        rhs[row + 2] += -( -con_fac * (dhdt + drdy) );
    }
    return 0;
}
int boundary_west(double* values, int row, int c_eq, int q_eq, int r_eq, Eigen::VectorXd& rhs, 
    double & dtinv, double & dxinv, double & theta, double & g, double eps_bc_corr, 
    bool stationary, bool do_convection, bool do_bed_shear_stress, bool do_viscosity, 
    double dx, double dy, int nx, int ny,
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp,
    std::vector<double>& htheta, std::vector<double>& qtheta, std::vector<double>& rtheta,
    std::vector<double>& zb, std::vector<std::string> bc_type, std::vector<std::string> bc_vars, int BC_WEST, std::vector<double> bc,
    std::vector<double>& w_nat, std::vector<double>& w_ess)
{
    memset(&values[c_eq], 0, 3 * 27 * sizeof(double));  // set all coefficients for one row of c-, q- and r-equation to zero

    int p_1 = c_eq/(3*27);  // node number of boundary point, ie west point of molucule
    int p_0 = p_1 - 1;
    int p_2 = p_1 + 1;
    int p_3 = p_1 + ny - 1;
    int p_4 = p_1 + ny;
    int p_5 = p_1 + ny + 1;
    int p_6 = p_1 + 2 * ny - 1;
    int p_7 = p_1 + 2 * ny;
    int p_8 = p_1 + 2 * ny + 1;

    double hn_b = hn[p_1];
    double hn_ip1 = hn[p_4];
    double hn_ip2 = hn[p_7];
    double hp_b = hp[p_1];
    double hp_ip1 = hp[p_4];
    double hp_ip2 = hp[p_7];
    double htheta_b   = htheta[p_1];
    double htheta_ip1 = htheta[p_4];
    double htheta_ip2 = htheta[p_7];

    double hn_ip12 = w_nat[0] * hn_b + w_nat[1] * hn_ip1 + w_nat[2] * hn_ip2;
    double hp_ip12 = w_nat[0] * hp_b + w_nat[1] * hp_ip1 + w_nat[2] * hp_ip2;
    double htheta_ip12 = w_nat[0] * htheta_b + w_nat[1] * htheta_ip1 + w_nat[2] * htheta_ip2;

    double qn_b = qn[p_1];
    double qn_ip1 = qn[p_4];
    double qn_ip2 = qn[p_7];
    double qp_b = qp[p_1];
    double qp_ip1 = qp[p_4];
    double qp_ip2 = qp[p_7];
    double qtheta_b   = qtheta[p_1];
    double qtheta_ip1 = qtheta[p_4];
    double qtheta_ip2 = qtheta[p_7];

    double qn_ip12 = w_nat[0] * qn_b + w_nat[1] * qn_ip1 + w_nat[2] * qn_ip2;
    double qp_ip12 = w_nat[0] * qp_b + w_nat[1] * qp_ip1 + w_nat[2] * qp_ip2;
    double qtheta_ip12 = w_nat[0] * qtheta_b + w_nat[1] * qtheta_ip1 + w_nat[2] * qtheta_ip2;

    double rn_b = rn[p_1];
    double rn_ip1 = rn[p_4];
    double rn_ip2 = rn[p_7];
    double rp_b = rp[p_1];
    double rp_ip1 = rp[p_4];
    double rp_ip2 = rp[p_7];
    double rtheta_b   = rtheta[p_1];
    double rtheta_ip1 = rtheta[p_4];
    double rtheta_ip2 = rtheta[p_7];

    double rn_ip12 = w_nat[0] * rn_b + w_nat[1] * rn_ip1 + w_nat[2] * rn_ip2;
    double rp_ip12 = w_nat[0] * rp_b + w_nat[1] * rp_ip1 + w_nat[2] * rp_ip2;
    double rtheta_ip12 = w_nat[0] * rtheta_b + w_nat[1] * rtheta_ip1 + w_nat[2] * rtheta_ip2;

    double zb_ip12 = w_nat[0] * zb[p_1] + w_nat[1] * zb[p_4] + w_nat[2] * zb[p_7];

    double h_given = bc[BC_WEST] - zb_ip12;
    double h_infty = -zb_ip12;
    double c_wave = std::sqrt(g * htheta_ip12);
    double con_fac = c_wave;

    // wwest
    if (bc_type[BC_WEST] == "dirichlet")
    {
        // Contribution Delta h
        int col_b  = c_eq + 1 * 3; // point of boundary, ie west point of molecule
        int col_e  = c_eq + 4 * 3;
        int col_ee = c_eq + 7 * 3;
        //
        values[col_b ] = -1.0 * theta;
        values[col_e ] =  1.0 * theta;
        values[col_ee] = 0.0;
        rhs[row] = -( htheta_ip1 + zb[p_4] - htheta_b - zb[p_1] );
        // Contribution Delta q
        col_b  = q_eq + 1 * 3;
        col_e  = col_b + 9;
        col_ee = col_e + 9;
        values[col_b  + 1] = 0.5 * theta;
        values[col_e  + 1] = 0.5 * theta;
        values[col_ee + 1] = 0.0;
        rhs[row + 1] = 0.0;
        // Contribution Delta r
        col_b  = r_eq + 1 * 3;
        col_e  = col_b + 9;
        col_ee = col_e + 9;
        values[col_b  + 2]  = -1.0 * theta;
        values[col_e  + 2]  =  1.0 * theta;
        values[col_ee + 2]  = 0.0;
        rhs[row + 2] = -( rtheta_ip1 - rtheta_b);
    }
    else
    {
        if (bc_type[BC_WEST] == "mooiman")
        {
            // essential boundary condition
            // first row
            int col_sb  = c_eq + 0 * 3;
            int col_se  = c_eq + 3 * 3;
            int col_see = c_eq + 6 * 3;
            int col_b  = c_eq + 1 * 3;
            int col_e  = c_eq + 4 * 3;
            int col_ee = c_eq + 7 * 3;
            int col_nb  = c_eq + 2 * 3;
            int col_ne  = c_eq + 5 * 3;
            int col_nee = c_eq + 8 * 3;
            //
            // face 0
            double face_fac = 0.5 * dy * 0.25;
            set_value(values, col_sb , face_fac * 1.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_se , face_fac * 1.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_see, face_fac * 1.0 * theta * w_ess[2] * c_wave);
            
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_e , face_fac * 3.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_ee, face_fac * 3.0 * theta * w_ess[2] * c_wave);
            
            //face 1
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_b , face_fac * 3.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_e , face_fac * 3.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_ee, face_fac * 3.0 * theta * w_ess[2] * c_wave);
            
            set_value(values, col_nb , face_fac * 1.0 * theta * w_ess[0] * c_wave);
            set_value(values, col_ne , face_fac * 1.0 * theta * w_ess[1] * c_wave);
            set_value(values, col_nee, face_fac * 1.0 * theta * w_ess[2] * c_wave);
            
            // Contribution Delta q
            // face 0
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_sb  + 1, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_se  + 1, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_see + 1, face_fac * 1.0 * theta * w_ess[2]);

            set_value(values, col_b  + 1, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_e  + 1, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ee + 1, face_fac * 3.0 * theta * w_ess[2]);

            // face 1
            face_fac = 0.5 * dy * 0.25;
            set_value(values, col_b  + 1, face_fac * 3.0 * theta * w_ess[0]);
            set_value(values, col_e  + 1, face_fac * 3.0 * theta * w_ess[1]);
            set_value(values, col_ee + 1, face_fac * 3.0 * theta * w_ess[2]);

            set_value(values, col_nb  + 1, face_fac * 1.0 * theta * w_ess[0]);
            set_value(values, col_ne  + 1, face_fac * 1.0 * theta * w_ess[1]);
            set_value(values, col_nee + 1, face_fac * 1.0 * theta * w_ess[2]);
            
            // Contribution Delta r
            set_value(values, col_b  + 2, 0.0);
            set_value(values, col_e  + 2, 0.0);
            set_value(values, col_ee + 2, 0.0);
            //
            // Right hand side
            face_fac = 0.5;
            double htheta_b   = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_2]);
            double htheta_ip1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_5]);
            double htheta_ip2 = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_8]);
            double ht_ip12 = face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2);
            htheta_b   = 0.25 * (3.0 * htheta[p_1] + 1.0 * htheta[p_0]);
            htheta_ip1 = 0.25 * (3.0 * htheta[p_4] + 1.0 * htheta[p_3]);
            htheta_ip2 = 0.25 * (3.0 * htheta[p_7] + 1.0 * htheta[p_7]);            
            ht_ip12 += face_fac * (w_ess[0] * htheta_b + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2);

            double qtheta_b   = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_2]);
            double qtheta_ip1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_5]);
            double qtheta_ip2 = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_8]);
            double qt_ip12 = face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2);
            qtheta_b   = 0.25 * (3.0 * qtheta[p_1] + 1.0 * qtheta[p_0]);
            qtheta_ip1 = 0.25 * (3.0 * qtheta[p_4] + 1.0 * qtheta[p_3]);
            qtheta_ip2 = 0.25 * (3.0 * qtheta[p_7] + 1.0 * qtheta[p_6]);
            qt_ip12 += face_fac * (w_ess[0] * qtheta_b + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2);

            rhs[row] = -dy * ( qt_ip12 + c_wave * (ht_ip12 - h_given) );
            if (bc_vars[BC_WEST] == "zeta")
            {
                rhs[row] += dy * 2. * c_wave * bc[BC_WEST];
            }
            if (bc_vars[BC_WEST] == "q")
            {
                rhs[row] += dy * 2. * bc[BC_WEST];
            }
        }
        else if (bc_type[BC_WEST] == "borsboom")
        {
            //
            // Essential boundary condition
            // first row
            // ----------------------------
            // momentum + c_wave * continuity
            int col_b  = c_eq + 1 * 3;
            int col_e  = c_eq + 4 * 3;
            int col_ee = c_eq + 7 * 3;
            //
            if (do_convection) { con_fac = c_wave + qp_ip12 / hp_ip12; }
            //
            set_value(values, col_b , dtinv * con_fac * w_ess[0]);
            set_value(values, col_e , dtinv * con_fac * w_ess[1]);
            set_value(values, col_ee, dtinv * con_fac * w_ess[2]);
            // Contribution Delta q
            set_value(values, col_b  + 1, dtinv * w_ess[0]);
            set_value(values, col_e  + 1, dtinv * w_ess[1]);
            set_value(values, col_ee + 1, dtinv * w_ess[2]);
            // No contribution for Delta r
            set_value(values, col_b  + 2, 0.0);
            set_value(values, col_e  + 2, 0.0);
            set_value(values, col_ee + 2, 0.0);
            //
            double dhdt = dtinv * (hp_b - hn_b) * w_ess[0]
                + dtinv * (hp_ip1 - hn_ip1) * w_ess[1]
                + dtinv * (hp_ip2 - hn_ip2) * w_ess[2];
            double dqdt = dtinv * (qp_b - qn_b) * w_ess[0]
                + dtinv * (qp_ip1 - qn_ip1) * w_ess[1]
                + dtinv * (qp_ip2 - qn_ip2) * w_ess[2];
            rhs[row] = -(dqdt + con_fac * dhdt);

            double corr_term = 0.0;
            if (bc_vars[BC_WEST] == "zeta")
            {
                set_value(values, col_b , + dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_e , + dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ee, + dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double htheta_ip12 = w_ess[0] * htheta_b + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2;
                corr_term = + dhdt + ( eps_bc_corr * ((bc[BC_WEST] - zb_ip12) - htheta_ip12) );
                rhs[row] += corr_term;
            }
            if (bc_vars[BC_WEST] == "q")
            {
                set_value(values, col_b  + 1, dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0]);
                set_value(values, col_e  + 1, dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1]);
                set_value(values, col_ee + 1, dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2]);

                double qtheta_ip12 = w_ess[0] * qtheta_b + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2;
                corr_term = - dqdt + eps_bc_corr * (bc[BC_WEST] - qtheta_ip12);
                rhs[row] += corr_term;
            }
        }
        //------------------------------------------------------------------
        // natural boundary condition
        // second row
        // -----------------------------------------------------------------
        // momentum - c_wave * continuity
        //
        int col_b  = q_eq + 1 * 3;
        int col_e  = q_eq + 4 * 3;
        int col_ee = q_eq + 7 * 3;
        //
        if (do_convection) { con_fac = c_wave - qp_ip12 / hp_ip12; }
        double dhdt = dtinv * (hp_b - hn_b) * w_nat[0]
            + dtinv * (hp_ip1 - hn_ip1) * w_nat[1]
            + dtinv * (hp_ip2 - hn_ip2) * w_nat[2];
        double dqdx = dxinv * (qtheta_ip1 - qtheta_b);
        double dqdt = dtinv * (qp_b - qn_b) * w_nat[0]
            + dtinv * (qp_ip1 - qn_ip1) * w_nat[1]
            + dtinv * (qp_ip2 - qn_ip2) * w_nat[2];
        double dzetadx = dxinv * (htheta_ip1 + zb[p_4] - htheta_b - zb[p_1]);
        //
        // momentum part dq/dt + gh d(zeta)/dx 
        //
        // Contribution Delta h
        values[col_b     ] = w_nat[0] * theta * g * dzetadx - dxinv * theta * g * htheta_ip12;
        values[col_e     ] = w_nat[1] * theta * g * dzetadx + dxinv * theta * g * htheta_ip12;
        values[col_ee    ] = w_nat[2] * theta * g * dzetadx;
        // Contribution Delta q
        values[col_b  + 1] = dtinv * w_nat[0];
        values[col_e  + 1] = dtinv * w_nat[1];
        values[col_ee + 1] = dtinv * w_nat[2];
        // Contribution Delta r
        values[col_b  + 2] = 0.0;
        values[col_e  + 2] = 0.0;
        values[col_ee + 2] = 0.0;
        //
        rhs[row + 1] = -(dqdt + g * htheta_ip12 * dzetadx);
        //
        if (do_convection)
        {
        }
        if (do_bed_shear_stress)
        {
        }
        if (do_viscosity)
        {
        }
        //
        // continuity part -c_wave * (dhdt + dq/dx + dr/dy)
        //
        // Contribution Delta h
        values[col_b     ] += -con_fac * dtinv * w_nat[0];
        values[col_e     ] += -con_fac * dtinv * w_nat[1];
        values[col_ee    ] += -con_fac * dtinv * w_nat[2];
        // Contribution Delta q
        values[col_b  + 1] += -con_fac * -dxinv * theta;
        values[col_e  + 1] += -con_fac * dxinv * theta;
        values[col_ee + 1] += 0.0;
        // Contribution Delta r
        values[col_b  + 2] += 0.0;
        values[col_e  + 2] += 0.0;
        values[col_ee + 2] += 0.0;
                            
        rhs[row + 1] += - ( -con_fac * (dhdt + dqdx));
        //------------------------------------------------------------------
        // third row 
        //------------------------------------------------------------------
        col_b  = r_eq + 1 * 3;
        col_e  = r_eq + 4 * 3;
        col_ee = r_eq + 7 * 3;
        // Contribution Delta h
        values[col_b     ] = 0.0;
        values[col_e     ] = 0.0;
        values[col_ee    ] = 0.0;
        // Contribution Delta q
        values[col_b  + 1] = 0.0;
        values[col_e  + 1] = 0.0;
        values[col_ee + 1] = 0.0;
        // Contribution Delta r
        values[col_b  + 2] = 1.0;
        values[col_e  + 2] = 0.0;
        values[col_ee + 2] = 0.0;
        //
        rhs[row + 2] = 0.0;
    }
    return 0;
}
void  molecule(std::vector<int>& p, int p_sw, int ny)
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
inline void set_value(double * values, int col, double data)
{ 
    values[col] += data; 
}
inline int ma_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}

