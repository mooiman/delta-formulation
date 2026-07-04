//---------------------------------------------------------------
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
//    Solving the 1D advection/diffusion equation, fully implicit with delta-formuation and Modified Newton iteration 
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
//
// Advection equation: dc/dt + u dc/dx - eps d^2 c/dx^2 = 0
//

//select 
  //  1 uniform u>0
  //  3 eq. 5.89: Wesseling 2001 (Principles of Computational Fluid Dynamics)
  //  4 eq. 5.99 : Wesseling 2001 (Principles of Computational Fluid Dynamics)
  //  5 ex. 2 : Ten Thye Boonkkamp and Anthonissen(nonstationary)
  //  
#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

// for bicgstab  solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <toml.h>

void GetArguments(long argc, char** argv, std::string* file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

#include "cfts.h"
#include "ugrid1d.h"
#include "adv_diff_source.h"
#include "adv_diff_boundary_condition.h"
#include "adv_diff_init_concentration.h"
#include "adv_diff_init_velocity.h"
#include "adv_diff_linear_operator.h"
#include "adv_diff_regularization.h"
#include "perf_timer.h"
#include "definition_map_file.h"


int main(int argc, char* argv[])
{
    bool stationary = false;
    std::string toml_file_name("--- not-defined ---");
    int status = -1;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path current_dir;
    std::filesystem::path output_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();

    REGULARIZATION* regular = new REGULARIZATION();

    toml::table tbl;
    if (argc == 3 )
    {
        (void)GetArguments(argc, argv, &toml_file_name);
    }
    else
    {
        toml_file_name = "../../test_cases/advection_diffusion/001/input_cpp.toml";
        // toml_file_name = "../../test_cases/advection_diffusion/002/input_cpp.toml";
        // toml_file_name = "../../test_cases/advection_diffusion/003/input_cpp.toml";
        // toml_file_name = "../../test_cases/advection_diffusion/004/input_cpp.toml";
        // 
        toml_file_name = "../../test_cases/advection_diffusion/006/input_cpp.toml";
        // 
        // toml_file_name = "../../test_cases/diffusion_only/001/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_only/002/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_only/003/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_only/004/input_cpp.toml";
        //
        // toml_file_name = "../../test_cases/diffusion_discontinue/001/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_discontinue/002/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_discontinue/003/input_cpp.toml";
        // toml_file_name = "../../test_cases/diffusion_discontinue/004/input_cpp.toml";
        // 
    }
    if (!std::filesystem::exists(toml_file_name))
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Input file: \n" << std::filesystem::absolute(toml_file_name) << std::endl;
        std::cout << "can not be opened." << std::endl;
        std::cout << "----------------------------" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
        exit(1);
    }
    tbl = toml::parse_file(toml_file_name);
    std::filesystem::path file_toml;
    file_toml = toml_file_name;
    current_dir = file_toml.parent_path();
    output_dir = current_dir.string() + "/output/";

    bool use_eq8 = tbl["Numerics"]["use_eq8"].value_or(bool(true));
    if (!use_eq8)
    {
        output_dir = current_dir.string() + "/output_eq8_not_used/";
    }
    std::filesystem::create_directory(output_dir);

    START_TIMERN(main);
    START_TIMER(Simulation initialization);

    std::string out_file;
    std::stringstream ss;
    ss << "adv_diff";
    out_file = output_dir.string() + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("SMART Numerics: advection-diffusion, BiCGSTab");

    std::ofstream log_file;
    log_file.open(log_filename);
    std::cout << "=== Input file =======================================" << std::endl;
    std::cout << std::filesystem::absolute(toml_file_name) << std::endl;
    std::cout << "======================================================" << std::endl;
    log_file << "======================================================" << std::endl;
    log_file << "Compile date and time executable: " << __DATE__ << ", " << __TIME__ << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;
    log_file << "=== Copy of the input file ============================" << std::endl;
    log_file << tbl << "\n";  // Write input TOML file to log_file
    log_file << "=======================================================" << std::endl;
    log_file << std::endl;
    log_file << "=== Used input variables ==============================" << std::endl;

    // Domain
    log_file << "[Domain]" << std::endl;
    double Lx = tbl["Domain"]["Lx"].value_or(double(1.0));
    log_file << "Lx = " << Lx << std::endl;
    double Interface = tbl["Domain"]["Interface"].value_or(double(0.5));  // Location of interface[-], fraction of Lx
    log_file << "Interface = " << Interface << std::endl;

    // Time
    log_file << std::endl << "[Time]" << std::endl;
    double tstart = tbl["Time"]["tstart"].value_or(double(0.0));
    log_file << "tstart = " << tstart << std::endl;
    double tstop = tbl["Time"]["tstop"].value_or(double(2.0 * 3600.));
    log_file << "tstop = " << tstop << std::endl;
    double dt = tbl["Time"]["dt"].value_or(double(0.0));  // default stationary
    if (dt == 0.0) { stationary = true;  }
    log_file << "dt = " << dt << std::endl;

    // boundary conditions
    log_file << std::endl << "[Boundary]" << std::endl;
    //std::vector<double> bc_vals{ 0., 0. };
    //status = get_toml_array(tbl, "bc_vals", bc_vals);
    double treg = tbl["Boundary"]["treg"].value_or(double(300.0));
    if (stationary)
    {
        log_file << "treg = " << 0.0 << std::endl;
    }
    else
    {
        log_file << "treg = " << treg << std::endl;
    }

    //Physics
    log_file << std::endl << "[Physics]" << std::endl;
    bool discon_diff = false;
    double eps_left;
    double eps_right;
    double u_const = tbl["Physics"]["u_const"].value_or(double(0.0));  // default, no velocitty
    log_file << "u_const = " << u_const << std::endl;
    double eps_const = tbl["Physics"]["eps_const"].value_or(double(0.0));
    log_file << "eps_const = " << eps_const << std::endl;
    if (eps_const == 0.0)
    {
        eps_left = tbl["Physics"]["eps_left"].value_or(double(0.0));
        log_file << "eps_left = " << eps_left << std::endl;
        eps_right = tbl["Physics"]["eps_right"].value_or(double(0.0));
        log_file << "eps_right = " << eps_right << std::endl;
        discon_diff = true;
    }

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    double dx = tbl["Numerics"]["dx"].value_or(double(0.01));
    log_file << "dx = " << dx << std::endl;
    double c_psi = tbl["Numerics"]["c_psi"].value_or(double(4.0));
    log_file << "c_psi = " << c_psi << std::endl;
    double theta = tbl["Numerics"]["theta"].value_or(double(0.501));
    if (stationary)
    {
        log_file << "theta = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "theta = " << theta << std::endl;
    }
    double alpha = tbl["Numerics"]["alpha"].value_or(double(0.125));  // 0.125 (Mass matrix)
    log_file << "alpha = " << alpha << std::endl;
    int iter_max = tbl["Numerics"]["iter_max"].value_or(int(150));
    log_file << "iter_max = " << iter_max << std::endl;
    double eps_newton = tbl["Numerics"]["eps_newton"].value_or(double(1.0e-12));
    log_file << "eps_newton = " << eps_newton << std::endl;
    double eps_bicgstab = tbl["Numerics"]["eps_bicgstab"].value_or(double(1.0e-12));
    log_file << "eps_bicgstab = " << eps_bicgstab << std::endl;
    //bool use_eq8 = tbl["Numerics"]["use_eq8"].value_or(bool(true));
    log_file << "use_eq8  = " << use_eq8 << std::endl;
    bool regularization_init = tbl["Numerics"]["regularization_init"].value_or(bool(true));
    log_file << "regularization_int  = " << regularization_init << std::endl;
    bool regularization_iter = tbl["Numerics"]["regularization_iter"].value_or(bool(true));
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    bool regularization = tbl["Numerics"]["regularization"].value_or(bool(true));
    log_file << "regularization = " << regularization << std::endl;

    //output
    log_file << std::endl << "[Output]" << std::endl;
    double dt_his = tbl["Output"]["dt_his"].value_or(double(1.0));
    if (stationary)
    {
        log_file << "dt_his = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_his = " << dt_his << std::endl;
    }
    double dt_map = tbl["Output"]["dt_map"].value_or(double(10.0));
    if (stationary)
    {
        log_file << "dt_map = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_map = " << dt_map << std::endl;
    }

    int select = 0;
    select = 2;  // virtual points random value (ie for boundary at right side c[0] = 1 and c[Right] = 0
    select = 1;  // uniform u>0
    // select = 3;  // eq. 5.89: Wesseling 2000 (Principles of Computational Fluid Dynamics)
    // select = 4;  // eq. 5.99: Wesseling 2000 (Principles of Computational Fluid Dynamics)
    // select = 5;  // ex. 2   : Ten Thye Boonkkamp and Anthonissen (nonstationary)

    int total_time_steps = int((tstop - tstart) / dt) + 1;  // [sec]
    double dxinv = 1. / dx;
    int nr_nodes = int(Lx / dx) + 1 + 2; // nr nodes; including 2 virtual points

    double dtinv;
    double dtinv_pseu;
    int wrimap;
    int wrihis;
    if (stationary)
    {
        dtinv = 0.0;                                      // stationary solution
        dtinv_pseu = 0.0;  //  educational guess: equal to velocity
        theta = 1.0;                                      // Stationary solution
        tstop = 1.;
        total_time_steps = 2;  // [sec]
        treg = 0.0;                                      // Thatcher-Harleman return time [s], when zero supply boundary value immediately
        wrihis = 1;      // write interval to his-file (every delta t)
        wrimap = 1;     // write interval to map-file (every 1 sec , or every delta t)
    }
    else
    {
        dtinv = 1. / dt;                                  // Inverse of dt [1/s]
        dtinv_pseu = 0.0;
        wrihis = std::max(int(dt * dtinv), int(dt_his * dtinv));      // write interval to his-file (every delta t)
        wrimap = std::max(int(dt * dtinv), int(dt_map * dtinv));     // write interval to map-file (every 1 sec , or every delta t)
    }

    double bc0;
    double bc1;

    std::vector<double> x(nr_nodes, 0.);  // x-coordinate
    std::vector<double> y(nr_nodes, 0.);  // y-coordinate
    std::vector<double> cn(nr_nodes, 0.0);  // constituent [-]
    std::vector<double> cp(nr_nodes, 0.);  // constituent [-]
    std::vector<double> q(nr_nodes, 0.);  // source
    std::vector<double> u(nr_nodes, 0.);  // advection velocity [m s-1], adv_diff_init_velocity
    std::vector<double> eps(nr_nodes, eps_const);  // diffusion [m2 s-1]
    std::vector<double> eps_given(nr_nodes, eps_const);  // constant in time, diffusion [m2 s-1]
    std::vector<double> eps_reg(nr_nodes, eps_const);  // regularized once at the beginning of the computation [m2 s-1]
    std::vector<double> pe(nr_nodes, 0.);  // peclet number [-]
    std::vector<double> psi(nr_nodes, 0.);  //
    std::vector<double> eq8(nr_nodes, 0.);  //
    std::vector<double> delta_c(nr_nodes, 0.);
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix
    std::vector<double> w_bc(3, 0.);  // weighting coefficients of the mass-matrix
    std::vector<double> tmp(nr_nodes, 0.);  // dummy vector

    Eigen::VectorXd solution(nr_nodes);    // solution vector [c]^{n}
    Eigen::VectorXd rhs(nr_nodes);        // RHS vector [c]^n

    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    double alpha_bc = (2. * alpha - 1. / 2.);
    alpha_bc = 0.0;
    w_bc[0] = 1. + alpha_bc;  // 0.75
    w_bc[1] = 1. - 2. * alpha_bc; // 1.5
    w_bc[2] = alpha_bc;  // -0.25

    //initialize x-coordinate
    for (int i = 0; i < nr_nodes; i++)
    {
        x[i] = double(i-1) * dx;
    }

    //initialize diffusion
    for (int i = 0; i < nr_nodes; ++i)
    {
        eps[i] = eps_const;
        eps_reg[i] = eps_const;
    }
    if (discon_diff)
    {
        double x_star = Interface * Lx;
        // double beta = 1.0 / (x_star - eps_left * x_star + eps_left);
        for (int i = 0; i < nr_nodes; ++i)
        {
            if (x[i] <= x_star)
            {
                // c_ana[i] = beta * x[i];
                eps_given[i] = eps_left;
            }
            else
            {
                // c_ana[i] = eps_left * beta * x[i] + eps_right - eps_left * beta;
                eps_given[i] = eps_right;
            }
            eps[i] = eps_given[i];
            eps_reg[i] = eps_given[i];
        }
    }

    // initial concentration
    (void)adv_diff_init_concentration(cn, select);
    (void)adv_diff_init_velocity(u, u_const, x, select);
    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        (void)regular->given_function(eps_reg, psi, eq8, eps_given, u_const, dx, c_psi, use_eq8);
        STOP_TIMER(Regularization_init);
        log_file << "=== eps_reg ===========================================" << std::endl;
        for (int i = 0; i < nr_nodes; ++i)
        {
            log_file << std::setprecision(15) << std::scientific << eps_reg[i] << std::endl;
        }
    }
    for (int i = 0; i < nr_nodes; ++i)
    {
        eps[i] = eps_reg[i];
        pe[i] = u[i] * dx / eps[i];
    }

    ////////////////////////////////////////////////////////////////////////////
    // Define map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    std::string map_c_name("cn_1d");
    std::string map_u_name("u_1d");
    std::string map_eps_name("eps_1d");
    std::string map_psi_name("psi_1d");
    std::string map_eq8_name("eq8_1d");
    std::string map_pe_name("pe_1d");
    UGRID1D* map_file = create_map_file(nc_mapfilename, x, map_c_name, map_u_name, map_eps_name, map_psi_name, map_eq8_name, map_pe_name);
    // Put data on map file
    int nst_map = 0;
    map_file->put_time(nst_map, double(0)* dt);
    map_file->put_time_variable(map_c_name, nst_map, cn);
    map_file->put_time_variable(map_u_name, nst_map, u);
    map_file->put_time_variable(map_eps_name, nst_map, eps);
    map_file->put_time_variable(map_psi_name, nst_map, psi);
    map_file->put_time_variable(map_eq8_name, nst_map, eq8);
    map_file->put_time_variable(map_pe_name, nst_map, pe);
    // End definition of map file
    ////////////////////////////////////////////////////////////////////////////
    // Define time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations (corrected for virtual points)
    int i_left = 1;
    int i_mid_left = nr_nodes / 4 + 1;
    int i_mid = nr_nodes / 2;
    int i_mid_right = 3 * nr_nodes / 4 - 1;
    int i_right = nr_nodes - 2;
    std::vector<double> x_obs = { x[i_left], x[i_mid_left], x[i_mid], x[i_mid_right], x[i_right] };
    std::vector<double> y_obs = { y[i_left], y[i_mid_left], y[i_mid], y[i_mid_right], y[i_right] };

    std::vector<std::string> obs_stations;
    obs_stations.push_back("West boundary");
    obs_stations.push_back("Halfway to west boundary");
    obs_stations.push_back("Centre");
    obs_stations.push_back("Halfway to east boundary");
    obs_stations.push_back("East boundary");
    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_cn_name("constituent");
    his_file->add_variable(his_cn_name, "", "Constituent", "-");

    std::string his_outer_iter_name("outer_iterations");
    his_file->add_variable_without_location(his_outer_iter_name, "", "Newton iterations", "-");
    std::string his_lin_solv_iter_name("lin_solver_iterations");
    his_file->add_variable_without_location(his_lin_solv_iter_name, "", "Lin. solver iterations", "-");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, double(0)* dt);

    std::vector<double> his_values = { cn[i_left], cn[i_mid_left], cn[i_mid], cn[i_mid_right], cn[i_right] };
    his_file->put_variable(his_cn_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_outer_iter_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);

    ////////////////////////////////////////////////////////////////////////////

    // source term q
    (void)adv_diff_source(q, select);

    double time = tstart + dt * double(0);
    //(void)adv_diff_boundary_condition(cn[0], cn[nr_nodes - 1], time, treg, select);
    // compute Peclet number
    for (int i = 0; i < nr_nodes; i++)
    {
        cp[i] = cn[i];  // place the old solution in the new one, needed for the iteration procedure
    }

    Eigen::SparseMatrix<double> A(nr_nodes, nr_nodes);
    for (int i = 0; i < nr_nodes; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
        solution[i] = 0.0;
    }

    // start time loop
    std::cout << "Start time-loop" << std::endl;
    std::cout << time << ";   " << dt << ";   " << tstart + tstop << std::endl;

    double dc_max = 0.0;
    int dc_maxi = 0;
    STOP_TIMER(Simulation initialization);
    START_TIMER(Time loop);
    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        int used_iter = 0;
        int used_lin_iter = 0;
        time = dt * double(nst);
#if defined(DEBUG)
        log_file << "=== Start time step ===================================" << std::endl;
        log_file << std::setprecision(4) << time << std::endl;
#endif
        adv_diff_boundary_condition(bc0, bc1, time, treg, select);
        for (int iter = 0; iter < iter_max; ++iter)
        {
#if defined(DEBUG)   
            //log_file << "====== Start iteration ================================" << std::endl;
#endif
            if (nst == 1 && iter == 0)
            {
                START_TIMER(Matrix initialization);
            }
            START_TIMER(Matrix set up);
            // interior nodes
            for (int i = 1; i < nr_nodes - 1; ++i)
            {
                double cn_im1 = cn[i - 1]; // = h^{n}_{i-1}
                double cn_i = cn[i];       // = h^{n}_{i}
                double cn_ip1 = cn[i + 1]; // = h^{n}_{i+1}
                double cn_im12 = 0.5 * (cn[i - 1] + cn[i]);
                double cn_ip12 = 0.5 * (cn[i] + cn[i + 1]);

                double cp_im1 = cp[i - 1];       // = h^{n+1,p}_{i-1}
                double cp_i = cp[i];       // = h^{n+1,p}_{i}
                double cp_ip1 = cp[i + 1];       // = h^{n+1,p}_{i-2}
                double cp_im12 = 0.5 * (cp[i - 1] + cp[i]);
                double cp_ip12 = 0.5 * (cp[i] + cp[i + 1]);

                double ctheta_i = theta * cp_i + (1.0 - theta) * cn_i;
                double ctheta_im1 = theta * cp_im1 + (1.0 - theta) * cn_im1;
                double ctheta_ip1 = theta * cp_ip1 + (1.0 - theta) * cn_ip1;
                double ctheta_im12 = theta * cp_im12 + (1.0 - theta) * cn_im12;
                double ctheta_ip12 = theta * cp_ip12 + (1.0 - theta) * cn_ip12;

                double u_im12 = 0.5 * (u[i - 1] + u[i]);
                double u_ip12 = 0.5 * (u[i] + u[i + 1]);

                double eps_im12 = 0.5 * (eps[i - 1] + eps[i]);
                double eps_ip12 = 0.5 * (eps[i] + eps[i + 1]);

                A.coeffRef(i, i - 1) = dx * dtinv * mass[0]
                    - 0.5 * u_im12 * theta
                    - eps_im12 * theta * dxinv;
                A.coeffRef(i, i) = dtinv_pseu + dx * dtinv * mass[1]
                    - 0.5 * u_im12 * theta + 0.5 * u_ip12 * theta
                    + eps_im12 * theta * dxinv
                    + eps_ip12 * theta * dxinv;
                A.coeffRef(i, i + 1) = dx * dtinv * mass[2]
                    + 0.5 * u_ip12 * theta
                    - eps_ip12 * theta * dxinv;
                rhs[i] = -(
                    dx * dtinv * mass[0] * (cp[i - 1] - cn[i - 1])
                    + dx * dtinv * mass[1] * (cp[i] - cn[i])
                    + dx * dtinv * mass[2] * (cp[i + 1] - cn[i + 1])
                    + u_ip12 * cp_ip12 - u_im12 * cp_im12
                    + eps_im12 * dxinv * (ctheta_i - ctheta_im1) - eps_ip12 * dxinv * (ctheta_ip1 - ctheta_i)
                    );
            }
            {
                //
                // west boundary
                //
                int i = 0;
                double cn_i = cn[i];       // = c^{n}_{i}
                double cn_ip1 = cn[i + 1];       // = c^{n}_{i+1}
                double cn_ip2 = cn[i + 2];       // = c^{n}_{i+2}
                double cn_ip12 = 0.5 * (cn_i + cn_ip1) + 0.5 * alpha_bc * (cn_i - 2. * cn_ip1 + cn_ip2);
                double cp_i = cp[i];       // = c^{n+1,p}_{i}
                double cp_ip1 = cp[i + 1];       // = c^{n+1,p}_{i+1}
                double cp_ip2 = cp[i + 2];       // = c^{n+1,p}_{i+2}
                double cp_ip12 = 0.5 * (cp_i + cp_ip1) + 0.5 * alpha_bc * (cp_i - 2. * cp_ip1 + cp_ip2);
                double ctheta_i = theta * cp_i + (1.0 - theta) * cn_i;
                double ctheta_ip1 = theta * cp_ip1 + (1.0 - theta) * cn_ip1;
                double ctheta_ip2 = theta * cp_ip2 + (1.0 - theta) * cn_ip2;
                double ctheta_ip12 = 0.5 * (ctheta_i + ctheta_ip1) + 0.5 * alpha_bc * (ctheta_i - 2.0 * ctheta_ip1 + ctheta_ip2);

                double eps_ip12 = 0.5 * (eps[i] + eps[i + 1]);

                if (false)
                {
                    // Homogeneous Neumann
                    A.coeffRef(i, i) = 1.0;
                    A.coeffRef(i, i + 1) = -1.0;
                    A.coeffRef(i, i + 2) = 0.0;
                    rhs[i] = 0.0;
                }
                else
                {
                    // Dirichlet
                    A.coeffRef(i, i) = 1.0 / 12.0 + 0.5 * theta * w_bc[0];
                    A.coeffRef(i, i + 1) = 10.0 / 12.0 + 0.5 * theta * w_bc[1];
                    A.coeffRef(i, i + 2) = 1.0 / 12.0 + 0.5 * theta * w_bc[2];
                    rhs[i] = +bc0 - cp_ip1 - 0.5 * (cp_i - 2. * cp_ip1 + cp_ip2);

                    A.coeffRef(i, i) = 1.0 / 12.0;
                    A.coeffRef(i, i + 1) = 10.0 / 12.0;
                    A.coeffRef(i, i + 2) = 1.0 / 12.0;
                    rhs[i] = +bc0 - cp_ip1;  // -0.5 * (cp_i - 2. * cp_ip1 + cp_ip2);

                    A.coeffRef(i, i) = 0.0;
                    A.coeffRef(i, i + 1) = 1.0;
                    A.coeffRef(i, i + 2) = 0.0;
                    rhs[i] = +bc0 - cp_ip1;

                    A.coeffRef(i + 1, i) = -1.0;
                    A.coeffRef(i + 1, i + 1) = 1.0;
                    A.coeffRef(i + 1, i + 2) = 0.0;
                    rhs[i + 1] = 0.0;

                }
            }
            {
                //
                // east boundary
                //
                int i = nr_nodes - 1;
                double cn_i = cn[i];       // = h^{n}_{i}
                double cn_im1 = cn[i - 1];       // = h^{n}_{i-1}
                double cn_im2 = cn[i - 2];       // = h^{n}_{i-2}
                double cn_im12 = 0.5 * (cn_i + cn_im1) + 0.5 * alpha_bc * (cn_im2 - 2. * cn_im1 + cn_i);
                double cp_i = cp[i];       // = h^{n+1,p}_{i}
                double cp_im1 = cp[i - 1];       // = h^{n+1,p}_{i-1}
                double cp_im2 = cp[i - 2];       // = h^{n+1,p}_{i-2}
                double cp_im12 = 0.5 * (cp_i + cp_im1) + 0.5 * alpha_bc * (cp_im2 - 2. * cp_im1 + cp_i);
                double ctheta_i = theta * cp_i + (1.0 - theta) * cn_i;
                double ctheta_im1 = theta * cp_im1 + (1.0 - theta) * cn_im1;
                double ctheta_im2 = theta * cp_im2 + (1.0 - theta) * cn_im2;
                double ctheta_im12 = 0.5 * (ctheta_i + ctheta_im1) + 0.5 * alpha_bc * (ctheta_im2 - 2.0 * ctheta_im1 + ctheta_i);

                double eps_im12 = 0.5 * (eps[i - 1] + eps[i]);

                if (false)
                {
                    // Homogeneous Neumann
                    A.coeffRef(i, i) = 1.0;
                    A.coeffRef(i, i - 1) = -1.0;
                    A.coeffRef(i, i - 2) = 0.0;
                    rhs[i] = 0.0;
                }
                else
                {
                    // Dirichlet
                    A.coeffRef(i, i) = 1.0 / 12.0 + 0.5 * theta * w_bc[0];
                    A.coeffRef(i, i - 1) = 10.0 / 12.0 + 0.5 * theta * w_bc[1];
                    A.coeffRef(i, i - 2) = 1.0 / 12.0 + 0.5 * theta * w_bc[2];
                    rhs[i] = bc1 - cp_im1 - 0.5 * (cp_i - 2. * cp_im1 + cp_im2);

                    A.coeffRef(i, i) = 1.0 / 12.0;
                    A.coeffRef(i, i - 1) = 10.0 / 12.0;
                    A.coeffRef(i, i - 2) = 1.0 / 12.0;
                    rhs[i] = bc1 - cp_im1;  // -0.5 * (cp_i - 2. * cp_im1 + cp_im2);

                    A.coeffRef(i, i) = 0.0;
                    A.coeffRef(i, i - 1) = 1.0;
                    A.coeffRef(i, i - 2) = 0.0;
                    rhs[i] = +bc1 - cp_im1;

                    A.coeffRef(i-1, i) = -1.0;
                    A.coeffRef(i-1, i - 1) = 1.0;
                    A.coeffRef(i-1, i - 2) = 0.0;
                    rhs[i-1] = 0.0;

                }
            }
            STOP_TIMER(Matrix set up);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(Matrix initialization);
            }
            if (nst == 1 && iter == 0)
            {
                START_TIMER(BiCGStab_initialization);
            }
            else
            {
                START_TIMER(BiCGStab);
            }

            Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
            solver.compute(A);
            //solution = solver.solve(rhs);
            solver.setTolerance(eps_bicgstab);
            solution = solver.solveWithGuess(rhs, solution);
            START_TIMER(Set solution);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(BiCGStab_initialization);
            }
            else
            {
                STOP_TIMER(BiCGStab);
            }
#if defined(DEBUG)
            //log_file << "time [sec]:" << dt * double(nst) 
            //         << "    iterations     :" << solver.iterations()
            //         << "    estimated error:" << solver.error() 
            //         << std::endl;
#endif

            dc_max = 0.0;
            dc_maxi = 0;
            for (int i = 0; i < nr_nodes; ++i)
            {
                cp[i] += solution[i];
                delta_c[i] = solution[i]; 
                if (dc_max < std::abs(delta_c[i]))
                {
                    dc_max = std::abs(delta_c[i]);
                    dc_maxi = i;
                }
            }

            if (nst <= -1)
            {
                log_file << "=== Matrix ============================================" << std::endl;
                log_file << std::setprecision(4) << std::scientific << Eigen::MatrixXd(A) << std::endl;
                log_file << "=== RHS ===============================================" << std::endl;
                log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
            }
            if (nst <= -1)
            {
                log_file << "=== cn, delta_c, cp ===================================" << std::endl;
                for (int i = 0; i < nr_nodes; ++i)
                {
                    log_file << std::setprecision(15) << std::scientific << cn[i] << "  -----  " << delta_c[i] << "  -----  " << cp[i]  << std::endl;
                }
                //log_file << "=== Eigen values ======================================" << std::endl;
                //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).eigenvalues() << std::endl;
                //log_file << "====== End iteration ==================================" << std::endl;
                //log_file << "time [sec]: " << dt * double(nst)
                //         << ";    iterations     : " << solver.iterations()
                //         << ";    estimated error: " << solver.error() 
                //         << std::endl;
            }

            used_iter = iter+1;
            used_lin_iter += (int) solver.iterations();
            STOP_TIMER(Set solution);
            if (dc_max < eps_newton || iter == iter_max - 1)
            {
                if (iter == iter_max - 1)
                {
                    log_file << "=======================================================" << std::endl;
                    log_file << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << ";   " << used_iter << std::endl;
                    log_file << std::scientific << std::setprecision(10) <<  "maximum delta c at index " << dc_maxi << " is " << dc_max << std::endl;
                }
                break;
            }
            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                (void) regular->given_function(cp, psi, eq8, cp, u_const, dx, c_psi, use_eq8);
                STOP_TIMER(Regularization_iter_loop);
                for (int i = 0; i < nr_nodes; ++i)
                {
                    cp[i] = tmp[i];
                    eps[i] = eps_reg[i];
                }
            }
        }
        std::cout << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << ";   " << used_iter << std::endl;
        if (nst <= -1)
        {
            log_file << "=======================================================" << std::endl;
            log_file << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << ";   " << used_iter << std::endl;
        }

        START_TIMER(Regularization_time_loop);
        (void)(void)regular->given_function(tmp, psi, eq8, cp, u_const, dx, c_psi, use_eq8);
        STOP_TIMER(Regularization_time_loop);
        if (regularization)
        {
            for (int i = 0; i < nr_nodes; ++i)
            {
                cn[i] = tmp[i];
                eps[i] = eps_reg[i] + std::abs(psi[i]);
                pe[i] = u[i] * dx / eps[i];
            }
        }
        else
        {
            for (int i = 0; i < nr_nodes; ++i)
            {
                cn[i] = cp[i];
            }
        }
        if (nst <= -1)
        {
            log_file << "=== cp, cn, cp-cn =====================================" << std::endl;
            for (int i = 0; i < nr_nodes; ++i)
            {
                log_file << std::setprecision(15) << std::scientific << cp[i] << ",   " << cn[i] << ",   " << cp[i] - cn[i] << std::endl;
            }
        }
        // Map-files
        if (std::fmod(nst, wrimap) == 0)
        {
            // Put data on time map file
            START_TIMER(Writing map - file);
            nst_map++;
            if (stationary)
            {
                map_file->put_time(nst_map, double(nst));
            }
            else
            {
                map_file->put_time(nst_map, double(nst)* dt);
            }
            map_file->put_time_variable(map_c_name, nst_map, cn);
            map_file->put_time_variable(map_u_name, nst_map, u);
            map_file->put_time_variable(map_eps_name, nst_map, eps);
            map_file->put_time_variable(map_psi_name, nst_map, psi);
            map_file->put_time_variable(map_eq8_name, nst_map, eq8);
            map_file->put_time_variable(map_pe_name, nst_map, pe);
            STOP_TIMER(Writing map - file);
        }
        // His-files
        if (std::fmod(nst, wrihis) == 0)
        {
            START_TIMER(Writing his - file);
            nst_his++;
            if (stationary)
            {
                his_file->put_time(nst_his, double(nst));
            }
            else
            {
                his_file->put_time(nst_his, double(nst) * dt);
            }
            std::vector<double> his_values = { cn[i_left], cn[i_mid_left], cn[i_mid], cn[i_mid_right],  cn[i_right] };
            his_file->put_variable(his_cn_name, nst_his, his_values);

            his_values = { double(used_iter) };
            his_file->put_variable(his_outer_iter_name, nst_his, his_values);
            his_values = { double(used_lin_iter) };
            his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);
            STOP_TIMER(Writing his - file);
        }
    } // End of the time loop
    STOP_TIMER(Time loop);
    log_file << "=======================================================" << std::endl;
    log_file.close();
    (void)map_file->close();
    (void)his_file->close();
    STOP_TIMER(Writing log - file);

    STOP_TIMER(main);
    PRINT_TIMER(timing_filename.data());
}
//------------------------------------------------------------------------------
/* @@ GetArguments
*
* Scan command-line arguments
*/
void GetArguments(long argc,   /* I Number of command line arguments */
    char** argv,
    std::string* file_name)
    /* Returns nothing */
{
    long  i;
    i = 0;
    while (i < argc)
    {
        if (strcmp(argv[i], "--toml") == 0)
        {
            i = i + 1;
            if (i < argc) *file_name = std::string(argv[i]);
        }
        else if (i != 0)
        {
            i = argc;
        }
        i = i + 1;
    }
    return;
}
int get_toml_array(toml::table tbl, std::string keyw, std::vector<std::string>& values)
{
    int status = 1;
    toml::array& arr = *tbl.get_as<toml::array>(keyw);
    if (&arr != nullptr)
    {
        values.clear();
        for (auto&& elem : arr)
        {
            auto e = elem.value_or("none");
            values.emplace_back(e);
        }
    }
    if (&arr != nullptr) status = 0;
    return status;
}
int get_toml_array(toml::table tbl, std::string keyw, std::vector<double>& values)
{
    int status = 1;
    toml::array& arr = *tbl.get_as<toml::array>(keyw);
    if (&arr != nullptr)
    {
        values.clear();
        for (auto&& elem : arr)
        {
            auto e = elem.value_or(0.0);
            values.emplace_back(e);
        }
    }
    if (&arr != nullptr) status = 0;
    return status;
}
int get_toml_array(toml::table tbl, std::string keyw, std::vector<bool>& values)
{
    int status = 1;
    toml::array& arr = *tbl.get_as<toml::array>(keyw);
    if (&arr != nullptr)
    {
        values.clear();
        for (auto&& elem : arr)
        {
            auto e = elem.value_or(false);
            values.emplace_back(e);
        }
    }
    if (&arr != nullptr) status = 0;
    return status;
}

