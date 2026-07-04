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
// Advection equation: dc/dt + d(uc)/dx = 0
//
#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <thread>

// for bicgstab  solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <toml.h>

void GetArguments(long argc, char** argv, std::string* file_name);
std::vector<double> bc_val(double x_loc);
std::vector<double> quadratic_interpolation_1d(std::vector<double> x_gp, std::vector<double> f_gp, int bc_east_west, double x_loc);
std::vector<double> hermite_interpolation_1d(std::vector<double> x_gp, std::vector<double> f_gp, int bc_east_west, double x_loc);
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
#include "regularization.h"
#include "perf_timer.h"
#include "definition_map_file.h"

int BC_WEST = 0;
int BC_EAST = 1;

double scv(double, double);
std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name);

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

    toml::table tbl;
    toml::table tbl_chp;  // table for a chapter
    if (argc == 3)
    {
        (void)GetArguments(argc, argv, &toml_file_name);
        if (!std::filesystem::exists(toml_file_name))
        {
            std::cout << "----------------------------" << std::endl;
            std::cout << "Input file \'" << toml_file_name << "\' can not be opened." << std::endl;
            std::cout << "Press Enter to finish";
            std::cin.ignore();
            exit(1);
        }
        tbl = toml::parse_file(toml_file_name);
        std::filesystem::path file_toml;
        file_toml = toml_file_name;
        current_dir = file_toml.parent_path();
        output_dir = current_dir.string() + "/output/";
        std::filesystem::create_directory(output_dir);
    }
    else
    {
        std::cout << "No \'toml\' file is read." << std::endl;
        current_dir = ".";
        output_dir = ".";
    }


    START_TIMERN(main);
    START_TIMER(Simulation initialization);
    SHAPE_CONC shape_conc = SHAPE_CONC::NONE;
    BND_TYPE bnd_type = BND_TYPE::NONE;
    int select_src = 0;  // 0: constant source == 0

    std::string out_file;
    std::stringstream ss;

    std::string logging = tbl["Logging"].value_or("None");

    ss << "advection_phi";
    out_file = output_dir.string() + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Advection exp(phi), 1D, BiCGSTab");

    std::ofstream log_file;
    log_file.open(log_filename);
    std::cout << "=== Input file =======================================" << std::endl;
    std::cout << "Executable compiled: " << __DATE__ << ", " << __TIME__ << std::endl;
    std::cout << std::filesystem::absolute(toml_file_name) << std::endl;
    std::cout << "======================================================" << std::endl;
    log_file << "======================================================" << std::endl;
    log_file << "Executable compiled: " << __DATE__ << ", " << __TIME__ << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;
    log_file << "=== Copy of the input file ============================" << std::endl;
    log_file << tbl << "\n";  // Write input TOML file to log_file
    log_file << "=======================================================" << std::endl;
    log_file << std::endl;

    // Domain
    tbl_chp = *tbl["Domain"].as_table();
    double Lx = tbl_chp["Lx"].value_or(double(1.0));

    // Time
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    double tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.));

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    std::string shape_type = tbl_chp["shape"].value_or("None");
    if (shape_type == "Constant") { shape_conc = SHAPE_CONC::Constant; }
    else if (shape_type == "Envelope") { shape_conc = SHAPE_CONC::Envelope; }
    else if (shape_type == "EnvelopePhi") { shape_conc = SHAPE_CONC::EnvelopePhi; }
    else { shape_conc = SHAPE_CONC::NONE; }

    // Boundary
    tbl_chp = *tbl["Boundary"].as_table();
    std::vector<std::string> bc_type;
    status = get_toml_array(tbl_chp, "bc_type", bc_type);
    if (bc_type[0] == "Constant") { bnd_type = BND_TYPE::Constant; }
    else if (bc_type[0] == "Sine") { bnd_type = BND_TYPE::Sine; }
    else { bnd_type = BND_TYPE::NONE; }
    std::vector<double> bc_vals;
    status = get_toml_array(tbl_chp, "bc_vals", bc_vals);
    double treg = tbl_chp["treg"].value_or(double(300.0));

    //Physics
    tbl_chp = *tbl["Physics"].as_table();
    double g = tbl_chp["g"].value_or(double(10.));  // Gravitational acceleration
    double u_const = tbl_chp["u_const"].value_or(double(0.0));  // default, no velocitty
    bool momentum_viscosity = tbl_chp["momentum_viscosity"].value_or(bool(false));
    double visc_const = tbl_chp["viscosity"].value_or(double(0.0));

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    double dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary
    if (dt == 0.0) { stationary = true; }
    double dx = tbl_chp["dx"].value_or(double(0.01));
    double c_psi = tbl_chp["c_psi"].value_or(double(4.0));
    double theta = tbl_chp["theta"].value_or(double(0.501));
    double alpha = tbl_chp["alpha"].value_or(double(0.125));  // 0.125 (Mass matrix)
    int iter_max = tbl_chp["iter_max"].value_or(int(150));
    double eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));
    double eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));
    bool regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    bool regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    bool regularization_time = tbl_chp["regularization_time"].value_or(bool(false));

    //Output
    log_file << std::endl << "[Output]" << std::endl;
    tbl_chp = *tbl["Output"].as_table();
    double dt_his = tbl_chp["dt_his"].value_or(double(1.0));
    double dt_map = tbl_chp["dt_map"].value_or(double(10.0));

    REGULARIZATION* regularization = new REGULARIZATION(iter_max, g);

    //  select 
    //      Constant c, u>0
    //      Sine c(0,t)= sine at west boundary, uniform u>0
    //      Envelope  uniform u>0, sine within cosine envelope

    double dxinv = 1. / dx;
    int nx = int(Lx * dxinv) + 1 + 2; // nr nodes; including 2 virtual points

    int total_time_steps = int((tstop - tstart) / dt) + 1;  // [sec]
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

    log_file << "=== Used input variables ==============================" << std::endl;
    log_file << "[Domain]" << std::endl;
    log_file << "Lx = " << Lx << std::endl;

    log_file << std::endl << "[Time]" << std::endl;
    log_file << "tstart = " << tstart << std::endl;
    log_file << "tstop = " << tstop << std::endl;

    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "shape = " << shape_type << std::endl;

    log_file << std::endl << "[Boundary]" << std::endl;
    log_file << "bc_type = ";
    for (int i = 0; i < bc_type.size() - 1; ++i) { log_file << bc_type[i] << ", "; }
    log_file << bc_type[bc_type.size() - 1] << std::endl;
    log_file << "bc_vals = ";
    for (int i = 0; i < bc_vals.size() - 1; ++i) { log_file << bc_vals[i] << ", "; }
    log_file << bc_vals[bc_vals.size() - 1] << std::endl;
    log_file << "treg = " << treg << std::endl;

    log_file << std::endl << "[Physics]" << std::endl;
    log_file << "g = " << g << std::endl;
    log_file << "u_const = " << u_const << std::endl;
    log_file << "momentum_viscosity = " << momentum_viscosity << std::endl;
    log_file << "viscosity = " << visc_const << std::endl;

    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "dt  = " << dt << std::endl;
    log_file << "dx  = " << dx << std::endl;
    log_file << "c_psi  = " << c_psi << std::endl;
    log_file << "iter_max  = " << iter_max << std::endl;
    log_file << "theta  = " << theta << std::endl;
    log_file << "alpha  = " << alpha << std::endl;
    log_file << "eps_newton  = " << eps_newton << std::endl;
    log_file << "eps_bicgstab  = " << eps_bicgstab << std::endl;
    log_file << "regularization_init = " << regularization_init << std::endl;
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    log_file << "regularization_time = " << regularization_time << std::endl;

    log_file << std::endl << "[Output]" << std::endl;
    log_file << "dt_his = " << dt_his << std::endl;
    log_file << "dt_map = " << dt_map << std::endl;

    log_file << "-------------------------------------------------------" << std::endl;
    log_file << "Nodes  : " << nx << std::endl;
    log_file << "Volumes: " << (nx - 1) << std::endl;
    log_file << "CFL    : " << u_const * dt / dx << std::endl;
    STOP_TIMER(Writing log-file);  // but two write statements are not timed
    if (status != 0) {
        log_file.close();
        exit(1);
    }
    log_file << "=======================================================" << std::endl;

    double bc0;
    double bc1;

    std::vector<double> x(nx, 0.);  // x-coordinate
    std::vector<double> y(nx, 0.);  // y-coordinate
    std::vector<double> zb(nx, -10.0);                  // regularized bed level
    std::vector<double> zb_ini(nx, -10.0);              // initial bed level
    std::vector<double> phin(nx, 0.0);  // ln(conc) [-]
    std::vector<double> phip(nx, 0.);  // ln(conc) [-]
    std::vector<double> conc(nx, 0.);  // constituent [-]; exp(ln(conc))
    std::vector<double> q(nx, 0.);  // source
    std::vector<double> u(nx, 10.);  // advection velocity [m s-1], adv_diff_init_velocity
    std::vector<double> psi(nx, 0.);  //
    std::vector<double> eps(nx, 0.);  //
    std::vector<double> pe(nx, 0.);
    std::vector<double> delta_phi(nx, 0.);
    std::vector<double> visc_given(nx, visc_const);  // Initialize viscosity array with given value
    std::vector<double> visc_reg(nx, visc_const);  // Initialize given viscosity array with regularized value
    std::vector<double> visc(nx, visc_const);  // Viscosity array used for computation, adjusted for cell peclet number
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix
    std::vector<double> w_gp(3, 0.);  // weighting coefficients on grid point
    std::vector<double> w_cc(3, 0.);  // weighting coefficients on cell center
    std::vector<double> tmp(nx, 0.);  // dummy vector

    Eigen::VectorXd solution(nx);    // solution vector [c]^{n+1}
    Eigen::VectorXd rhs(nx);        // RHS vector [c]^n

    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    double alpha_bc = 2. * alpha - 1. / 2.;
    std::vector<double> w_nat(3, 0.0);
    w_nat[0] = 0.5 * (1.0 + alpha_bc);
    w_nat[1] = 0.5 * (1.0 - 2.0 * alpha_bc);
    w_nat[2] = 0.5 * alpha_bc;
    std::vector<double> w_ess(3, 0.0);
    w_ess[0] = w_nat[0];
    w_ess[1] = w_nat[1];
    w_ess[2] = w_nat[2];

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i-1) * dx;
    }

    // initial concentration
    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        (void)regularization->given_function(visc_reg, psi, visc_given, dx, c_psi);
        for (int i = 0; i < nx; ++i)
        {
            visc[i] = visc_reg[i] + std::abs(psi[i]);
        }
        STOP_TIMER(Regularization_init);
    }

    (void)adv_diff_init_concentration(mass, x, Lx, shape_conc, phin);
    (void)adv_diff_init_velocity(u, g, zb, x, shape_conc);  // g = 10., zb = -10 -> u = sqrt(g*zb) = 10.

    double time = tstart + dt * double(0);
    ////////////////////////////////////////////////////////////////////////////
    // Define map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    UGRID1D* map_file = new UGRID1D();

    std::vector<std::string> map_names;
    map_names.push_back("phin_1d");
    map_names.push_back("u_1d");
    map_names.push_back("visc_1d");
    map_names.push_back("psi_1d");
    map_names.push_back("pe_1d");
    map_names.push_back("conc_1d");
    map_file = create_map_file(nc_mapfilename, model_title, x, map_names);

    // Put data on map file
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    for (size_t i = 0; i < phin.size(); ++i)
    {
        conc[i] = std::exp(phin[i]);
    }
    map_file->put_time_variable(map_names[0], nst_map, phin);
    map_file->put_time_variable(map_names[1], nst_map, u);
    map_file->put_time_variable(map_names[2], nst_map, visc);
    map_file->put_time_variable(map_names[3], nst_map, psi);
    map_file->put_time_variable(map_names[4], nst_map, pe);
    map_file->put_time_variable(map_names[5], nst_map, conc);

    ////////////////////////////////////////////////////////////////////////////
    // Define time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations (corrected for virtual points)
    int p_a = 1;
    int p_b = nx / 4 + 1;
    int p_c = nx / 2;
    int p_d = 3 * nx / 4 - 1;
    int p_e = nx - 2;
    std::vector<double> x_obs = { x[p_a], x[p_b], x[p_c], x[p_d], x[p_e] };
    std::vector<double> y_obs = { y[p_a], y[p_b], y[p_c], y[p_d], y[p_e] };

    int nsig = 0;
    for (int i = 0; i < x_obs.size(); ++i)
    {
        nsig = std::max(nsig, (int)std::log10(x_obs[i]));
    }
    nsig += 1;
    std::vector<std::string> obs_stations;
    obs_stations.push_back(setup_obs_name(x[p_a], y[p_a], nsig, ": West boundary"));
    obs_stations.push_back(setup_obs_name(x[p_b], y[p_b], nsig, ": Halfway to west boundary"));
    obs_stations.push_back(setup_obs_name(x[p_c], y[p_c], nsig, ": Centre"));
    obs_stations.push_back(setup_obs_name(x[p_d], y[p_d], nsig, ": Halfway to east boundary"));
    obs_stations.push_back(setup_obs_name(x[p_e], y[p_e], nsig, ": East boundary"));

    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_phi_name("lnc");
    his_file->add_variable(his_phi_name, "", "Phi, ln(c)", "-");
    std::string his_expphi_name("exp(lnc)");
    his_file->add_variable(his_expphi_name, "", "Concentration", "-");

    std::string his_newton_iter_name("newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "", "Newton iterations", "-");
    std::string his_lin_solv_iter_name("lin_solver_iterations");
    his_file->add_variable_without_location(his_lin_solv_iter_name, "", "Lin. solver iterations", "-");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, time);

    std::vector<double> his_values = { phin[p_a], phin[p_b], phin[p_c], phin[p_d], phin[p_e] };
    his_file->put_variable(his_phi_name, nst_his, his_values);
    for (size_t i = 0; i < his_values.size(); ++i)
    {
        his_values[i] = std::exp(his_values[i]);
    }
    his_file->put_variable(his_expphi_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);
    his_values = { 0 };
    his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);

    ////////////////////////////////////////////////////////////////////////////

    (void)adv_diff_source(q, select_src);

    //(void)adv_diff_boundary_condition(phin[0], phin[nx - 1], time, treg, select_bc);
    for (int i = 0; i < nx; i++)
    {
        phip[i] = phin[i];  // place the old solution in the new one, needed for the iteration procedure
    }

    Eigen::SparseMatrix<double> A(nx, nx);
    for (int i = 0; i < nx; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
        solution[i] = 0.0;
    }

    // start time loop
    std::cout << "Start time-loop" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "tstart= " << tstart + time << ";   tstop= " << tstart + tstop << ";   dt= " << dt << ";   dx= " << dx << std::endl;

    double dc_max = 0.0;
    int dc_maxi = 0;
    STOP_TIMER(Simulation initialization);
    START_TIMER(Time loop);
    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        int used_newton_iter = 0;
        int used_lin_solv_iter = 0;
        time = dt * double(nst);
#if defined(DEBUG)
        log_file << "=== Start time step ===================================" << std::endl;
        log_file << std::setprecision(4) << time << std::endl;
#endif
        adv_diff_boundary_condition(bc0, bc1, bc_vals[BC_WEST], bc_vals[BC_EAST], time, treg, bnd_type);
        if (regularization_time)
        {
            START_TIMER(Regularization_time_loop);
            if (momentum_viscosity)
            {
                //for (size_t i = 0; i < phin.size(); ++i)
                //{
                //    conc[i] = std::exp(phip[i]);
                //}
                (void)regularization->artificial_viscosity(psi, u, phip, c_psi, dx);
                for (int i = 0; i < nx; ++i)
                {
                    psi[i] = std::max(1e-25, psi[i]);
                    visc[i] = visc_reg[i] + std::abs(psi[i]);
                }
            }
            STOP_TIMER(Regularization_time_loop);
        }
        START_TIMER(Newton iteration);
        for (int iter = 0; iter < iter_max; ++iter)
        {
            used_newton_iter += 1;
            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                if (momentum_viscosity)
                {
                    //for (size_t i = 0; i < phin.size(); ++i)
                    //{
                    //    conc[i] = std::exp(phip[i]);
                    //}
                    (void)regularization->artificial_viscosity(psi, u, phip, c_psi, dx);
                    for (int i = 0; i < nx; ++i)
                    {
                        psi[i] = std::max(1e-25, psi[i]);
                        visc[i] = visc_reg[i] + std::abs(psi[i]);
                    }
                }
                STOP_TIMER(Regularization_iter_loop);
            }
            if (nst == 1 && iter == 0)
            {
                START_TIMER(Matrix initialization);
            }
            START_TIMER(Matrix set up);
            //
            // interior nodes
            //
            for (int i = 1; i < nx - 1; ++i)
            {
                double phin_im1 = phin[i - 1];
                double phin_i = phin[i];
                double phin_ip1 = phin[i + 1];
                double phin_im12 = 0.5 * (phin[i - 1] + phin[i]);
                double phin_ip12 = 0.5 * (phin[i] + phin[i + 1]);

                double phip_im1 = phip[i - 1];
                double phip_i = phip[i];
                double phip_ip1 = phip[i + 1];
                double phip_im12 = 0.5 * (phip[i - 1] + phip[i]);
                double phip_ip12 = 0.5 * (phip[i] + phip[i + 1]);

                double phitheta_im12 = theta * phip_im12 + (1.0 - theta) * phin_im12;
                double phitheta_ip12 = theta * phip_ip12 + (1.0 - theta) * phin_ip12;

                // d(c)/dt
                double phin_im14 = scv(phin_i, phin_im1);
                double phin_ip14 = scv(phin_i, phin_ip1);
                double phip_im14 = scv(phip_i, phip_im1);
                double phip_ip14 = scv(phip_i, phip_ip1);
                A.coeffRef(i, i - 1) = 0.5 * dx * dtinv * std::exp(phip_im14) * 0.25;
                A.coeffRef(i, i) = 0.5 * dx * dtinv * 0.75 * (std::exp(phip_im14) + std::exp(phip_ip14));
                A.coeffRef(i, i + 1) = 0.5 * dx * dtinv * std::exp(phip_ip14) * 0.25;
                rhs[i] = -(
                    0.5 * dx * dtinv * (std::exp(phip_im14) - std::exp(phin_im14))
                    + 0.5 * dx * dtinv * (std::exp(phip_ip14) - std::exp(phin_ip14))
                    );
                //
                // d(uc)/dx
                double u_im12 = 0.5 * (u[i - 1] + u[i]);
                double u_ip12 = 0.5 * (u[i] + u[i + 1]);
                A.coeffRef(i, i - 1) += -0.5 * u_im12 * theta * std::exp(phip_im12);
                A.coeffRef(i, i) +=
                    -0.5 * u_im12 * theta * std::exp(phip_im12)
                    + 0.5 * u_ip12 * theta * std::exp(phip_ip12);
                A.coeffRef(i, i + 1) += 0.5 * u_ip12 * theta * std::exp(phip_ip12);
                rhs[i] += -(
                    u_ip12 * std::exp(phip_ip12) - u_im12 * std::exp(phip_im12)
                    );
                //
                // d(visc(c/dx))/dx
                if (momentum_viscosity)
                {
                    double visc_im12 = 0.5 * (visc[i - 1] + visc[i]);
                    double visc_ip12 = 0.5 * (visc[i] + visc[i + 1]);
                    A.coeffRef(i, i - 1) += - theta * visc_im12 * dxinv * std::exp(phip_im1);
                    A.coeffRef(i, i    ) += + theta * visc_im12 * dxinv * std::exp(phip_im1) +  theta * visc_ip12 * dxinv * std::exp(phip_ip1);
                    A.coeffRef(i, i + 1) += - theta * visc_ip12 * dxinv * std::exp(phip_ip1);
                    rhs[i] += -(
                        - visc_ip12 * dxinv * (std::exp(phip_ip1) - std::exp(phip_i))
                        + visc_im12 * dxinv * (std::exp(phip_i) - std::exp(phip_im1))
                        );
                }
            }
            {
                //
                // wwest boundary (u>0; essential boundary)
                //
                int i = 0;

                double phip_i = phip[i];             // = c^{n+1,p}_{i}
                double phip_ip1 = phip[i + 1];       // = c^{n+1,p}_{i+1}
                double phip_ip2 = phip[i + 2];       // = c^{n+1,p}_{i+2}

                w_ess[0] = 1. / 12.;
                w_ess[1] = 10. / 12.;
                w_ess[2] = 1. / 12.;
                w_ess[0] = w_nat[0];
                w_ess[1] = w_nat[1];
                w_ess[2] = w_nat[2];
                A.coeffRef(i, i    ) = w_ess[0];
                A.coeffRef(i, i + 1) = w_ess[1];
                A.coeffRef(i, i + 2) = w_ess[2];
                rhs[i] = -25. - (w_ess[0] * phip_i + w_ess[1] * phip_ip1 + w_ess[2] * phip_ip2);  // if u>0 this is upwind
            }
            {
                //
                // eeast boundary (u>0; natural boundary)
                //
                int i = nx - 1;

                double phin_i = phin[i];             // = h^{n}_{i}
                double phin_im1 = phin[i - 1];       // = h^{n}_{i-1}
                double phin_im2 = phin[i - 2];       // = h^{n}_{i-2}
                double phip_i = phip[i];             // = h^{n+1,p}_{i}
                double phip_im1 = phip[i - 1];       // = h^{n+1,p}_{i-1}
                double phip_im2 = phip[i - 2];       // = h^{n+1,p}_{i-2}
                double phitheta_i = theta * phip_i + (1.0 - theta) * phin_i;
                double phitheta_im1 = theta * phip_im1 + (1.0 - theta) * phin_im1;

                // Outflow boundary (natural boundary)
                double dcdt = dtinv * (
                    w_nat[0] * (phip_i - phin_i) +
                    w_nat[1] * (phip_im1 - phin_im1) +
                    w_nat[2] * (phip_im2 - phin_im2)
                    );

                double u_im12 = 0.5 * (u[i] + u[i - 1]);
                double udcdx =  u_im12 * (phitheta_i - phitheta_im1) * dxinv;

                A.coeffRef(i, i    ) = dtinv * w_nat[0] + theta * dxinv * u_im12;
                A.coeffRef(i, i - 1) = dtinv * w_nat[1] - theta * dxinv * u_im12;
                A.coeffRef(i, i - 2) = dtinv * w_nat[2];
                rhs[i] = - (dcdt + udcdx);
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
            solver.setTolerance(eps_bicgstab);
            //solution = solver.solve(rhs);
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
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(2) << std::scientific << time
                         << "    BiCGstab iterations: " << solver.iterations()
                         << "    estimated error:" << solver.error() 
                         << std::endl;
            }

            // The new solution is the previous iterant plus the delta
            dc_max = 0.0;
            dc_maxi = 0;
            for (int i = 0; i < nx; ++i)
            {
                phip[i] += solution[i];
                delta_phi[i] = solution[i]; 
                if (dc_max < std::abs(delta_phi[i]))
                {
                    dc_max = std::abs(delta_phi[i]);
                    dc_maxi = i;
                }
            }

            if (logging == "matrix")
            {
                log_file << "=== Matrix ============================================" << std::endl;
                log_file << std::setprecision(4) << std::scientific << Eigen::MatrixXd(A) << std::endl;
                log_file << "=== RHS ===============================================" << std::endl;
                log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
            }
            if (logging == "matrix")
            {
                log_file << "=== cn, delta_c, cp ===================================" << std::endl;
                for (int i = 0; i < nx; ++i)
                {
                    log_file << std::setprecision(15) << std::scientific << phin[i] << "  -----  " << delta_phi[i] << "  -----  " << phip[i]  << std::endl;
                }
                //log_file << "=== Eigen values ======================================" << std::endl;
                //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).eigenvalues() << std::endl;
                log_file << "====== End iteration ==================================" << std::endl;
                log_file << "time [sec]: " << dt * double(nst)
                         << ";    iterations     : " << solver.iterations()
                         << ";    estimated error: " << solver.error() 
                         << std::endl;
            }

            STOP_TIMER(Set solution);
            used_lin_solv_iter = std::max(used_lin_solv_iter, (int)solver.iterations());
            if (dc_max < eps_newton)
            {
                break;
            }
        }
        STOP_TIMER(Newton iteration);
        if (stationary)
        {
            std::cout << "stationary solution " << std::endl;
            log_file << "stationary solution " << std::endl;
            log_file << std::setprecision(8) << std::scientific
                << "    Iter: " << used_newton_iter
                << "    Delta c^{n + 1,p + 1}: " << dc_max << " at: " << dc_maxi
                << std::endl;
        }
        else
        {
            if (std::fmod(nst, 100) == 0)
            {
                std::cout << std::fixed << std::setprecision(2) << tstart + time << ";   " << tstart + tstop << std::endl;
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(2) << std::scientific << time
                    << std::setprecision(8) << std::scientific
                    << "    Newton iterations  : " << used_newton_iter
                    << "    Delta c^{n + 1,p + 1}: " << dc_max << " at: " << dc_maxi
                    << std::endl;
            }
        }
        if (used_newton_iter == iter_max)
        {
            if (dc_max > eps_newton)
            {
                log_file << "    ----    maximum number of iterations reached, probably not converged" << std::endl;
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            phin[i] = phip[i];
        }
        if (logging == "matrix")
        {
            log_file << "=== cp, cn, cp-cn =====================================" << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                log_file << std::setprecision(15) << std::scientific << phip[i] << ",   " << phin[i] << ",   " << phip[i] - phin[i] << std::endl;
            }
        }

        // Map-files
        if (std::fmod(nst, wrimap) == 0)
        {
            // Put data on time map file
            START_TIMER(Writing map-file);
            nst_map++;
            if (stationary)
            {
                map_file->put_time(nst_map, double(nst));
            }
            else
            {
                map_file->put_time(nst_map, double(nst)* dt);
            }
            for (size_t i = 0; i < phin.size(); ++i)
            {
                conc[i] = std::exp(phin[i]);
            }
            map_file->put_time_variable(map_names[0], nst_map, phin);
            map_file->put_time_variable(map_names[1], nst_map, u);
            map_file->put_time_variable(map_names[2], nst_map, eps);
            map_file->put_time_variable(map_names[3], nst_map, psi);
            map_file->put_time_variable(map_names[4], nst_map, pe);
            map_file->put_time_variable(map_names[5], nst_map, conc);
            STOP_TIMER(Writing map-file);
        }

        // His-files
        if (std::fmod(nst, wrihis) == 0)
        {
            START_TIMER(Writing his-file);
            nst_his++;
            if (stationary)
            {
                his_file->put_time(nst_his, double(nst));
            }
            else
            {
                his_file->put_time(nst_his, double(nst) * dt);
            }

            std::vector<double> his_values = { phin[p_a], phin[p_b], phin[p_c], phin[p_d],  phin[p_e] };
            his_file->put_variable(his_phi_name, nst_his, his_values);
            for (size_t i = 0; i < his_values.size(); ++i)
            {
                his_values[i] = std::exp(his_values[i]);
            }
            his_file->put_variable(his_expphi_name, nst_his, his_values);

            his_values = { double(used_newton_iter) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            his_values = { double(used_lin_solv_iter) };
            his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);
            STOP_TIMER(Writing his-file);
        }
    } // End of the time loop
    STOP_TIMER(Time loop);

    log_file.close();
    (void)map_file->close();
    (void)his_file->close();
    STOP_TIMER(Writing log-file);

    STOP_TIMER(main);
    PRINT_TIMER(timing_filename.data());

    std::chrono::duration<int, std::milli> timespan(1000);
    std::this_thread::sleep_for(timespan);
    return 0;
}
//------------------------------------------------------------------------------
std::vector<double> bc_val(double xi)
{ 
    // from document poisseuille flow and email 14 dec 2023
    // if boundary at c_1 then xi=1

    std::vector<double> ret_val(3, 0.0);
    ret_val[0] = 13. / 12. - 1.5 * xi + 0.5 * xi * xi;
    ret_val[1] = -1. / 6. + 2. * xi - xi * xi;
    ret_val[2] = 1. / 12. - 0.5 * xi + 0.5 * xi * xi;

    return ret_val;
}
std::vector<double> quadratic_interpolation_1d(std::vector<double> x, std::vector<double> f, int bc_east_west, double x_loc)
{
    std::vector<double> x_gp(3, 0.0);
    std::vector<double> f_gp(3, 0.0);
    if (bc_east_west == BC_WEST)
    {
        x_gp[0] = x[0];
        x_gp[1] = x[1];
        x_gp[2] = x[2];
        f_gp[0] = f[0];
        f_gp[1] = f[1];
        f_gp[2] = f[2];
        x_loc = (1.0 - x_loc) * x[0] + x_loc * x[1];
    }
    if (bc_east_west == BC_EAST)
    {
        int nx = (int)x.size();
        x_gp[0] = x[nx - 1];
        x_gp[1] = x[nx - 2];
        x_gp[2] = x[nx - 3];
        f_gp[0] = f[0];
        f_gp[1] = f[1];
        f_gp[2] = f[2];
        x_loc = (1.0 - x_loc) * x[nx - 1] + x_loc * x[nx - 2];
    }
    //double b_1 = (f_gp[1] - f_gp[0]) / (x_gp[1] - x_gp[0]);
    //double b_2 = ((f_gp[2] - f_gp[1]) / (x_gp[2] - x_gp[1]) - (f_gp[1] - f_gp[0]) / (x_gp[1] - x_gp[0])) / (x_gp[2] - x_gp[0]);
    //double f_xloc = f_gp[0] + b_1 * (x_loc - x_gp[0]) + b_2 * (x_loc - x_gp[0]) * (x_loc - x_gp[1]);
    //double dfdx_xloc = b_1 + b_2 * (x_loc - x_gp[0]) + b_2 * (x_loc - x_gp[1]);
    //double d2fdx2_xloc = 2 * b_2;

    // Calculate interpolation weights at x_loc such that
    // f(x_loc) = w[0] * f_gp[0] + w[1] * f_gp[1] + w[2] * f_gp[2] in non - equidistant grids
    std::vector<double> w(3, 0.0);
    w[0] = 1 - (x_loc - x_gp[0]) / (x_gp[1] - x_gp[0]) + ((x_loc - x_gp[0]) * (x_loc - x_gp[1])) / ((x_gp[1] - x_gp[0]) * (x_gp[2] - x_gp[0]));
    w[1] = (x_loc - x_gp[0]) / (x_gp[1] - x_gp[0]) - ((x_loc - x_gp[0]) * (x_loc - x_gp[1])) / ((x_gp[1] - x_gp[0]) * (x_gp[2] - x_gp[0])) - ((x_loc - x_gp[0]) * (x_loc - x_gp[1])) / ((x_gp[2] - x_gp[1]) * (x_gp[2] - x_gp[0]));
    w[2] = ((x_loc - x_gp[0]) * (x_loc - x_gp[1])) / ((x_gp[2] - x_gp[1]) * (x_gp[2] - x_gp[0]));

    return w;
}
std::vector<double> hermite_interpolation_1d(std::vector<double> x, std::vector<double> f, int bc_east_west, double x_loc)
{
    std::vector<double> x_gp(3, 0.0);
    std::vector<double> f_gp(3, 0.0);
    if (bc_east_west == BC_WEST)
    {
        x_gp[0] = x[0];
        x_gp[1] = x[1];
        x_gp[2] = x[2];
        f_gp[0] = f[0];
        f_gp[1] = f[1];
        f_gp[2] = f[2];
        x_loc = (1.0 - x_loc) * x[0] + x_loc * x[1];
    }
    if (bc_east_west == BC_EAST)
    {
        int nx = (int)x.size();
        x_gp[0] = x[nx - 1];
        x_gp[1] = x[nx - 2];
        x_gp[2] = x[nx - 3];
        f_gp[0] = f[nx - 1];
        f_gp[1] = f[nx - 2];
        f_gp[2] = f[nx - 3];
        x_loc = (1.0 - x_loc) * x[nx - 1] + x_loc * x[nx - 2];
    }

    std::vector<double> x_cc(2, 0.0);
    std::vector<double> f_cc(2, 0.0);
    std::vector<double> dfdx_cc(2, 0.0);

    // Using cell - centered values and derivatives
    x_cc[0] = (x_gp[0] + x_gp[1]) * 0.5;
    x_cc[1] = (x_gp[1] + x_gp[2]) * 0.5;
    f_cc[0] = (f_gp[0] + f_gp[1]) * 0.5;
    f_cc[1] = (f_gp[1] + f_gp[2]) * 0.5;
    dfdx_cc[0] = (f_gp[1] - f_gp[0]) / (x_gp[1] - x_gp[0]);
    dfdx_cc[1] = (f_gp[2] - f_gp[1]) / (x_gp[2] - x_gp[1]);

    // Define coordinate t on which basis functions are constructed
    double t = (x_loc - x_cc[0]) / (x_cc[1] - x_cc[0]);

    // Define basis functions
    double h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
    double h10 = t * (1.0 - t) * (1.0 - t);
    double h01 = (3.0 - 2.0 * t) * t * t;
    double h11 = (t - 1.0) * t * t;

    // Define interpolation polynomial
    //double p = h00 * f_cc[0] + h10 * (x_cc[1] - x_cc[0]) * dfdx_cc[0] + h01 * f_cc[1] + h11 * (x_cc[1] - x_cc[0]) * dfdx_cc[1];

    // Calculate weights such that p = w[0] * f[0] + w[1] * f[1] + w[0] * f[0] in non - equidistant grids
    // f_cc = 0.5 * (f_gp[0] + f_gp[1]), f_cc[1] = 0.5 * (f_gp[1] + f_gp[0])
    std::vector<double> w(3, 0.0);
    w[0] = h00 * 0.5 - h10 * (x_cc[1] - x_cc[0]) / (x_gp[1] - x_gp[0]);
    w[1] = h00 * 0.5 + h10 * (x_cc[1] - x_cc[0]) / (x_gp[1] - x_gp[0]) + h01 * 0.5 - h11 * (x_cc[1] - x_cc[0]) / (x_gp[0] - x_gp[1]);
    w[2] = h01 * 0.5 + h11 * (x_cc[1] - x_cc[0]) / (x_gp[0] - x_gp[1]);

    return w;
}
double scv(double h1, double h2)
{
    // value at subcontrol volume
    return 0.25 * (3. * h1 + h2);
}
std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name)
{
    std::stringstream ss_x;
    std::stringstream ss_y;
    ss_x << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << x_obs;
    ss_y << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << y_obs;
    return( "x=" + ss_x.str() + obs_name);  // ss_y not used
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
        status = 0;
    }
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
        status = 0;
    }
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
        status = 0;
    }
    return status;
}

