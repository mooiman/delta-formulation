//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D Burgers equation, fully implicit with delta-formulation and Modified Newton iteration 
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
#include <Eigen/Eigenvalues> 
#include <include/KDtree.hpp>
#include <toml.h>

#include "boundary_condition.h"
#include "adv_diff_init_velocity.h"
#include "adv_diff_linear_operator.h"
#include "cfts.h"
#include "compile_date_and_time.h"
#include "data_input_struct.h"
#include "definition_map_file.h"
#include "main_version.h"
#include "observation_stations.h"
#include "perf_timer.h"
#include "print_matrix.h"
#include "read_input_toml_file.h"
#include "regularization.h"
#include "ugrid1d.h"

void GetArguments(long argc, char** argv, std::filesystem::path & file_name);
int write_used_input(struct _data_input data, std::ofstream & log_file);
int set_his_values(std::vector<_ObservationPoint>& obs_points, std::vector<double> & array, std::vector<double>& his_values);
std::vector<double> bc_val(double x_loc);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);
//
// Solve the linear Burgers equation
//
// dc/dt + d(uc)/dx - nu d2(u)dx2 = 0
//

int main(int argc, char* argv[])
{
    bool stationary = false;
    std::filesystem::path toml_file_name("---not-defined---");

    int status = -1;

    int BC_WEST = 0;
    int BC_EAST = 1;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path start_dir;
    std::filesystem::path output_dir;
    std::filesystem::path input_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();
    start_dir = std::filesystem::current_path();

    if (argc == 3)
    {
        (void)GetArguments(argc, argv, toml_file_name);
        if (!std::filesystem::is_regular_file(toml_file_name))
        {
            std::cout << "----------------------------" << std::endl;
            std::cout << "Input file \'" << toml_file_name << "\' can not be opened." << std::endl;
            std::cout << "Press Enter to finish";
            //std::cin.ignore();
            std::chrono::duration<int, std::milli> timespan(3000);
            std::this_thread::sleep_for(timespan);
            exit(1);
        }
        input_dir = std::filesystem::absolute(toml_file_name);
        input_dir.remove_filename();
        output_dir = input_dir / "output" / "";
        std::filesystem::create_directory(output_dir);
    }
    else
    {
        std::cout << "======================================================" << std::endl;
        std::cout << "Executable compiled : " << compileDateTime() << std::endl;
        std::cout << "Git commit time/hash: " << getbuildstring_main() << std::endl;
        std::cout << std::endl;
        std::cout << "usage: wave_1d.exe --toml <input_file>" << std::endl;
        std::cout << "======================================================" << std::endl;
        //std::cin.ignore();
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        exit(1);
    }
    const std::chrono::zoned_time now{ std::chrono::current_zone(), std::chrono::system_clock::now() };
    auto start_date_time = std::format("{:%F %H:%M:%OS %Oz}", now);
    std::cout << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Start time          : " << start_date_time << std::endl;
    std::cout << "Executable compiled : " << compileDateTime() << std::endl;
    std::cout << "Git commit time/hash: " << getbuildstring_main() << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Start directory     : " << start_dir << std::endl;
    std::cout << "Input directory     : " << input_dir << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;

    START_TIMERN(Main);
    START_TIMER(Writing log-file);
    START_TIMER(Initialization);

    std::string out_file;
    std::stringstream ss;

    struct _data_input input_data;
    try
    {
        input_data = read_toml_file(input_dir, toml_file_name);
    }
    catch (const toml::parse_error& err) 
    {
        std::cout << "\n+++++++++++++++++++++\nTOML parse error\n"
                  << err 
                  << "\n+++++++++++++++++++++\n\n";
        std::chrono::duration<int, std::milli> timespan(5000);
        std::this_thread::sleep_for(timespan);
        exit(1);
    }

    // convert given time to seconds
    if (input_data.time.tunit != "s")
    {
        if (input_data.time.tunit == "m")
        {
            input_data.time.tstart = 60.0 * input_data.time.tstart;
            input_data.time.tstop  = 60.0 * input_data.time.tstop;
        }
        if (input_data.time.tunit == "h")
        {
            input_data.time.tstart = 3600.0 * input_data.time.tstart;
            input_data.time.tstop  = 3600.0 * input_data.time.tstop;
        }
        if (input_data.time.tunit == "d")
        {
            input_data.time.tstart = 24.0 * 3600.0 * input_data.time.tstart;
            input_data.time.tstop  = 24.0 * 3600.0 * input_data.time.tstop;
        }
    }
    if (input_data.numerics.dt == 0.0) { stationary = true;  }
    if (stationary) { 
        input_data.boundary.treg = 0.0; 
        input_data.numerics.theta = 1.0;
    }

    ss << "_dx" << input_data.numerics.dx << "_dt" << input_data.numerics.dt;

    out_file = output_dir.string() + "burgers_eq" + ss.str();

    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");

    std::ofstream log_file;
    log_file.open(log_filename);

    std::cout << "=== Input file =======================================" << std::endl;
    std::cout << std::filesystem::absolute(toml_file_name) << std::endl;
    std::cout << "======================================================" << std::endl;

    log_file << "======================================================" << std::endl;
    log_file << "Start time          : " << start_date_time << std::endl;
    log_file << "Executable compiled : " << compileDateTime() << std::endl;
    log_file << "Git commit time/hash: " << getbuildstring_main() << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;

    //  select 
    //      Constant c, u>0
    //      Sine c(0,t)= sine at west boundary, uniform u>0
    //      Envelope  uniform u>0, sine within cosine envelope

    double dxinv = 1./input_data.numerics.dx;                                 // invers grid size [m]
    size_t nx = int(input_data.domain.Lx * dxinv) + 1 + 2 ;                    // number of nodes; including 2 virtual points
    size_t ny = 1; 
    size_t nxny = nx * ny; 

    int total_time_steps = int((input_data.time.tstop - input_data.time.tstart) / input_data.numerics.dt) + 1;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    double dtinv_pseu;                                    // Inverse of pseudo time step, if dt==0 then stationary solution [1/s]
    int wrihis;                                           // write interval to his-file
    int wrimap;                                           // write interval to map-file

    if (stationary)
    {
        input_data.numerics.dt = 0.0;                     // Time step size [s]
        dtinv = 0.0;                                      // stationary solution
        dtinv_pseu = 0.0;                                 // Inverse of pseudo time step, if dt==0 then stationary solution [1/s]
        input_data.boundary.eps_bc_corr = 1.0;
        input_data.numerics.iter_max = 2 * input_data.numerics.iter_max;
        input_data.numerics.theta = 1.0;                                      // Stationary solution
        input_data.time.tstop = 1.;
        total_time_steps = 2;                             // initial step (step 1), stationary result (step 2)
        input_data.boundary.treg = 0.0;                   // Thatcher-Harleman return time [s], when zero supply boundary value immediately
        wrihis = 1;                                       // write interval to his-file
        wrimap = 1;                                       // write interval to map-file
    }
    else
    {
        dtinv = 1. / input_data.numerics.dt;               // Inverse of dt [1/s]
        dtinv_pseu = 0.0;                                  // Inverse of pseudo time step, if dt==0 then stationary solution [1/s]
        wrihis = std::max(int(input_data.numerics.dt * dtinv), int(input_data.output.dt_his * dtinv));      // write interval to his-file (every delta t)
        if (input_data.output.dt_map == 0.0)
        {
            wrimap = total_time_steps - 1;  // write only first and last time step
        }
        else
        {
            wrimap = std::max(int(input_data.numerics.dt * dtinv), int(input_data.output.dt_map * dtinv));     // write interval to map-file (every 1 sec , or every delta t)
        }
    }
    std::string solver_name("bicgstab");

    std::string model_title("not defined");
    if (input_data.physics.do_viscosity) { 
        model_title = "Burgers (viscid)";
    }
    else 
    { 
        model_title = "Burgers (inviscid), 1D";
    }

    model_title += ", " + solver_name;

    REGULARIZATION* regularization = new REGULARIZATION(input_data.numerics.iter_max, input_data.physics.g, input_data.log.logging);

    log_file << "=== Used input variables ==============================" << std::endl;
    status = write_used_input(input_data, log_file);
    log_file << "=======================================================" << std::endl;

    double bc0;
    double bc1;

    // Copy input data to loccal data
    std::string logging = input_data.log.logging;

    double Lx       = input_data.domain.Lx;
    double x_origin = input_data.domain.x_origin;

    double eps_bc_corr = input_data.boundary.eps_bc_corr;
    double treg = input_data.boundary.treg;
    std::vector<std::string> bc_type = input_data.boundary.bc_type;
    std::vector<std::string> bc_vars = input_data.boundary.bc_vars;
    std::vector<std::string> bc_signals = input_data.boundary.bc_signals;
    std::vector<double> bc_vals = input_data.boundary.bc_vals;

    double u_initial = input_data.initial.u_initial;
    std::vector<std::string> ini_vars = input_data.initial.ini_vars;

    double dt = input_data.numerics.dt;
    double dx = input_data.numerics.dx;
    double c_psi = input_data.numerics.c_psi;
    double eps_bicgstab = input_data.numerics.eps_bicgstab;
    double eps_newton = input_data.numerics.eps_newton;
    double theta = input_data.numerics.theta;
    int iter_max = input_data.numerics.iter_max;
    std::string linear_solver = input_data.numerics.linear_solver;
    bool regularization_init = input_data.numerics.regularization_init;
    bool regularization_iter = input_data.numerics.regularization_iter;
    bool regularization_time = input_data.numerics.regularization_time;

    double visc_const = input_data.physics.visc_const;
    bool do_convection = input_data.physics.do_convection;
    bool do_viscosity = input_data.physics.do_viscosity;
    bool do_source =input_data.physics.do_source;
    std::string src_type = input_data.physics.src_type;

    double tstart = input_data.time.tstart;
    double tstop = input_data.time.tstop;

    std::vector<double> x(nx, 0.);  // x-coordinate
    std::vector<double> y(nx, 0.);  // y-coordinate
    std::vector<double> un(nx, u_initial);  // velocity [m s-1]
    std::vector<double> up(nx, 0.0);  // velocity [m s-1]
    std::vector<double> utheta(nx, 0.0);  // velocity at time level theta [(1-theta) * un + theta * up] [m s-1]
    std::vector<double> psi(nx, 0.);  //
    std::vector<double> delta_u(nx, 0.);
    std::vector<double> pe(nx, 0.);
    std::vector<double> cfl(nx, 0.);
    std::vector<double> visc_given(nx, visc_const);  // Initial viscosity array with given value
    std::vector<double> visc_reg(nx, visc_const);  // Initial given viscosity array with regularized value
    std::vector<double> visc(nx, visc_const);  // Viscosity array used for computation, adjusted with artificial viscosity
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix
    std::vector<double> w_bc(3, 0.);  // weighting coefficients on boundary
    std::vector<double> w_gp(3, 0.);  // weighting coefficients on grid point
    std::vector<double> w_cc(3, 0.);  // weighting coefficients on cell center
    std::vector<double> tmp(nx, 0.);  // dummy vector

    Eigen::VectorXd solution(nx);    // solution vector [c]^{n+1}
    Eigen::VectorXd rhs(nx);        // RHS vector [c]^n

    double alpha = 1./8.;
    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    double alpha_bc = 2. * alpha - 1. / 2.;
    std::vector<double> w_nat(3, 0.0);
    w_nat[0] = 0.5 * (1.0 + alpha_bc);
    w_nat[1] = 0.5 * (1.0 - 2.0 * alpha_bc);
    w_nat[2] = 0.5 * alpha_bc;
    std::vector<double> w_ess(3, 0.0);
    w_ess[0] = 1./12.;
    w_ess[1] = 10./12.;
    w_ess[2] = 1./12.;
    w_ess[0] = 11./24.;
    w_ess[1] = 14./24.;
    w_ess[2] = -1./24.;

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i - 1) * dx + x_origin;
    }
    adv_diff_init_velocity(un, u_initial, x, ini_vars[0]);
    //  Create kdtree, needed to locate the observation points
    std::vector<std::vector<double>> xy_points;
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::vector<double> point = {x[i], y[i]};
        xy_points.push_back(point);
    }
    KDTree xy_tree(xy_points);

    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        regularization->given_function(visc_reg, psi, visc_given, dx, c_psi, log_file);
        for (size_t i = 0; i < nx; ++i)
        {
            visc[i] = visc_reg[i] + std::abs(psi[i]);
        }
        STOP_TIMER(Regularization_init);
    }
    for (int i = 0; i < nx; ++i)
    {
        cfl[i] = up[i] * dt / dx;
    }
    if (do_viscosity)
    {
        for (int i = 0; i < nx; ++i)
        {
            pe[i] = up[i] * dx / visc[i];
        }
    }

    double time = tstart + dt * double(0);
    ////////////////////////////////////////////////////////////////////////////
    // Create map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    UGRID1D* map_file = new UGRID1D();

    std::vector<std::string> map_names;
    map_names.push_back("u_1d");
    map_names.push_back("visc_1d");
    map_names.push_back("psi_1d");
    map_names.push_back("pe_1d");
    map_names.push_back("cfl_1d");
    map_file = create_map_file(nc_mapfilename, model_title, x, map_names);
    
    // Put data on map file
    START_TIMER(Writing map-file);
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    map_file->put_time_variable(map_names[0], nst_map, un);
    map_file->put_time_variable(map_names[1], nst_map, visc);
    map_file->put_time_variable(map_names[2], nst_map, psi);
    map_file->put_time_variable(map_names[3], nst_map, pe);
    map_file->put_time_variable(map_names[4], nst_map, cfl);
    STOP_TIMER(Writing map-file);

    // End define map file
    ////////////////////////////////////////////////////////////////////////////
    // Create time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    std::vector<std::string> obs_station_names;
    status = def_observation_stations(obs_station_names, input_data.obs_points, xy_tree, x, y, nx, ny);

    std::vector<double> x_obs;
    std::vector<double> y_obs;
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        _ObservationPoint obs = input_data.obs_points[i];
        x_obs.push_back(obs.x);
        y_obs.push_back(obs.y);
    }
    his_file->add_stations(obs_station_names, x_obs, y_obs);

    his_file->add_time_series();

    std::string his_u_name("u_velocity");
    std::string his_visc_name("his_visc_name");
    std::string his_psi_name("his_psi_name");
    std::string his_peclet_name("his_peclet_name");
    std::string his_cfl_name("his_cfl_name");

    his_file->add_variable(his_u_name, "", "Velocity", "m s-1");
    his_file->add_variable(his_cfl_name, "", "CFL (u.dt/dx)", "-");
    if (do_viscosity) 
    {
        his_file->add_variable(his_visc_name, "", "Viscosity (used)", "m2 s-1");
        his_file->add_variable(his_psi_name, "", "Psi", "m2 s-1");
        his_file->add_variable(his_peclet_name, "", "Peclet (u.dx/nu)", "-");
    }

    std::string his_newton_iter_name("newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "", "Newton iterations", "-");
    std::string his_lin_solv_iter_name("linenar_solver_iterations");
    his_file->add_variable_without_location(his_lin_solv_iter_name, "iterations", "LinSolver iterations", "-");

    // Put data on time history file
    START_TIMER(Writing his-file);
    int nst_his = 0;
    his_file->put_time(nst_his, time);
    std::vector<double> his_values;

    status = set_his_values(input_data.obs_points, un, his_values);
    his_file->put_variable(his_u_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, cfl, his_values);
    his_file->put_variable(his_cfl_name, nst_his, his_values);
    if (do_viscosity) 
    {
        status = set_his_values(input_data.obs_points, pe, his_values);
        his_file->put_variable(his_peclet_name, nst_his, his_values);
        status = set_his_values(input_data.obs_points, visc, his_values);
        his_file->put_variable(his_visc_name, nst_his, his_values);
        status = set_his_values(input_data.obs_points, psi, his_values);
        his_file->put_variable(his_psi_name, nst_his, his_values);
    }

    his_values.clear();
    his_values = { 0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    his_values.clear();
    his_values = { 0.0 };
    his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);

    STOP_TIMER(Writing his-file);
    // End define history file
    ////////////////////////////////////////////////////////////////////////////

    if (total_time_steps <= 1)
    {
        std::cout << "No time loop performed, due to total_time_steps <= 1" << std::endl;
        std::cout << "Press Enter to finish";
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
    }
    Eigen::SparseMatrix<double> A(nx, nx);
    if (logging == "pattern")
    {
        std::string header_text = "=== Matrix build matrix pattern =======================";
        print_matrix_pattern(A, 2, nx, 1, header_text, log_file);
    }

    STOP_TIMER(Initialization);

    for (size_t i = 0; i < nx; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
    }

    // start time loop
    std::cout << "Start time-loop" << std::endl;
    if (stationary) { std::cout << "Stationary solution" << std::endl; }
    else { std::cout << "Time dependent simulation" << std::endl; }
    std::cout << std::fixed << std::setprecision(3) << "tstart= " << tstart + time << ";   tstop= " << tstart + tstop << ";   dt= " << dt << ";   dx= " << dx << std::endl;
 
   // Start time loop
    double du_max = 0.0;
    size_t du_maxi = 0;

    for (int i = 0; i < nx; i++)
    {
        up[i] = un[i];  // place the old solution in the new one, needed for the iteration procedure
    }

    for (int i = 0; i < nx; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
        solution[i] = 0.0;
    }


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
        std::vector<double> bc(2, 0.0);
        boundary_condition(bc[BC_WEST ], bc_vals[BC_WEST ], time, treg, bc_signals[BC_WEST], u_initial);
        boundary_condition(bc[BC_EAST ], bc_vals[BC_EAST ], time, treg, bc_signals[BC_EAST], u_initial);

        if (regularization_time)
        {
            START_TIMER(Regularization_time_loop);
            if (do_viscosity)
            {
                if (time == 7200){
                    int a = 1;
                }
                regularization->artificial_viscosity(psi, un, c_psi, dx, visc_reg);
                for (int i = 0; i < nx; ++i)
                {
                    visc[i] = visc_reg[i] + std::abs(psi[i]);
                }
            }
            STOP_TIMER(Regularization_time_loop);
        }

        START_TIMER(Newton iteration);
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
        for (int iter = 0; iter < iter_max; ++iter)
        {
            used_newton_iter += 1;
            if (logging == "iteration")
            {
                log_file << "Iteration: " << used_newton_iter << std::endl;
            }
            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                if (do_viscosity)
                {
                    regularization->artificial_viscosity(psi, up, c_psi, dx, visc_reg);
                    for (int i = 0; i < nx; ++i)
                    {
                        visc[i] = visc_reg[i] + std::abs(psi[i]);
                    }
                }
                STOP_TIMER(Regularization_iter_loop);
            }

            for (size_t k = 0; k < nx; ++k)
            {
                utheta[k] = theta * up[k] + (1.0 - theta) * un[k];
            }

            if (nst == 1 && iter == 0)
            {
                START_TIMER(Matrix initialization);
            }
            //
            // interior nodes
            //
            for (int i = 1; i < nx - 1; ++i)
            {
                double un_i   = un[i];
                double un_im1 = un[i - 1];
                double un_ip1 = un[i + 1];
                double up_i   = up[i];
                double up_im1 = up[i - 1];
                double up_ip1 = up[i + 1];

                double un_im12 = 0.5 * (un[i - 1] + un[i]);
                double un_ip12 = 0.5 * (un[i] + un[i + 1]);

                double up_im12 = 0.5 * (up[i - 1] + up[i]);
                double up_ip12 = 0.5 * (up[i] + up[i + 1]);

                double utheta_im12 = 0.5 * ( utheta[i - 1] + utheta[i] );
                double utheta_ip12 = 0.5 * ( utheta[i] + utheta[i + 1] );

                // dx * d(u)/dt
                A.coeffRef(i, i - 1) = dx * dtinv * mass[0];
                A.coeffRef(i, i) = dtinv_pseu + dx * dtinv * mass[1];
                A.coeffRef(i, i + 1) = dx * dtinv * mass[2];
                tmp[i] = -(
                      dx * dtinv * mass[0] * (up[i - 1] - un[i - 1])
                    + dx * dtinv * mass[1] * (up[i    ] - un[i    ])
                    + dx * dtinv * mass[2] * (up[i + 1] - un[i + 1])
                    );
                rhs[i] = tmp[i];
                //
                // 0.5 * d(uu)/dx = 0.5 * uu{i+1/2} - 0.5 uu_{i-1/2}
                utheta_im12 = 0.5 * (utheta[i - 1] + utheta[i]);
                utheta_ip12 = 0.5 * (utheta[i] + utheta[i + 1]);

                A.coeffRef(i, i - 1) += - 0.5 * utheta_im12 * theta;
                A.coeffRef(i, i)     += - 0.5 * utheta_im12 * theta;
                A.coeffRef(i, i)     += + 0.5 * utheta_ip12 * theta;
                A.coeffRef(i, i + 1) += + 0.5 * utheta_ip12 * theta;
                tmp[i] += -(
                    + 0.5 * utheta_ip12 * utheta_ip12 - 0.5 * utheta_im12 * utheta_im12
                    );
                rhs[i] = tmp[i];
                //
                //  d(visc(du/dx))/dx
                if (do_viscosity)
                {
                    double visc_im12 = - 0.5 * (visc[i - 1] + visc[i]);
                    double visc_ip12 = - 0.5 * (visc[i] + visc[i + 1]);

                    A.coeffRef(i, i - 1) += + visc_im12 * theta * dxinv;
                    A.coeffRef(i, i    ) += - visc_im12 * theta * dxinv;
                    A.coeffRef(i, i    ) += - visc_ip12 * theta * dxinv;
                    A.coeffRef(i, i + 1) += + visc_ip12 * theta * dxinv;
                    rhs[i] += -(
                        visc_ip12 * dxinv * (utheta[i + 1] - utheta[i]) - visc_im12 * dxinv * (utheta[i] - utheta[i - 1])
                        );
                }
                //
                // added source term
                // 
                if (do_source && src_type == "quadratic")
                {
                    A.coeffRef(i, i - 1) += theta * 2.* utheta[i-1] * 1./4.;
                    A.coeffRef(i, i    ) += theta * 2.* utheta[i  ] * 3./4.;
                    A.coeffRef(i, i    ) += theta * 2.* utheta[i  ] * 3./4.;
                    A.coeffRef(i, i + 1) += theta * 2.* utheta[i+1] * 1./4.;
                    double u_0 = 0.25 * (utheta[i-1] + 3. * utheta[i]);
                    double u_1 = 0.25 * (utheta[i+1] + 3. * utheta[i]);
                    rhs[i] += (
                        0.5 * dx * (u_0 * u_0 + u_1 * u_1)
                        );
                }
            }
            {
                //
                // wwest boundary (u>0; essential boundary)
                //
                int i = 0;

                double un_i   = un[i];           // = u^{n}_{i}
                double un_ip1 = un[i + 1];       // = u^{n}_{i+1}
                double un_ip2 = un[i + 2];       // = u^{n}_{i+2}

                double up_i   = up[i];           // = u^{n+1,p}_{i}
                double up_ip1 = up[i + 1];       // = u^{n+1,p}_{i+1}
                double up_ip2 = up[i + 2];       // = u^{n+1,p}_{i+2}

                double corr_term = 0.0;
                if (bc_type[BC_WEST] == "dirichlet")
                {
                    double un_b = w_ess[0] * un_i + w_ess[1] * un_ip1 + w_ess[2] * un_ip2;
                    double up_b = w_ess[0] * up_i + w_ess[1] * up_ip1 + w_ess[2] * up_ip2;
                    double dudt = dtinv * (up_b - un_b);
                    double u_bnd = w_ess[0] * up_i + w_ess[1] * up_ip1 + w_ess[2] * up_ip2;
                    A.coeffRef(i, i    ) = -dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0];
                    A.coeffRef(i, i + 1) = -dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1];
                    A.coeffRef(i, i + 2) = -dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2];
                    corr_term = dudt - eps_bc_corr * (bc[BC_WEST] - u_bnd);
                    rhs[i] = corr_term;
                }
                else
                {
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "West boundary: boundary type \"" << bc_type[BC_WEST] << "\" not yet supported" << std::endl;
                    std::cout << "Press Enter to finish";
                    std::cin.ignore();
                    exit(1);
                }
            }
            {
                //
                // eeast boundary (u>0; natural boundary)
                //
                int i = nx - 1;

                double un_i   = un[i];             // = u^{n}_{i}
                double un_im1 = un[i - 1];         // = u^{n}_{i-1}
                double un_im2 = un[i - 2];         // = u^{n}_{i-2}

                double up_i   = up[i];             // = u^{n+1,p}_{i}
                double up_im1 = up[i - 1];         // = u^{n+1,p}_{i-1}
                double up_im2 = up[i - 2];         // = u^{n+1,p}_{i-2}

                double utheta_i = utheta[i];
                double utheta_im1 = utheta[i - 1];
                double utheta_im2 = utheta[i - 2];

                double corr_term = 0.0;
                if (bc_type[BC_EAST] == "dirichlet")
                {
                    double un_b = w_nat[0] * un_i + w_nat[1] * un_im1 + w_nat[2] * un_im2;
                    double up_b = w_nat[0] * up_i + w_nat[1] * up_im1 + w_nat[2] * up_im2;

                    double dhdt = dtinv * (up_b - un_b);
                    A.coeffRef(i, i    ) = -dtinv * w_nat[0] - eps_bc_corr * theta * w_nat[0];
                    A.coeffRef(i, i - 1) = -dtinv * w_nat[1] - eps_bc_corr * theta * w_nat[1];
                    A.coeffRef(i, i - 2) = -dtinv * w_nat[2] - eps_bc_corr * theta * w_nat[2];
                    corr_term = dhdt - eps_bc_corr * (bc[BC_EAST] - up_b);
                    rhs[i] = corr_term;
                }
                else if (bc_type[BC_EAST] == "borsboom")
                {
                    // Outflow boundary (natural boundary)
                    //    double un_b = w_nat[0] * un_i + w_nat[1] * un_im1 + w_nat[2] * un_im2;
                    //    double up_b = w_nat[0] * up_i + w_nat[1] * up_im1 + w_nat[2] * up_im2;
                    //    double dudt = dtinv * (up_b - un_b);
                    //    
                    //    double up_im12 = 0.5 * (up[i] + up[i - 1]);
                    //    double du2dx = 0.5 * (up[i] + up[i - 1])* (up[i] - up[i - 1]) * dxinv;
                    //    double d2udx2 = 0.0;
                    //    
                    //    A.coeffRef(i, i    ) = dtinv * w_nat[0] + theta * dxinv * up_im12;
                    //    A.coeffRef(i, i - 1) = dtinv * w_nat[1] - theta * dxinv * up_im12;
                    //    A.coeffRef(i, i - 2) = dtinv * w_nat[2];
                    //    
                    //    rhs[i] = - (dudt + du2dx + d2udx2);
 
             
                
                    // Outflow boundary (natural boundary)
                    double dudt = dtinv * (
                        w_nat[0] * (up_i   - un_i) +
                        w_nat[1] * (up_im1 - un_im1) +
                        w_nat[2] * (up_im2 - un_im2)
                        );

                    double utheta_im12 = 0.5 * (utheta[i] + utheta[i - 1]);
                    double ududx =  utheta_im12 * (utheta_i - utheta_im1) * dxinv;

                    A.coeffRef(i, i    ) = dtinv * w_nat[0] + theta * dxinv * utheta_im12;
                    A.coeffRef(i, i - 1) = dtinv * w_nat[1] - theta * dxinv * utheta_im12;
                    A.coeffRef(i, i - 2) = dtinv * w_nat[2];
                    rhs[i] = - (dudt + ududx);
                
                    // viscosity contribution
                    if (do_viscosity)
                    {
                        double visc_im12 = -0.5 * (visc[i - 1] + visc[i-1]);
                        A.coeffRef(i, i - 1) += + visc_im12 * theta * dxinv;
                        A.coeffRef(i, i    ) += - visc_im12 * theta * dxinv;
                        rhs[i] += -(
                            visc_im12 * dxinv * dxinv * (utheta[i] - 2. * utheta[i - 1] + utheta[i - 2])
                            );
                    }
                }
                else
                {
                    std::cout << "----------------------------" << std::endl;
                    std::cout << "East boundary: boundary type \"" << bc_type[BC_EAST] << "\" not yet supported" << std::endl;
                    std::cout << "Press Enter to finish";
                    std::cin.ignore();
                    exit(1);
                }
            }
            STOP_TIMER(Matrix set up);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(Matrix initialization);
            }
            if (nst == 1 && iter == 0)
            {
                START_TIMER(BiCGstab_initialization);
            }
            else
            {
                START_TIMER(BiCGstab);
            }

            Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
            solver.compute(A);
            solver.setTolerance(eps_bicgstab);
            //solution = solver.solve(rhs);
            solution = solver.solveWithGuess(rhs, solution);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(BiCGstab_initialization);
            }
            else
            {
                STOP_TIMER(BiCGstab);
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
            for (size_t i = 0; i < nx; ++i)
            {
                up[i] += solution[i];
                delta_u[i] = solution[i]; 
                if ( du_max < std::abs(delta_u[i]) )
                {
                    du_max = std::abs(delta_u[i]);
                    du_maxi = i;
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
                log_file << "=== cn, delta_c, up ===================================" << std::endl;
                for (int i = 0; i < nx; ++i)
                {
                    log_file << std::setw(17) << std::setprecision(15) << std::scientific 
                        << un[i] << "  -----  " 
                        << delta_u[i] << "  -----  " 
                        << up[i] << std::endl;

                }
            }
            if (logging == "eigen_values")
            {
                log_file << "=== Eigen values ======================================" << std::endl;
                log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).eigenvalues() << std::endl;
            }
            if (logging == "matrix")
            {
                log_file << "====== End iteration ==================================" << std::endl;
                log_file << "time [sec]: " << std::setw(8) << std::setprecision(3) << time
                         << std::setprecision(8) << std::scientific
                         << ";    iterations     : " << solver.iterations()
                         << ";    estimated error: " << solver.error() 
                         << std::endl;
            }

            used_lin_solv_iter = std::max(used_lin_solv_iter, (int)solver.iterations());
            if (du_max < eps_newton)
            {
                break;
            }
        }
        STOP_TIMER(Newton iteration);

        if (stationary)
        {
            log_file << "stationary solution " << std::endl;
            log_file << std::setprecision(8) << std::scientific
                << "    Newton iterations  : " << used_newton_iter << std::endl
                << "    Delta u^{n + 1,p + 1}: " << du_max << " at: " << du_maxi
                << std::endl;
        }
        else
        {
            if (std::fmod(time, input_data.output.dt_screen) == 0)
            {
                std::cout << std::fixed << std::setprecision(2) << tstart + time << ";   " << tstart + tstop << std::endl;
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(2) << std::scientific << time
                    << std::setprecision(8) << std::scientific
                    << "    Newton iterations  : " << used_newton_iter << std::endl
                    << "    Delta u^{n + 1,p + 1}: " << du_max << " at: " << du_maxi
                    << std::endl;
            }
        }
        if (used_newton_iter == iter_max)
        {
            if (du_max > eps_newton)
            {
                log_file.setf(std::ios::fixed, std::ios::floatfield);
                log_file << "    ----    maximum number of iterations reached, probably not converged, at time: " <<  time << " [sec]" << std::endl;

            }
        }
        for (size_t i = 0; i < nx; ++i)
        {
            un[i] = up[i];
        }

        for (int i = 0; i < nx; ++i)
        {
            cfl[i] = up[i] * dt / dx;
        }
        if (do_viscosity)
        {
            for (int i = 0; i < nx; ++i)
            {
                pe[i] = up[i] * dx / visc[i];
            }
        }
        if (logging == "matrix")
        {
            log_file << "=== up, cn, up-cn =====================================" << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                log_file << std::setprecision(15) << std::scientific 
                    << up[i] << ",   " 
                    << un[i] << ",   " 
                    << up[i] - un[i] 
                    << std::endl;
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
                map_file->put_time(nst_map, double(nst) * dt);
            }
            map_file->put_time_variable(map_names[0], nst_map, un);
            map_file->put_time_variable(map_names[4], nst_map, cfl);
            if (do_viscosity)
            {
                map_file->put_time_variable(map_names[1], nst_map, visc);
                map_file->put_time_variable(map_names[2], nst_map, psi);
                map_file->put_time_variable(map_names[3], nst_map, pe);
            }
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
            status = set_his_values(input_data.obs_points, un, his_values);
            his_file->put_variable(his_u_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, cfl, his_values);
            his_file->put_variable(his_cfl_name, nst_his, his_values);
            if (do_viscosity)
            {
                status = set_his_values(input_data.obs_points, pe, his_values);
                his_file->put_variable(his_peclet_name, nst_his, his_values);
                status = set_his_values(input_data.obs_points, visc, his_values);
                his_file->put_variable(his_visc_name, nst_his, his_values);
                status = set_his_values(input_data.obs_points, psi, his_values);
                his_file->put_variable(his_psi_name, nst_his, his_values);
            }

            his_values.clear();
            his_values = { double(used_newton_iter) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            his_values.clear();
            his_values = { double(used_lin_solv_iter) };
            his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);
            STOP_TIMER(Writing his-file);
        }
    } // End of the time loop
    STOP_TIMER(Time loop);

    // Finalization
    log_file.close();
    status = his_file->close();
    status = map_file->close();

    STOP_TIMER(Main);
    PRINT_TIMER(timing_filename.data());

    std::cout << "--- Klaar ---" << std::endl;
    std::chrono::duration<int, std::milli> timespan(1000);
    std::this_thread::sleep_for(timespan);
    return 0;
}
//------------------------------------------------------------------------------
int set_his_values(std::vector<_ObservationPoint>& obs_points, std::vector<double> & array, std::vector<double>& his_values)
{
    his_values.clear();
    for (size_t i = 0; i < obs_points.size(); ++i)
    {
        his_values.push_back(array[obs_points[i].idx]);
    }
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
//------------------------------------------------------------------------------
/* @@ GetArguments
*
* Scan command-line arguments
*/
void GetArguments(long argc,   /* I Number of command line arguments */
    char** argv,
    std::filesystem::path & file_name)
    /* Returns nothing */
{
    long  i;
    i = 0;
    while (i < argc)
    {
        if (strcmp(argv[i], "--toml") == 0)
        {
            i = i + 1;
            if (i < argc) file_name = std::string(argv[i]);
        }
        else if (i != 0)
        {
            i = argc;
        }
        i = i + 1;
    }
    return;
}
