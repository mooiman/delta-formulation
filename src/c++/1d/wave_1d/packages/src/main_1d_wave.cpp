//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D wave equation, fully implicit with delta-formulation and Modified Newton iteration 
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
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <thread>
#include <format>
#include <string_view>

#include <include/KDtree.hpp>
#include <toml.h>
#include <include/KDtree.hpp>
// for BiCGstab  solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "bed_level.h"
#include "bed_shear_stress.h"
#include "boundary_condition.h"
#include "cfts.h"
#include "compile_date_and_time.h"
#include "data_input_struct.h"
#include "definition_map_file.h"
#include "initial_conditions.h"
#include "main_version.h"
#include "observation_stations.h"
#include "perf_timer.h"
#include "print_matrix.h"
#include "read_input_toml_file.h"
#include "regularization.h"
#include "ugrid1d.h"

double Fabs(double, double);
void GetArguments(long argc, char** argv, std::filesystem::path & file_name);
inline size_t p_index(size_t i, size_t j, size_t nx);
int write_used_input(struct _data_input data, std::ofstream & log_file);
int set_his_values(std::vector<_ObservationPoint>& obs_points, std::vector<double> & array, std::vector<double>& his_values);

int read_bed_level(std::string filename, std::vector<double> & value);
int initialize_scalar(double, std::vector<double>&, std::vector<double>&);
double scv(double, double);

// Solve the linear wave equation
// Continuity equation: d(h)/dt + d(q)/dx = 0
// Momentum equation  : d(q)/dt + gh d(zeta)/dx + convection + bed_shear_stress - diffusivity = 0

int main(int argc, char* argv[])
{
    bool stationary = false;
    double sign = 1.0;
    int solver_iterations = -1;
    std::filesystem::path toml_file_name("---not-defined---");

    std::map<BED_LEVEL_ENUM, std::string> bed_level_name;
    bed_level_name[BED_LEVEL_ENUM::FLAT] = "flat";
    bed_level_name[BED_LEVEL_ENUM::SHOAL] = "shoal";
    bed_level_name[BED_LEVEL_ENUM::SLOPED] = "sloped";
    bed_level_name[BED_LEVEL_ENUM::WAVY] = "wavy";
    bed_level_name[BED_LEVEL_ENUM::WAVY_SLOPED] = "wavy_sloped";
    bed_level_name[BED_LEVEL_ENUM::WEIR] = "weir";

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
    BED_LEVEL_ENUM bed_level_type = BED_LEVEL_ENUM::NONE;

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

    if (input_data.numerics.dt == 0.0) { stationary = true;  }
    if (stationary) { 
        input_data.boundary.treg = 0.0; 
        input_data.numerics.theta = 1.0;
    }
    if (input_data.initial.gauss_mu != -INFINITY)
    {
        input_data.initial.gauss_mu_x = input_data.initial.gauss_mu; 
    }
    if (input_data.initial.gauss_sigma != -INFINITY)
    {
        input_data.initial.gauss_sigma_x = input_data.initial.gauss_sigma; 
    }

    ss << "_dx" << input_data.numerics.dx << "_dt" << input_data.numerics.dt;

    //
    // Bed level
    //
    std::string geometry_type = input_data.domain.geometry_type;
    BED_LEVEL * bed = new BED_LEVEL();
    status = bed->set_bed_level_type(geometry_type, bed_level_type);
    if (bed_level_type == BED_LEVEL_ENUM::FLAT)
    {
        std::stringstream depth_strm;
        depth_strm << std::fixed << std::setprecision(0) << input_data.domain.depth;
        bed_level_name[BED_LEVEL_ENUM::FLAT] = "flat" + depth_strm.str();
    }

    std::string inp_bed = "bed_level_" + bed_level_name[bed_level_type];

    out_file = output_dir.string() + inp_bed + ss.str();
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

    double dxinv = 1./input_data.numerics.dx;                                 // invers grid size [m]
    size_t nx = int(input_data.domain.Lx * dxinv) + 1 + 2 ;                    // number of nodes; including 2 virtual points
    size_t ny = 1; 
    size_t nxny = nx * ny; 

    int total_time_steps = int((input_data.time.tstop - input_data.time.tstart) / input_data.numerics.dt) + 1;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    int wrihis;                                           // write interval to his-file
    int wrimap;                                           // write interval to map-file
    if (stationary)
    {
        input_data.numerics.dt = 0.0;                                         // Time step size [s]
        dtinv = 0.0;                                      // stationary solution
        input_data.boundary.eps_bc_corr = 1.0;
        input_data.numerics.iter_max = 2 * input_data.numerics.iter_max;
        input_data.numerics.theta = 1.0;                                      // Stationary solution
        input_data.time.tstop = 1.;
        total_time_steps = 2;                             // initial step (step 1), stationary result (step 2)
        input_data.boundary.treg = 0.0;                                       // Thatcher-Harleman return time [s], when zero supply boundary value immediately
        wrihis = 1;                                       // write interval to his-file
        wrimap = 1;                                       // write interval to map-file
    }
    else
    {
        dtinv = 1. / input_data.numerics.dt;                                  // Inverse of dt [1/s]
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

    std::string model_title("Linear wave equation");

    if (input_data.physics.do_q_equation)
    {
        model_title = "Linear wave equation";
        if (input_data.physics.do_bed_shear_stress)
        {
            model_title += " + bed shear stress";
        }
        if (input_data.physics.do_convection)
        {
            model_title += " + convection";
        }
        if (input_data.physics.do_viscosity)
        {
            model_title += " + viscosity";
        }
    }
    else if (input_data.physics.do_q_equation)
    {
        model_title = "Linear wave equation, only q-equation";
    }
    else
    {
        model_title = "No q-equation (=> no waves computed)";
    }

    model_title += ", " + solver_name;

    REGULARIZATION* regularization = new REGULARIZATION(input_data.numerics.iter_max, input_data.physics.g, input_data.log.logging);

    log_file << "=== Used input variables ==============================" << std::endl;
    status = write_used_input(input_data, log_file);
    log_file << "=======================================================" << std::endl;

    // Copy input data to loccal data
    std::string logging = input_data.log.logging;

    double Lx = input_data.domain.Lx;
    double depth = input_data.domain.depth;

    double eps_bc_corr = input_data.boundary.eps_bc_corr;
    double treg = input_data.boundary.treg;
    std::vector<std::string> bc_type = input_data.boundary.bc_type;
    std::vector<std::string> bc_vars = input_data.boundary.bc_vars;
    std::vector<double> bc_vals = input_data.boundary.bc_vals;

    double gauss_amp = input_data.initial.gauss_amp;
    double gauss_mu = input_data.initial.gauss_mu;
    double gauss_mu_x = input_data.initial.gauss_mu_x;
    double gauss_mu_y = input_data.initial.gauss_mu_y;
    double gauss_sigma = input_data.initial.gauss_sigma;
    double gauss_sigma_x = input_data.initial.gauss_sigma_x;
    double gauss_sigma_y = input_data.initial.gauss_sigma_y;
    std::vector<std::string> ini_vars = input_data.initial.ini_vars;

    double dt = input_data.numerics.dt;
    double dx = input_data.numerics.dx;
    double c_psi = input_data.numerics.c_psi;
    double eps_bicgstab = input_data.numerics.eps_bicgstab;
    double eps_newton = input_data.numerics.eps_newton;
    double eps_abs = input_data.numerics.eps_abs;
    double theta = input_data.numerics.theta;
    int iter_max = input_data.numerics.iter_max;
    std::string linear_solver = input_data.numerics.linear_solver;
    bool regularization_init = input_data.numerics.regularization_init;
    bool regularization_iter = input_data.numerics.regularization_iter;
    bool regularization_time = input_data.numerics.regularization_time;

    double g = input_data.physics.g;
    double chezy_coefficient = input_data.physics.chezy_coefficient;
    double visc_const = input_data.physics.visc_const;
    bool do_linear_waves = input_data.physics.do_linear_waves;
    bool do_bed_shear_stress = input_data.physics.do_bed_shear_stress;
    bool do_convection = input_data.physics.do_convection;
    bool do_viscosity = input_data.physics.do_viscosity;

    double tstart = input_data.time.tstart;
    double tstop = input_data.time.tstop;

    // Declare arrays
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix in x-direction

    std::vector<double> x(nx, 0.);                      // x-coordinate
    std::vector<double> y(nx, 0.);                      // y-coordinate, also needed in 1D
    std::vector<double> zb(nx, -10.0);                  // regularized bed level
    std::vector<double> zb_given(nx, -10.0);            // initial given bed level

    std::vector<double> rhs_viscosity(nx, 0.);      // rhs of the momentum equation for viscosity
    std::vector<double> tmp1(nx, 0.);                   // help array
    std::vector<double> tmp2(nx, 0.);                   // help array

    std::vector<double> qn(nx, 0.);                     // flow at (n)
    std::vector<double> hn(nx, 0.);                     // total depth at (n)
    std::vector<double> qp(nx, 0.);                     // flow at (n+1,p), previous iteration
    std::vector<double> hp(nx, 0.);                     // total depth at (n+1,p), previous iteration
    std::vector<double> delta_h(nx, 0.);
    std::vector<double> delta_q(nx, 0.);
    std::vector<double> s_given(nx, 0.);  // given initial water level
    std::vector<double> u_given(nx, 0.);  // given initial velocity
    std::vector<double> s(nx, 0.);  // smoothed water level
    std::vector<double> u(nx, 0.);  // water velocity, needed for postprocessing
    std::vector<double> riemann_pos(nx, 0.);  // Riemann invarinat going to the right, needed for postprocessing
    std::vector<double> riemann_neg(nx, 0.);  // Riemann invarinat going to the left, needed for postprocessing
    std::vector<double> visc_given(nx, visc_const);  // Initialize viscosity array with given value
    std::vector<double> visc_reg(nx, visc_const);  // Initialize given viscosity array with regularized value
    std::vector<double> visc(nx, visc_const);  // Viscosity array used for computation
    std::vector<double> pe(nx, 0.);  // peclet number [-]
    std::vector<double> psi(nx, 0.);  //
    std::vector<double> froude(nx, 0.);  //

    std::vector<double> htheta(nx, 0.);
    std::vector<double> qtheta(nx, 0.);

    Eigen::VectorXd solution(2 * nx);               // solution vector [h, q]^{n}
    Eigen::VectorXd rhs(2 * nx);                // RHS vector [h, q]^n

    solution.setZero(); // Delta h and Delta q
    rhs.setZero();

    double cf_given = g / (chezy_coefficient * chezy_coefficient);  // bed friction coefficient
    std::vector<double> cf(nx, cf_given);  // Bed shear stress coefficient

    double alpha = 1. / 8.;                                   // Linear (spatial) interpolation coefficient
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

    //initialize water level
    std::cout << "    Initialisation" << std::endl;

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i - 1) * dx - Lx / 2;
    }
    //  Create kdtree, needed to locate the observation points
    std::vector<std::vector<double>> xy_points;
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::vector<double> point = {x[i], y[i]};
        xy_points.push_back(point);
    }
    KDTree xy_tree(xy_points);

    status = 0;
    double min_zb = *std::min_element(zb_given.begin(), zb_given.end());
    double s_offset = 0.0;

    log_file << "Nodes   : " << nx << std::endl;
    log_file << "Elements: " << (nx - 1) <<  std::endl;
    log_file << "Volumes : " << (nx - 2) << std::endl;
    log_file << "=======================================================" << std::endl;
    std::cout << "    Nodes: " << nx << std::endl;
    std::cout << "======================================================" << std::endl;

    STOP_TIMER(Writing log-file);  // but two write statements are not timed
    if (status != 0) {
        log_file.close();
        exit(1);
    }

    initial_conditions(x, nx, s_given, u_given, ini_vars, 
        gauss_amp, gauss_mu_x, gauss_sigma_x);

    status = bed->initialize_bed_level(bed_level_type, x, zb_given, model_title, depth);

    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        regularization->given_function(zb, psi, zb_given, dx, c_psi);
        regularization->given_function(visc_reg, psi, visc_given, dx, c_psi);
        for (size_t i = 0; i < nx; ++i)
        {
            visc[i] = visc_reg[i] + std::abs(psi[i]);
        }
        STOP_TIMER(Regularization_init);
    }
    else
    {
        for (size_t i = 0; i < zb_given.size(); ++i)
        {
            zb[i] = zb_given[i];
        }
    }
    for (int i = 0; i < nx; ++i)
    {
        froude[i] = qp[i] / hp[i] / std::sqrt(g * hp[i]);
        pe[i] = qp[i] / hp[i] * dx / visc[i];
    }

    for (size_t k = 0; k < zb_given.size(); k++)
    {
        hn[k] = s_given[k] - zb[k];  // Initial water depth
        qn[k] = hn[k] * u_given[k];  // initial flow flux
        //
        hp[k] = hn[k];  // initial water depth
        qp[k] = qn[k];  // initial flow flux
        //
        riemann_pos[k] = std::abs(qn[k] + std::sqrt(g * hn[k]) * hn[k]);
        riemann_neg[k] = std::abs(qn[k] - std::sqrt(g * hn[k]) * hn[k]);
    }
    for (size_t i = 0; i < zb_given.size(); ++i)
    {
        s[i] = hn[i] + zb[i];
        u[i] = qn[i] / hn[i];
    }

    double time = double(0) * dt;
    ////////////////////////////////////////////////////////////////////////////
    // Create map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    UGRID1D* map_file = new UGRID1D();

    std::vector<std::string> map_names;
    map_names.push_back("hn_1d");
    map_names.push_back("qn_1d");
    map_names.push_back("s_1d");
    map_names.push_back("u_1d");
    map_names.push_back("zb_1d");
    map_names.push_back("delta_h");
    map_names.push_back("delta_q");
    map_names.push_back("visc_reg_1d");
    map_names.push_back("visc_1d");
    map_names.push_back("psi_1d");
    map_names.push_back("froude_1d");
    map_names.push_back("pe_1d");
    map_names.push_back("visc_rhs");

    map_file = create_map_file(nc_mapfilename, model_title, x, map_names, do_viscosity);

    // Put data on map file
    START_TIMER(Writing map-file);
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    map_file->put_time_variable(map_names[0], nst_map, hn);
    map_file->put_time_variable(map_names[1], nst_map, qn);
    map_file->put_time_variable(map_names[2], nst_map, s);
    map_file->put_time_variable(map_names[3], nst_map, u);
    map_file->put_time_variable(map_names[4], nst_map, zb);
    map_file->put_time_variable(map_names[5], nst_map, delta_h);
    map_file->put_time_variable(map_names[6], nst_map, delta_q);
    map_file->put_time_variable(map_names[10], nst_map, froude);
    if (do_viscosity)
    {
        map_file->put_time_variable(map_names[7], nst_map, visc_reg);
        map_file->put_time_variable(map_names[8], nst_map, visc);
        map_file->put_time_variable(map_names[9], nst_map, psi);
        map_file->put_time_variable(map_names[11], nst_map, pe);
        map_file->put_time_variable(map_names[12], nst_map, rhs_viscosity);
    }
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

    std::string his_h_name("hn_1d");
    std::string his_q_name("qn_1d");
    std::string his_s_name("water_level");
    std::string his_u_name("u_velocity");
    std::string his_zb_name("bed_level");
    std::string his_froude_name("froude");
    std::string his_peclet_name("his_peclet_name");
    std::string his_visc_name("his_visc_name");
    std::string his_psi_name("his_psi_name");
    std::string his_visc_rhs_name("his_visc_rhs_name");
    std::string his_riemann_pos_name("Riemann_right_going");
    std::string his_riemann_neg_name("Riemann_left_going");
    his_file->add_variable(his_h_name, "", "Water depth", "m");
    his_file->add_variable(his_q_name, "", "Water flux", "m2 s-1");
    his_file->add_variable(his_s_name, "", "Water level", "m");
    his_file->add_variable(his_u_name, "", "Water velocity", "m s-1");
    his_file->add_variable(his_zb_name, "", "Bed level", "m");
    his_file->add_variable(his_froude_name, "", "Froude", "-");
    if (do_viscosity)
    {
        his_file->add_variable(his_peclet_name, "", "Peclet", "-");
        his_file->add_variable(his_visc_name, "", "Viscosity (used)", "m2 s-1");
        his_file->add_variable(his_psi_name, "", "Psi", "m2 s-1");
        his_file->add_variable(his_visc_rhs_name, "", "Viscosity (rhs)", "m2 s-1");
    }

    //his_file->add_variable(his_riemann_pos_name, "", "Rieman(+)", "m2 s-1");
    //his_file->add_variable(his_riemann_neg_name, "", "Rieman(-)", "m2 s-1");

    std::string his_newton_iter_name("newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "", "Newton iterations", "-");
    std::string his_LinSolver_iter_name("his_LinSolver_iterations");
    his_file->add_variable_without_location(his_LinSolver_iter_name, "iterations", "LinSolver iterations", "-");
    std::string his_LinSolver_iter_error_name("LinSolver_iteration_error");
    his_file->add_variable_without_location(his_LinSolver_iter_error_name, "iteration_error", "LinSolver iteration error", "-");

    // Put data on time history file
    START_TIMER(Writing his-file);
    int nst_his = 0;
    his_file->put_time(nst_his, time);
    std::vector<double> his_values;

    status = set_his_values(input_data.obs_points, hn, his_values);
    his_file->put_variable(his_h_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, qn, his_values);
    his_file->put_variable(his_q_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, s, his_values);
    his_file->put_variable(his_s_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, u, his_values);
    his_file->put_variable(his_u_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, zb, his_values);
    his_file->put_variable(his_zb_name, nst_his, his_values);

    status = set_his_values(input_data.obs_points, froude, his_values);
    his_file->put_variable(his_froude_name, nst_his, his_values);

    if (do_viscosity)
    {
        status = set_his_values(input_data.obs_points, pe, his_values);
        his_file->put_variable(his_peclet_name, nst_his, his_values);
        status = set_his_values(input_data.obs_points, visc, his_values);
        his_file->put_variable(his_visc_name, nst_his, his_values);
        status = set_his_values(input_data.obs_points, psi, his_values);
        his_file->put_variable(his_psi_name, nst_his, his_values);
        status = set_his_values(input_data.obs_points, rhs_viscosity, his_values);
        his_file->put_variable(his_visc_rhs_name, nst_his, his_values);
    }

    //his_values = { riemann_pos[p_a], riemann_pos[p_b], riemann_pos[p_c], riemann_pos[p_d], riemann_pos[p_e] };
    //his_file->put_variable(his_riemann_pos_name, nst_his, his_values);
    //his_values = { riemann_neg[p_a], riemann_neg[p_b], riemann_neg[p_c], riemann_neg[p_d], riemann_neg[p_e] };
    //his_file->put_variable(his_riemann_neg_name, nst_his, his_values);

    his_values.clear();
    his_values = { 0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    his_values.clear();
    his_values = { 0.0 };
    his_file->put_variable(his_LinSolver_iter_name, nst_his, his_values);

    his_values.clear();
    his_values = { 0.0 };
    his_file->put_variable(his_LinSolver_iter_error_name, nst_his, his_values);

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
    Eigen::SparseMatrix<double> A(2 * nx, 2 * nx);
    if (logging == "pattern")
    {
        std::string header_text = "=== Matrix build matrix pattern =======================";
        print_matrix_pattern(A, 2, nx, 1, header_text, log_file);
    }

    STOP_TIMER(Initialization);

    for (int i = 0; i < 2 * nx; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
    }

    std::cout << "Start time-loop" << std::endl;
    if (stationary) { std::cout << "Stationary solution" << std::endl; }
    else { std::cout << "Time dependent simulation" << std::endl; }
    std::cout << std::fixed << std::setprecision(3) << "tstart= " << tstart + time << ";   tstop= " << tstart + tstop << ";   dt= " << dt << std::endl;
 
   // Start time loop

    START_TIMER(Time loop);
    double dh_max = 0.0;
    double dq_max = 0.0;
    size_t dh_maxi = 0;
    size_t dq_maxi = 0;

    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        time = dt * double(nst);

        int select = 1;  // Constant boundary condition
        std::vector<double> bc(2, 0.0);
        boundary_condition(bc[BC_EAST ], bc_vals[BC_EAST ], time, treg, select);
        boundary_condition(bc[BC_WEST ], bc_vals[BC_WEST ], time, treg, select);

        // Set the right-hand side (rhs) of the equations (at t = nst*dt)
        //std::cout << std::setprecision(5) << "Time: " << dt*double(nst) << std::endl;

        double hn_im2, hn_im1, hn_0, hn_ip1, hn_ip2;  // previous timestep (n)
        double hn_im12, hn_ip12;
        double hp_im2, hp_im1, hp_0, hp_ip1, hp_ip2;  // previous iteration (p)
        double hp_im12, hp_ip12;

        double qn_im2, qn_im1, qn_0, qn_ip1, qn_ip2;  // previous timestep (n)
        double qn_im12, qn_ip12;
        double qp_im2, qp_im1, qp_0, qp_ip1, qp_ip2;  // previous iteration (p)
        double qp_im12, qp_ip12;

        double htheta_0;
        double htheta_im1, htheta_im2, htheta_ip1, htheta_ip2;
        double htheta_im12, htheta_ip12;
        double qtheta_0;
        double qtheta_im1, qtheta_im2, qtheta_ip1, qtheta_ip2;
        double qtheta_im12, qtheta_ip12;

        int used_newton_iter = 0;
        int used_lin_solv_iter = 0;
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
                    regularization->artificial_viscosity(psi, hp, qp, zb, c_psi, dx);
                    for (int i = 0; i < nx; ++i)
                    {
                        visc[i] = visc_reg[i] + std::abs(psi[i]);
                        pe[i] = qp[i] / hp[i] * dx / visc[i];
                    }
                }
                STOP_TIMER(Regularization_iter_loop);
            }

            for (size_t k = 0; k < nx; ++k)
            {
                htheta[k] = theta * hp[k] + (1.0 - theta) * hn[k];
                qtheta[k] = theta * qp[k] + (1.0 - theta) * qn[k];
            }

            if (nst == 1 && iter == 0)
            {
                START_TIMER(Matrix initialization);
            }
            //
            // interior nodes
            //
            for (int i = 1; i < nx - 1; i++)
            {
                hn_im1 = hn[i - 1]; // = h^{n}_{i-1}
                hn_0 = hn[i];       // = h^{n}_{i}
                hn_ip1 = hn[i + 1]; // = h^{n}_{i+1}

                qn_im1 = qn[i - 1]; // = q^{n}_{i-1}
                qn_0 = qn[i];       // = q^{n}_{i}
                qn_ip1 = qn[i + 1]; // = q^{n}_{i+1}

                hn_im12 = 0.5 * (hn_0 + hn_im1);
                hn_ip12 = 0.5 * (hn_ip1 + hn_0);

                qn_im12 = 0.5 * (qn_0 + qn_im1);
                qn_ip12 = 0.5 * (qn_ip1 + qn_0);

                hp_im1 = hp[i - 1]; // = h^{n+1,p}_{i-1}
                hp_0 = hp[i];       // = h^{n+1,p}_{i}
                hp_ip1 = hp[i + 1]; // = h^{n+1,p}_{i+1}

                qp_im1 = qp[i - 1]; // = q^{n+1,p}_{i-1}
                qp_0 = qp[i];  // = q^{n+1,p}_{i}
                qp_ip1 = qp[i + 1]; // = q^{n+1,p}_{i+1}

                hp_im12 = 0.5 * (hp_0 + hp_im1);  // = h^{n+1,p}_{i-1/2}
                hp_ip12 = 0.5 * (hp_ip1 + hp_0);  // = h^{n+1,p}_{i+1/2}

                qp_im12 = 0.5 * (qp_0 + qp_im1);  // = q^{n+1,p}_{i-1/2}
                qp_ip12 = 0.5 * (qp_ip1 + qp_0);  // = q^{n+1,p}_{i+1/2}

                htheta_0   = htheta[i];
                htheta_im1 = htheta[i - 1];
                htheta_ip1 = htheta[i + 1];
                htheta_im12 = 0.5 * (htheta[i] + htheta[i - 1]);
                htheta_ip12 = 0.5 * (htheta[i] + htheta[i + 1]);

                qtheta_0   = qtheta[i];
                qtheta_im1 = qtheta[i - 1];
                qtheta_ip1 = qtheta[i + 1];
                qtheta_im12 = 0.5 * (qtheta[i] + qtheta[i - 1]);
                qtheta_ip12 = 0.5 * (qtheta[i] + qtheta[i + 1]);

                int ph   = 2 * p_index(i    , 0, nx);  // continuity equation
                int ph_e = 2 * p_index(i + 1, 0, nx);  // continuity equation
                int ph_w = 2 * p_index(i - 1, 0, nx);  // continuity equation

                int c_eq = ph;
                int q_eq = c_eq + 1;
                //
                // continuity equation (dh/dt ... = 0)
                //
                rhs[c_eq] = 0.0;
                // mass-matrix continuity
                A.coeffRef(c_eq, ph_w) = dtinv * dx * mass[0];
                A.coeffRef(c_eq, ph) = dtinv * dx * mass[1];
                A.coeffRef(c_eq, ph_e) = dtinv * dx * mass[2];
                // flow flux
                A.coeffRef(c_eq, ph_w + 1) = -0.5 * theta;
                A.coeffRef(c_eq, ph + 1) = 0.5 * theta - 0.5 * theta;
                A.coeffRef(c_eq, ph_e + 1) = 0.5 * theta;
                //
                rhs[c_eq] -=
                    dtinv * dx * mass[0] * (hp_im1 - hn_im1) +
                    dtinv * dx * mass[1] * (hp_0 - hn_0) +
                    dtinv * dx * mass[2] * (hp_ip1 - hn_ip1);
                // flux
                rhs[c_eq] += -(qtheta_ip12 - qtheta_im12);
                //
                // momentum (dq/dt ... = 0)
                //
                rhs[q_eq] = 0.0;
                // mass-matrix momentum
                A.coeffRef(q_eq, ph_w + 1) = dtinv * dx * mass[0];
                A.coeffRef(q_eq, ph + 1) = dtinv * dx * mass[1];
                A.coeffRef(q_eq, ph_e + 1) = dtinv * dx * mass[2];
                // 
                rhs[q_eq] -=
                    dtinv * dx * mass[0] * (qp_im1 - qn_im1) +
                    dtinv * dx * mass[1] * (qp_0 - qn_0) +
                    dtinv * dx * mass[2] * (qp_ip1 - qn_ip1);
                //
                // pressure
                //
                double ztheta_im1 = htheta_im1 + zb[i - 1];
                double ztheta_0 = htheta_0 + zb[i];
                double ztheta_ip1 = htheta_ip1 + zb[i + 1];
                // Contribution due to Delta h^{n+1, p+1}
                A.coeffRef(q_eq, ph_w) = 0.125 * theta * g * (ztheta_0 - ztheta_im1);
                A.coeffRef(q_eq, ph) = 0.375 * theta * (ztheta_0 - ztheta_im1) + 0.375 * theta * g * (ztheta_ip1 - ztheta_0);
                A.coeffRef(q_eq, ph_e) = 0.125 * theta * g * (ztheta_ip1 - ztheta_0);
                // contribution due to Delta zeta^{n+1, p+1} 
                A.coeffRef(q_eq, ph_w) += -0.125 * theta * g * (htheta_im1 + 3. * htheta_0);
                A.coeffRef(q_eq, ph) += 0.125 * theta * g * (htheta_im1 + 3. * htheta_0) - 0.125 * theta * g * (htheta_ip1 + 3. * htheta_0);
                A.coeffRef(q_eq, ph_e) += 0.125 * theta * g * (htheta_ip1 + 3. * htheta_0);

                double rhs_pressure = -0.5 * g * 0.25 * (
                    (1.0 * htheta_im1 + 3.0 * htheta_0) * (ztheta_0 - ztheta_im1) +
                    (3.0 * htheta_0 + 1.0 * htheta_ip1) * (ztheta_ip1 - ztheta_0)
                    );
                rhs[q_eq] += rhs_pressure;
                //
                if (do_convection)
                {
                    //
                    // convection
                    //
                    A.coeffRef(q_eq, ph_w) += +0.5 * theta * qtheta_im12 * qtheta_im12 / (htheta_im12 * htheta_im12);
                    A.coeffRef(q_eq, ph) += -0.5 * theta * qtheta_ip12 * qtheta_ip12 / (htheta_ip12 * htheta_ip12) + 0.5 * qtheta_im12 * qtheta_im12 / (htheta_im12 * htheta_im12);
                    A.coeffRef(q_eq, ph_e) += -0.5 * theta * qtheta_ip12 * qtheta_ip12 / (htheta_ip12 * htheta_ip12);
                    //
                    A.coeffRef(q_eq, ph_w + 1) += -theta * qtheta_im12 / htheta_im12;
                    A.coeffRef(q_eq, ph + 1) += -theta * qtheta_im12 / htheta_im12 + theta * qtheta_ip12 / htheta_ip12;
                    A.coeffRef(q_eq, ph_e + 1) += +theta * qtheta_ip12 / htheta_ip12;
                    //
                    double rhs_convection = -(qtheta_ip12 * qtheta_ip12 / htheta_ip12 - qtheta_im12 * qtheta_im12 / htheta_im12);
                    rhs[q_eq] += rhs_convection;
                }
                if (do_viscosity)
                {
                    //
                    // viscosity
                    // 
                    double visc_im12 = -0.5 * (visc[i - 1] + visc[i]);
                    double visc_ip12 = -0.5 * (visc[i] + visc[i + 1]);

                    double A_im12 = visc_im12 * theta;
                    double A_ip12 = visc_ip12 * theta;
                    double B_im12 = -visc_im12 * theta / (htheta_im12) * dxinv * (htheta_0 - htheta_im1);
                    double B_ip12 = -visc_ip12 * theta / (htheta_ip12) * dxinv * (htheta_ip1 - htheta_0);
                    double C_im12 = visc_im12 * theta * qtheta_im12 / (htheta_im12 * htheta_im12) * dxinv * (htheta_0 - htheta_im1);
                    double C_ip12 = visc_ip12 * theta * qtheta_ip12 / (htheta_ip12 * htheta_ip12) * dxinv * (htheta_ip1 - htheta_0);
                    double D_im12 = -visc_im12 * theta * qtheta_im12 / htheta_im12;
                    double D_ip12 = -visc_ip12 * theta * qtheta_ip12 / htheta_ip12;
                    // 
                    // nu theta d(dzeta)/dx
                    A.coeffRef(q_eq, ph_w + 1) +=  dxinv * A_im12;
                    A.coeffRef(q_eq, ph + 1)   += -dxinv * A_im12;
                    A.coeffRef(q_eq, ph + 1)   += -dxinv * A_ip12;
                    A.coeffRef(q_eq, ph_e + 1) +=  dxinv * A_ip12;
                    //
                    // - nu theta q/h dq
                    A.coeffRef(q_eq, ph_w + 1) += -0.5 * B_im12;
                    A.coeffRef(q_eq, ph + 1)   += -0.5 * B_im12;
                    A.coeffRef(q_eq, ph + 1)   +=  0.5 * B_ip12;
                    A.coeffRef(q_eq, ph_e + 1) +=  0.5 * B_ip12;
                    //
                    // nu theta q/h^2 d(h)/dx dh
                    A.coeffRef(q_eq, ph_w) += -0.5 * C_im12;
                    A.coeffRef(q_eq, ph)   += -0.5 * C_im12;
                    A.coeffRef(q_eq, ph)   +=  0.5 * C_ip12;
                    A.coeffRef(q_eq, ph_e) +=  0.5 * C_ip12;
                    //
                    // - nu theta q/h d(dh)/dx
                    A.coeffRef(q_eq, ph_w) +=  dxinv * D_im12;
                    A.coeffRef(q_eq, ph)   += -dxinv * D_im12;
                    A.coeffRef(q_eq, ph)   += -dxinv * D_ip12;
                    A.coeffRef(q_eq, ph_e) +=  dxinv * D_ip12;
                    //
                    rhs_viscosity[i] = -(
                           visc_ip12 * dxinv * (qtheta_ip1 - qtheta_0) - visc_ip12 * qtheta_ip12 / htheta_ip12 * dxinv * (htheta_ip1 - htheta_0)
                        - (visc_im12 * dxinv * (qtheta_0 - qtheta_im1) - visc_im12 * qtheta_im12 / htheta_im12 * dxinv * (htheta_0 - htheta_im1)
                          )
                        );
                    //rhs_viscosity[i] = std::abs(rhs_viscosity[i]);
                    rhs[q_eq] += rhs_viscosity[i];
                }
            }
            if (do_bed_shear_stress)
            {
                START_TIMER(Bed shear stress);
                status = bed_stress_matrix_rhs(A, rhs, hp, qp, hn, qn, cf_given, dx, theta);
                STOP_TIMER(Bed shear stress);
            }
            {
                //==============================================================
                // wwest boundary
                //==============================================================
                int i = 0;
                int ph = 2 * p_index(i, 0, nx);  // continuity equation
                int ph_e = 2 * p_index(i + 1, 0, nx);  // continuity equation (east)
                int ph_ee = 2 * p_index(i + 2, 0, nx);  // continuity equation (east-east)

                int c_eq = ph; 
                int q_eq = c_eq + 1;  // u-momentum equation 

                hn_0   = hn[i];       // = h^{n}_{i}
                hn_ip1 = hn[i + 1];       // = h^{n}_{i+1}
                hn_ip2 = hn[i + 2];       // = h^{n}_{i+2}
                hp_0   = hp[i];       // = h^{n+1,p}_{i}
                hp_ip1 = hp[i + 1];       // = h^{n+1,p}_{i+1}
                hp_ip2 = hp[i + 2];       // = h^{n+1,p}_{i+2}
                htheta_0   = htheta[i];
                htheta_ip1 = htheta[i + 1];
                htheta_ip2 = htheta[i + 2];

                qn_0   = qn[i];       // = q^{n}_{i}
                qn_ip1 = qn[i + 1];       // = q^{n}_{i+1}
                qn_ip2 = qn[i + 2];       // = q^{n}_{i+2}
                qp_0   = qp[i];       // = q^{n+1,p}_{i}
                qp_ip1 = qp[i + 1];       // = q^{n+1,p}_{i+1}
                qp_ip2 = qp[i + 2];       // = q^{n+1,p}_{i+1}
                qtheta_0   = qtheta[i];
                qtheta_ip1 = qtheta[i + 1];
                qtheta_ip2 = qtheta[i + 2];

                double htheta_b = w_ess[0] * htheta_0 + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2;
                double zb_b = w_ess[0] * zb[i] + w_ess[1] * zb[i + 1] + w_ess[2] * zb[i + 2];

                A.coeffRef(c_eq, ph) = 0.0;
                A.coeffRef(c_eq, ph_e) = 0.0;
                A.coeffRef(c_eq, ph_ee) = 0.0;
                //
                A.coeffRef(c_eq, ph + 1) = 0.0;
                A.coeffRef(c_eq, ph_e + 1) = 0.0;
                A.coeffRef(c_eq, ph_ee + 1) = 0.0;
                //
                rhs[c_eq] = 0.0;
                //
                double h_given = bc[BC_WEST] - zb_b;
                double h_infty = h_given;  // s_offset - zb_b;
                double c_wave = std::sqrt(g * htheta_b);
                if (bc_type[BC_WEST] == "mooiman")
                {
                    // q
                    A.coeffRef(c_eq, ph + 1) = w_ess[0];
                    A.coeffRef(c_eq, ph_e + 1) = w_ess[1];
                    A.coeffRef(c_eq, ph_ee + 1) = w_ess[2];
                    //  c_wave * h
                    A.coeffRef(c_eq, ph) = w_ess[0] * c_wave;
                    A.coeffRef(c_eq, ph_e) = w_ess[1] * c_wave;
                    A.coeffRef(c_eq, ph_ee) = w_ess[2] * c_wave;
                    //
                    rhs[c_eq] = -(qp_ip12  + c_wave * (hp_ip12 - h_infty));
                    if (bc_vars[BC_WEST] == "zeta" ) { rhs[c_eq] += 2. * c_wave * bc[BC_WEST]; }
                    if (bc_vars[BC_WEST] == "u" ) { rhs[c_eq] += 2. * h_infty * bc[BC_WEST]; }
                    if (bc_vars[BC_WEST] == "q" ) { rhs[c_eq] += 2. * bc[BC_WEST]; }
                    if (bc_vars[BC_WEST] == "w" ) { rhs[c_eq] += 2. * bc[BC_WEST]; }
                }
                else if (bc_type[BC_WEST] == "borsboom")
                {
                    //----------------------------------------------------------------------
                    // Essential boundary condition
                    //----------------------------------------------------------------------
                    double hn_b = w_ess[0] * hn_0 + w_ess[1] * hn_ip1 + w_ess[2] * hn_ip2;
                    double hp_b = w_ess[0] * hp_0 + w_ess[1] * hp_ip1 + w_ess[2] * hp_ip2;
                    double htheta_b = w_ess[0] * htheta_0 + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2;

                    double qn_b = w_ess[0] * qn_0 + w_ess[1] * qn_ip1 + w_ess[2] * qn_ip2;
                    double qp_b = w_ess[0] * qp_0 + w_ess[1] * qp_ip1 + w_ess[2] * qp_ip2;
                    double qtheta_b = w_ess[0] * qtheta_0 + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2;
                    //
                    // momentum + c_wave * continuity (ingoing signal)
                    // 
                    double con_fac = c_wave;
                    if (do_convection) { con_fac = c_wave - qp_b / hp_b; }
                    //
                    // continuity part (added and multiplied by c_wave)
                    //
                    A.coeffRef(c_eq, ph   ) = dtinv * con_fac * w_ess[0];
                    A.coeffRef(c_eq, ph_e ) = dtinv * con_fac * w_ess[1];
                    A.coeffRef(c_eq, ph_ee) = dtinv * con_fac * w_ess[2];
                    // 
                    // momentum part
                    // 
                    A.coeffRef(c_eq, ph    + 1) = dtinv * w_ess[0];
                    A.coeffRef(c_eq, ph_e  + 1) = dtinv * w_ess[1];
                    A.coeffRef(c_eq, ph_ee + 1) = dtinv * w_ess[2];
                    //
                    double dhdt = dtinv * (hp_b - hn_b);
                    double dqdt = dtinv * (qp_b - qn_b);
                    rhs[c_eq] = - ( dqdt + con_fac * dhdt );

                    double corr_term = 0.0;
                    if (bc_vars[BC_WEST] == "zeta")
                    {
                        A.coeffRef(c_eq, ph   ) += -dtinv * w_ess[0] - eps_bc_corr * theta * w_ess[0];
                        A.coeffRef(c_eq, ph_e ) += -dtinv * w_ess[1] - eps_bc_corr * theta * w_ess[1];
                        A.coeffRef(c_eq, ph_ee) += -dtinv * w_ess[2] - eps_bc_corr * theta * w_ess[2];
                        corr_term = ( dhdt + eps_bc_corr * ((bc[BC_WEST] - zb_b) - htheta_b) );
                        rhs[c_eq] += corr_term;
                    }
                    if (bc_vars[BC_WEST] == "q")
                    {
                        A.coeffRef(c_eq, ph    + 1) += dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0];
                        A.coeffRef(c_eq, ph_e  + 1) += dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1];
                        A.coeffRef(c_eq, ph_ee + 1) += dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2];
                        corr_term = -dqdt + eps_bc_corr * (bc[BC_WEST] - qtheta_b);
                        rhs[c_eq] += corr_term;
                    }
                }
                //----------------------------------------------------------------------
                // Natural boundary condition (west)
                //----------------------------------------------------------------------
                double hn_b = w_nat[0] * hn_0 + w_nat[1] * hn_ip1 + w_nat[2] * hn_ip2;
                double hp_b = w_nat[0] * hp_0 + w_nat[1] * hp_ip1 + w_nat[2] * hp_ip2;
                htheta_b = w_nat[0] * htheta_0 + w_nat[1] * htheta_ip1 + w_nat[2] * htheta_ip2;

                double qn_b = w_nat[0] * qn_0 + w_nat[1] * qn_ip1 + w_nat[2] * qn_ip2;
                double qp_b = w_nat[0] * qp_0 + w_nat[1] * qp_ip1 + w_nat[2] * qp_ip2;
                double qtheta_b = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_ip1 + w_nat[2] * qtheta_ip2;

                zb_b = w_nat[0] * zb[i] + w_nat[1] * zb[i + 1] + w_nat[2] * zb[i + 2];

                A.coeffRef(q_eq, ph   ) = 0.0;
                A.coeffRef(q_eq, ph_e ) = 0.0;
                A.coeffRef(q_eq, ph_ee) = 0.0;
                //
                A.coeffRef(q_eq, ph    + 1) = 0.0;
                A.coeffRef(q_eq, ph_e  + 1) = 0.0;
                A.coeffRef(q_eq, ph_ee + 1) = 0.0;
                //
                rhs[q_eq] = 0.0;
                //
                // momentum - c_wave * continuity
                // 
                double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                    + dtinv * (hp_ip1 - hn_ip1) * w_nat[1]
                    + dtinv * (hp_ip2 - hn_ip2) * w_nat[2];
                double dqdx = dxinv * (qtheta_ip1 - qtheta_0);
                double dqdt = dtinv * (qp_0 - qn_0) * w_nat[0]
                    + dtinv * (qp_ip1 - qn_ip1) * w_nat[1]
                    + dtinv * (qp_ip2 - qn_ip2) * w_nat[2];
                double dzetadx = dxinv * (htheta_ip1 + zb[i + 1] - htheta_0 - zb[i]);
                // 
                // momentum part dq/dt + gh d(zeta)/dx
                // 
                A.coeffRef(q_eq, ph   ) += w_nat[0] * theta * g * dzetadx - dxinv * theta * g * htheta_b;
                A.coeffRef(q_eq, ph_e ) += w_nat[1] * theta * g * dzetadx + dxinv * theta * g * htheta_b;
                A.coeffRef(q_eq, ph_ee) += w_nat[2] * theta * g * dzetadx;

                A.coeffRef(q_eq, ph + 1   ) += dtinv * w_nat[0];
                A.coeffRef(q_eq, ph_e + 1 ) += dtinv * w_nat[1];
                A.coeffRef(q_eq, ph_ee + 1) += dtinv * w_nat[2];
                rhs[q_eq] += -( dqdt + g * htheta_b * dzetadx );
                if (do_convection) // 
                {
                    double aa = - dxinv * 2. * qtheta_b / (htheta_b * htheta_b) * (qtheta_ip1 - qtheta_0)
                        + dxinv * 2. * (qtheta_b * qtheta_b) / (htheta_b * htheta_b * htheta_b) * (htheta_ip1 - htheta_0);
                    double bb = dxinv * 2. / htheta_b * (qtheta_ip1 - qtheta_0) - dxinv * 2. * qtheta_b / (htheta_b * htheta_b) * (htheta_ip1 - htheta_0);
                    double cc = -(qtheta_b * qtheta_b) / (htheta_b * htheta_b);
                    double dd = 2. * qtheta_b / htheta_b;

                    A.coeffRef(q_eq, ph   ) += theta * aa * w_nat[0] + dxinv * theta * cc;
                    A.coeffRef(q_eq, ph_e ) += theta * aa * w_nat[1] - dxinv * theta * cc;
                    A.coeffRef(q_eq, ph_ee) += theta * aa * w_nat[2];
                    A.coeffRef(q_eq, ph    + 1) += theta * bb * w_nat[0] - dxinv * theta * dd;
                    A.coeffRef(q_eq, ph_e  + 1) += theta * bb * w_nat[1] + dxinv * theta * dd;
                    A.coeffRef(q_eq, ph_ee + 1) += theta * bb * w_nat[2];
                    rhs[q_eq] += -(
                          dxinv * dd * (qtheta_ip1 - qtheta_0) +
                        + dxinv * cc * (htheta_ip1 - htheta_0)
                        );
                }
                if (do_bed_shear_stress) // 
                {
                    double cf_i = cf[i];
                    double cf_ip1 = cf[i + 1];

                    double cf_ip12 = scv(cf_i, cf_ip1);
                    double htheta_b = scv(htheta_0, htheta_ip1);
                    double qtheta_b = scv(qtheta_0, qtheta_ip1);
                    double abs_qtheta_b = scv(Fabs(qtheta_0, eps_abs), Fabs(qtheta_ip1, eps_abs));

                    A.coeffRef(q_eq, ph   ) += -0.5 * theta * cf_ip12 * 2. * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph_e ) += -0.5 * theta * cf_ip12 * 2. * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph_ee) += 0.0;

                    double J1_ip12 = cf_ip12 * qtheta_b / (htheta_b * htheta_b) * dxinv * (Fabs(qtheta_ip1, eps_abs) - Fabs(qtheta_0, eps_abs));
                    double J2_ip12 = cf_ip12 * abs_qtheta_b / (htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph    + 1) += 0.5 * theta * (J1_ip12 + J2_ip12);
                    A.coeffRef(q_eq, ph_e  + 1) += 0.5 * theta * (J1_ip12 + J2_ip12);
                    A.coeffRef(q_eq, ph_ee + 1) += 0.0;
                    rhs[q_eq] += -(
                        cf_ip12 * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b)
                        );
                }
                if (do_viscosity) // 
                {
                    double visc_0   = visc[i];
                    double visc_ip1 = visc[i + 1];
                    double visc_ip2 = visc[i + 2];
                    double visc_b = w_nat[0] * visc_0 + w_nat[1] * visc_ip1 + w_nat[2] * visc_ip2;
                    double htheta_b = w_nat[0] * htheta_0 + w_nat[1] * htheta_ip1 + w_nat[2] * htheta_ip2;
                    double qtheta_b = w_nat[0] * htheta_0 + w_nat[1] * qtheta_ip1 + w_nat[2] * qtheta_ip2;
                    double dviscdx = dxinv * (visc_ip1 - visc_0);
                    double dhdx = dxinv * (htheta_ip1 - htheta_0);
                    double dqdx = dxinv * (qtheta_ip1 - qtheta_0);
                    double d2qdx2 = dxinv * dxinv * (qtheta_0 - 2. * qtheta_ip1 + qtheta_ip2);
                    double d2hdx2 = dxinv * dxinv * (htheta_0 - 2. * htheta_ip1 + htheta_ip2);
                    double viscos = (
                        dviscdx * (dqdx - qtheta_b / htheta_b * dhdx)
                        + visc_b * d2qdx2 
                        - visc_b * 1.0 / htheta_b * dhdx * dqdx
                        + visc_b * qtheta_b / (htheta_b * htheta_b) * dhdx * dhdx 
                        - visc_b * qtheta_b / htheta_b * d2hdx2
                        );
                    rhs[q_eq] += std::abs(viscos);
                }
                //
                // continuity part (added and multiplied by -c_wave)
                //
                double con_fac = c_wave;
                if (do_convection) { con_fac = c_wave + qp_ip12 / hp_ip12; }
                A.coeffRef(q_eq, ph   ) += -con_fac * dtinv * w_nat[0];
                A.coeffRef(q_eq, ph_e ) += -con_fac * dtinv * w_nat[1];
                A.coeffRef(q_eq, ph_ee) += -con_fac * dtinv * w_nat[2];

                A.coeffRef(q_eq, ph    + 1) += -con_fac * dxinv * -theta;
                A.coeffRef(q_eq, ph_e  + 1) += -con_fac * dxinv *  theta;
                A.coeffRef(q_eq, ph_ee + 1) += 0.0;
                rhs[q_eq] += con_fac * (dhdt + dqdx);
            }
            {
                //==============================================================
                // eeast boundary
                //==============================================================
                int i = nx - 1;
                int ph = 2 * p_index(i, 0, nx);  // continuity equation
                int ph_w = 2 * p_index(i - 1, 0, nx);  // continuity equation
                int ph_ww = 2 * p_index(i - 2, 0, nx);  // continuity equation

                int c_eq = ph;
                int q_eq = c_eq + 1;  // u-momentum equation

                hn_0   = hn[i];       // = h^{n}_{i}
                hn_im1 = hn[i - 1];       // = h^{n}_{i-1}
                hn_im2 = hn[i - 2];       // = h^{n}_{i-2}
                hp_0   = hp[i];       // = h^{n+1,p}_{i}
                hp_im1 = hp[i - 1];       // = h^{n+1,p}_{i-1}
                hp_im2 = hp[i - 2];       // = h^{n+1,p}_{i-2}
                htheta_0   = htheta[i];
                htheta_im1 = htheta[i - 1];
                htheta_im2 = htheta[i - 2];

                qn_0   = qn[i];       // = q^{n}_{i}
                qn_im1 = qn[i - 1];       // = q^{n}_{i-1}
                qn_im2 = qn[i - 2];       // = q^{n}_{i-1}
                qp_0   = qp[i];       // = q^{n+1,p}_{i}
                qp_im1 = qp[i - 1];       // = q^{n+1,p}_{i-1}
                qp_im2 = qp[i - 2];       // = q^{n+1,p}_{i-2}
                qtheta_0   = qtheta[i];
                qtheta_im1 = qtheta[i - 1];
                qtheta_im2 = qtheta[i - 2];

                double htheta_b = w_ess[0] * htheta_0 + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2;
                double zb_b = w_ess[0] * zb[i] + w_ess[1] * zb[i - 1] + w_ess[2] * zb[i - 2];

                A.coeffRef(c_eq, ph) = 0.0;
                A.coeffRef(c_eq, ph_w) = 0.0;
                A.coeffRef(c_eq, ph_ww) = 0.0;
                //
                A.coeffRef(c_eq, ph + 1) = 0.0;
                A.coeffRef(c_eq, ph_w + 1) = 0.0;
                A.coeffRef(c_eq, ph_ww + 1) = 0.0;
                //
                rhs[c_eq] = 0.0;

                double h_given = bc[BC_EAST] - zb_b;
                double h_infty = h_given;  // s_offset - zb_b;
                double c_wave = std::sqrt(g * htheta_b);

                if (bc_type[BC_EAST] == "mooiman")
                {
                    // q
                    A.coeffRef(c_eq, ph + 1) = w_ess[0];
                    A.coeffRef(c_eq, ph_w + 1) = w_ess[1];
                    A.coeffRef(c_eq, ph_ww + 1) = w_ess[2];
                    // c_wave * h
                    A.coeffRef(c_eq, ph) = -w_ess[0] * c_wave;
                    A.coeffRef(c_eq, ph_w) = -w_ess[1] * c_wave;
                    A.coeffRef(c_eq, ph_ww) = -w_ess[2] * c_wave;
                    //
                    rhs[c_eq] = -(qp_im12 - c_wave * (hp_im12 - h_infty));
                    if (bc_vars[BC_EAST] == "zeta" ) { rhs[c_eq] -= 2. * c_wave * bc[BC_EAST]; }
                    if (bc_vars[BC_EAST] == "u" ) { rhs[c_eq] -= 2. * h_infty * bc[BC_EAST]; }
                    if (bc_vars[BC_EAST] == "q" ) { rhs[c_eq] -= 2. * bc[BC_EAST]; }
                    if (bc_vars[BC_EAST] == "w" ) { rhs[c_eq] -= 2. * bc[BC_EAST]; }
                }
                else if (bc_type[BC_EAST] == "borsboom")
                {
                    //----------------------------------------------------------------------
                    // Essential boundary condition
                    //----------------------------------------------------------------------
                    double hn_b = w_ess[0] * hn_0 + w_ess[1] * hn_im1 + w_ess[2] * hn_im2;
                    double hp_b = w_ess[0] * hp_0 + w_ess[1] * hp_im1 + w_ess[2] * hp_im2;
                    double htheta_b = w_ess[0] * htheta_0 + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2;

                    double qn_b = w_ess[0] * qn_0 + w_ess[1] * qn_im1 + w_ess[2] * qn_im2;
                    double qp_b = w_ess[0] * qp_0 + w_ess[1] * qp_im1 + w_ess[2] * qp_im2;
                    double qtheta_b = w_ess[0] * qtheta_0 + w_ess[1] * qtheta_im1 + w_ess[2] * qtheta_im2;

                    //
                    // momentum - c_wave * continuity (ingoing signal)
                    // 
                    double con_fac = c_wave;
                    if (do_convection) { con_fac = c_wave + qp_b / hp_b; }
                    //
                    // continuity part (added and multiplied by -c_wave)
                    //
                    A.coeffRef(c_eq, ph   ) = dtinv * -con_fac * w_ess[0];
                    A.coeffRef(c_eq, ph_w ) = dtinv * -con_fac * w_ess[1];
                    A.coeffRef(c_eq, ph_ww) = dtinv * -con_fac * w_ess[2];
                    // 
                    // momentum part
                    // 
                    A.coeffRef(c_eq, ph    + 1) = dtinv * w_ess[0];
                    A.coeffRef(c_eq, ph_w  + 1) = dtinv * w_ess[1];
                    A.coeffRef(c_eq, ph_ww + 1) = dtinv * w_ess[2];
                    //
                    double dhdt = dtinv * (hp_b - hn_b);
                    double dqdt = dtinv * (qp_b - qn_b);
                    rhs[c_eq] = - (dqdt - con_fac * dhdt);

                    double corr_term = 0.0;
                    if (bc_vars[BC_EAST] == "zeta")
                    {
                        if (stationary) { sign = -1.0; }
                        A.coeffRef(c_eq, ph   ) += dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0];
                        A.coeffRef(c_eq, ph_w ) += dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1];
                        A.coeffRef(c_eq, ph_ww) += dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2];
                        corr_term = - (dhdt + sign * eps_bc_corr * ((bc[BC_EAST] - zb_b) - htheta_b));
                        rhs[c_eq] += corr_term;
                        sign = 1.0;
                    }
                    if (bc_vars[BC_EAST] == "q")
                    {
                        if (stationary) { sign = -1.0; }
                        A.coeffRef(c_eq, ph    + 1) += dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0];
                        A.coeffRef(c_eq, ph_w  + 1) += dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1];
                        A.coeffRef(c_eq, ph_ww + 1) += dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2];
                        corr_term = (-dqdt + sign * eps_bc_corr * (bc[BC_EAST] - qtheta_b));
                        rhs[c_eq] += corr_term;
                        sign = 1.0;
                    }
                }
                //----------------------------------------------------------------------
                // Natural boundary condition (east)
                //----------------------------------------------------------------------
                double hn_b = w_nat[0] * hn_0 + w_nat[1] * hn_im1 + w_nat[2] * hn_im2;
                double hp_b = w_nat[0] * hp_0 + w_nat[1] * hp_im1 + w_nat[2] * hp_im2;
                htheta_b = w_nat[0] * htheta_0 + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;

                double qn_b = w_nat[0] * qn_0 + w_nat[1] * qn_im1 + w_nat[2] * qn_im2;
                double qp_b = w_nat[0] * qp_0 + w_nat[1] * qp_im1 + w_nat[2] * qp_im2;
                double qtheta_b = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;

                zb_b = w_nat[0] * zb[i] + w_nat[1] * zb[i - 1] + w_nat[2] * zb[i - 2];

                //            
                A.coeffRef(q_eq, ph   ) = 0.0;
                A.coeffRef(q_eq, ph_w ) = 0.0;
                A.coeffRef(q_eq, ph_ww) = 0.0;
                //
                A.coeffRef(q_eq, ph    + 1) = 0.0;
                A.coeffRef(q_eq, ph_w  + 1) = 0.0;
                A.coeffRef(q_eq, ph_ww + 1) = 0.0;
                //
                rhs[q_eq] = 0.0;
                //
                // momentum + c_wave * continuity
                // 
                double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                    + dtinv * (hp_im1 - hn_im1) * w_nat[1]
                    + dtinv * (hp_im2 - hn_im2) * w_nat[2];
                double dqdx = dxinv * (qtheta_0 - qtheta_im1);
                double dqdt = dtinv * (qp_0 - qn_0) * w_nat[0]
                    + dtinv * (qp_im1 - qn_im1) * w_nat[1]
                    + dtinv * (qp_im2 - qn_im2) * w_nat[2];
                double dzetadx = dxinv * (htheta_0 + zb[i] - htheta_im1 - zb[i - 1]);
                // 
                // momentum part dq/dt + gh d(zeta)/dx
                // 
                A.coeffRef(q_eq, ph   ) += w_nat[0] * theta * g * dzetadx + dxinv * theta * g * htheta_b;
                A.coeffRef(q_eq, ph_w ) += w_nat[1] * theta * g * dzetadx - dxinv * theta * g * htheta_b;
                A.coeffRef(q_eq, ph_ww) += w_nat[2] * theta * g * dzetadx;
                A.coeffRef(q_eq, ph    + 1) += dtinv * w_nat[0];
                A.coeffRef(q_eq, ph_w  + 1) += dtinv * w_nat[1];
                A.coeffRef(q_eq, ph_ww + 1) += dtinv * w_nat[2];
                rhs[q_eq] += - ( dqdt + g * htheta_b * dzetadx );
                if (do_convection)
                {
                    double aa = - dxinv * 2. * qtheta_b / (htheta_b * htheta_b) * (qtheta_0 - qtheta_im1)
                        + dxinv * 2. * (qtheta_b * qtheta_b) / (htheta_b * htheta_b * htheta_b) * (htheta_0 - htheta_im1);
                    double bb = dxinv * 2. / htheta_b * (qtheta_0 - qtheta_im1) - dxinv * 2. * qtheta_b / (htheta_b * htheta_b) * (htheta_0 - htheta_im1);
                    double cc = -(qtheta_b * qtheta_b) / (htheta_b * htheta_b);
                    double dd = 2. * qtheta_b / htheta_b;

                    A.coeffRef(q_eq, ph   ) += theta * aa * w_nat[0] - dxinv * theta * cc;
                    A.coeffRef(q_eq, ph_w ) += theta * aa * w_nat[1] + dxinv * theta * cc;
                    A.coeffRef(q_eq, ph_ww) += theta * aa * w_nat[2];
                    A.coeffRef(q_eq, ph    + 1) += theta * bb * w_nat[0] + dxinv * theta * dd;
                    A.coeffRef(q_eq, ph_w  + 1) += theta * bb * w_nat[1] - dxinv * theta * dd;
                    A.coeffRef(q_eq, ph_ww + 1) += theta * bb * w_nat[2];
                    rhs[q_eq] += -(
                        dxinv * dd * (qtheta_0 - qtheta_im1)
                        + dxinv * cc * (htheta_0 - htheta_im1)
                        );
                }
                if (do_bed_shear_stress) // 
                {
                    double cf_im1 = cf[i - 1];
                    double cf_i = cf[i];
                    double htheta_b = w_nat[0] * htheta_0 + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;
                    double qtheta_b = w_nat[0] * htheta_0 + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;
                    double abs_qtheta_b = scv(Fabs(qtheta_0, eps_abs), Fabs(qtheta_im1, eps_abs));
                    double cf_b = scv(cf_i, cf_im1);
                    //
                    // bed_shear_stress
                    //
                    A.coeffRef(q_eq, ph   ) += -0.5 * theta * cf_b * 2. * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph_w ) += -0.5 * theta * cf_b * 2. * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph_ww) += 0.0;
                    //
                    double J1_b = cf_b * qtheta_b / (htheta_b * htheta_b) * dxinv * (Fabs(qtheta_0, eps_abs) - Fabs(qtheta_im1, eps_abs));
                    double J2_b = cf_b * abs_qtheta_b / (htheta_b * htheta_b);
                    A.coeffRef(q_eq, ph    + 1) += 0.5 * theta * (J1_b + J2_b);
                    A.coeffRef(q_eq, ph_w  + 1) += 0.5 * theta * (J1_b + J2_b);
                    A.coeffRef(q_eq, ph_ww + 1) += 0.0;
                    //
                    double rhs_bed_stress = -(
                        cf_b * qtheta_b * abs_qtheta_b / (htheta_b * htheta_b)
                        );
                    rhs[q_eq] += rhs_bed_stress;
                }
                if (do_viscosity) // 
                {
                    double visc_0   = visc[i];
                    double visc_im1 = visc[i - 1];
                    double visc_im2 = visc[i - 2];
                    double visc_b = w_nat[0] * visc_0 + w_nat[1] * visc_im1 + w_nat[2] * visc_im2;
                    double htheta_b = w_nat[0] * htheta_0 + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;
                    double qtheta_b = w_nat[0] * htheta_0 + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;
                    double dviscdx = dxinv * (visc_0 - visc_im1);
                    double dhdx = dxinv * (htheta_0 - htheta_im1);
                    double dqdx = dxinv * (qtheta_0 - qtheta_im1);
                    double d2qdx2 = dxinv * dxinv * (qtheta_0 - 2. * qtheta_im1 + qtheta_im2);
                    double d2hdx2 = dxinv * dxinv * (htheta_0 - 2. * htheta_im1 + htheta_im2);
                    double viscos = (
                        dviscdx * (dqdx - qtheta_b / htheta_b * dhdx)
                        + visc_b * d2qdx2 
                        - visc_b * 1.0 / htheta_b * dhdx * dqdx
                        + visc_b * qtheta_b / (htheta_b * htheta_b) * dhdx * dhdx 
                        - visc_b * qtheta_b / htheta_b * d2hdx2
                        );
                    //viscos = dxinv * dxinv * (qtheta_0/htheta_0 - 2. * qtheta_im1/htheta_im1 + qtheta_im2/htheta_im2);
                    //viscos = dxinv * dxinv * (qtheta_im1/htheta_im1 - 2. * qtheta_im1/htheta_im1 + qtheta_im2/htheta_im2);
                    rhs[q_eq] += std::abs(viscos);
                }
                //
                // continuity part (added and multiplied by +c_wave)
                //
                double con_fac = c_wave;
                if (do_convection) { con_fac = c_wave - qp_b / hp_b; }
                A.coeffRef(q_eq, ph   ) += con_fac * dtinv * w_nat[0];
                A.coeffRef(q_eq, ph_w ) += con_fac * dtinv * w_nat[1];
                A.coeffRef(q_eq, ph_ww) += con_fac * dtinv * w_nat[2];
                A.coeffRef(q_eq, ph    + 1) += con_fac * dxinv *  theta;
                A.coeffRef(q_eq, ph_w  + 1) += con_fac * dxinv * -theta;
                A.coeffRef(q_eq, ph_ww + 1) += 0.0;
                rhs[q_eq] += -con_fac * (dhdt + dqdx);
            }
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
            if (logging == "matrix" && (nst == 1 || nst == total_time_steps-1) && iter == 0)
            {
                std::string header_text = "=== Matrix ============================================";
                print_matrix(A, 2, nx, ny, header_text, log_file);
                header_text = "=== RHS ===============================================";
                print_vector(rhs, 2, nx, ny, header_text, log_file);
            }

            solver.compute(A);
            solver.setTolerance(eps_bicgstab);
            solution = solver.solve(rhs);
            //solution = solver.solveWithGuess(rhs, solution);
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
            // 
            // The new solution is the previous iterant plus the delta
            //
            dh_max = 0.0;
            dq_max = 0.0;
            dh_maxi = 0;
            dq_maxi = 0;
            for (size_t i = 0; i < nx; ++i)
            {
                hp[i] += solution[2 * i]; // h, continuity-eq
                qp[i] += solution[2 * i + 1];  // q, momentum-eq
                delta_h[i] = solution[2 * i];  // h, continuity-eq
                delta_q[i] = solution[2 * i + 1];  // q, momentum-eq
                if (dh_max < std::abs(delta_h[i]))
                {
                    dh_max = std::abs(delta_h[i]);
                    dh_maxi = i;
                }
                if (dq_max < std::abs(delta_q[i]))
                {
                    dq_max = std::abs(delta_q[i]);
                    dq_maxi = i;
                }
            }

            if (logging == "matrix" && (nst == 1 || nst == total_time_steps-1) && iter == 0)
            {
                std::string header_text = "=== Solution ==========================================";
                print_vector(solution, 2, nx, ny, header_text, log_file);

                log_file << "=== hp, qp, zeta ======================================" << std::endl;
                for (size_t i = 0; i < nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << hp[i] << ", " << qp[i] << ", " << hp[i] + zb[i] << std::endl;
                }
                log_file << "=======================================================" << std::endl;
            }
            used_lin_solv_iter = std::max(used_lin_solv_iter, (int)solver.iterations());
            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                //for (size_t i = 0; i < nr_nodes; ++i)
                //{
                //    u[i] = qp[i] / hp[i];
                //}
                //(void)regular->first_derivative(psi, visc_reg, u, dx);
                STOP_TIMER(Regularization_iter_loop);
                //for (size_t i = 0; i < nr_nodes; ++i)
                //{
                    //visc[i] = visc_reg[i] * std::abs(psi[i]);
                    //pe[i] = qp[i] / hp[i] * dx / visc[i];
                //}
            }
            if (dh_max < eps_newton && dq_max < eps_newton)
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
                << "    Delta h^{n + 1,p + 1}: " << dh_max << " at index: " << dh_maxi << std::endl
                << "    Delta q^{n + 1,p + 1}: " << dq_max << " at index: " << dq_maxi << std::endl;
        }
        else
        {
            if (std::fmod(time, input_data.output.dt_screen) == 0)
            {
                std::cout << std::fixed << std::setprecision(2) << tstart + time << ";   " << tstart + tstop << std::endl;
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(4) << std::scientific << time
                    << std::setprecision(8) << std::scientific << std::endl
                    << "    Newton iterations  : " << used_newton_iter << std::endl
                    << "    Delta h^{n + 1,p + 1}: " << dh_max << " at index: " << dh_maxi << std::endl
                    << "    Delta q^{n + 1,p + 1}: " << dq_max << " at index: " << dq_maxi << std::endl;
            }
        }
        if (used_newton_iter == iter_max)
        {
            if (dh_max > eps_newton || dq_max > eps_newton)
            {
                log_file << "    ----    maximum number of iterations reached, probably not converged, at time: " <<  time << " [sec]" << std::endl;
            }
        }
        for (size_t i = 0; i < nx; ++i)
        {
            hn[i] = hp[i];  // h, continuity-eq
            qn[i] = qp[i];  // q, momentum-eq
            s[i] = hn[i] + zb[i];
            u[i] = qn[i]/hn[i];
            riemann_pos[i] = std::abs(qn[i] + std::sqrt(g * hn[i]) * hn[i]);
            riemann_neg[i] = std::abs(qn[i] - std::sqrt(g * hn[i]) * hn[i]);
        }

        if (regularization_time)
        {
            if (do_viscosity)
            {
                START_TIMER(Regularization_time_loop);
                (void)regularization->artificial_viscosity(psi, hp, qp, zb, c_psi, dx);
                for (size_t i = 0; i < nx; ++i)
                {
                    visc[i] = visc_reg[i] + std::abs(psi[i]);
                }
                STOP_TIMER(Regularization_time_loop);
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            froude[i] = qp[i] / hp[i] / std::sqrt(g * hp[i]);
        }
        if (do_viscosity)
        {
            for (int i = 0; i < nx; ++i)
            {
                pe[i] = qp[i] / hp[i] * dx / visc[i];
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
            map_file->put_time_variable(map_names[0], nst_map, hn);
            map_file->put_time_variable(map_names[1], nst_map, qn);
            map_file->put_time_variable(map_names[2], nst_map, s);
            map_file->put_time_variable(map_names[3], nst_map, u);
            map_file->put_time_variable(map_names[4], nst_map, zb);
            map_file->put_time_variable(map_names[5], nst_map, delta_h);
            map_file->put_time_variable(map_names[6], nst_map, delta_q);
            map_file->put_time_variable(map_names[10], nst_map, froude);
            if (do_viscosity)
            {
                map_file->put_time_variable(map_names[7], nst_map, visc_reg);
                map_file->put_time_variable(map_names[8], nst_map, visc);
                map_file->put_time_variable(map_names[9], nst_map, psi);
                map_file->put_time_variable(map_names[11], nst_map, pe);
                map_file->put_time_variable(map_names[12], nst_map, rhs_viscosity);
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

            status = set_his_values(input_data.obs_points, hn, his_values);
            his_file->put_variable(his_h_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, qn, his_values);
            his_file->put_variable(his_q_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, s, his_values);
            his_file->put_variable(his_s_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, u, his_values);
            his_file->put_variable(his_u_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, zb, his_values);
            his_file->put_variable(his_zb_name, nst_his, his_values);

            status = set_his_values(input_data.obs_points, froude, his_values);
            his_file->put_variable(his_froude_name, nst_his, his_values);

            if (do_viscosity)
            {
                status = set_his_values(input_data.obs_points, pe, his_values);
                his_file->put_variable(his_peclet_name, nst_his, his_values);
                status = set_his_values(input_data.obs_points, visc, his_values);
                his_file->put_variable(his_visc_name, nst_his, his_values);
                status = set_his_values(input_data.obs_points, psi, his_values);
                his_file->put_variable(his_psi_name, nst_his, his_values);
                status = set_his_values(input_data.obs_points, rhs_viscosity, his_values);
                his_file->put_variable(his_visc_rhs_name, nst_his, his_values);
            }

            //his_values = { riemann_pos[p_a], riemann_pos[p_b], riemann_pos[p_c], riemann_pos[p_d], riemann_pos[p_e] };
            //his_file->put_variable(his_riemann_pos_name, nst_his, his_values);
            //his_values = { riemann_neg[p_a], riemann_neg[p_a], riemann_neg[p_c], riemann_neg[p_d], riemann_neg[p_e] };
            //his_file->put_variable(his_riemann_neg_name, nst_his, his_values);

            his_values.clear();
            his_values = { double(used_newton_iter) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            his_values.clear();
            his_values = { double(solver.iterations()) };
            his_file->put_variable(his_LinSolver_iter_name, nst_his, his_values);
            his_values.clear();
            his_values = { solver.error() };
            his_file->put_variable(his_LinSolver_iter_error_name, nst_his, his_values);

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
size_t p_index(size_t i, size_t j, size_t nx)
{
    return j * nx + i;
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

int initialize_scalar(double alpha, std::vector<double>& value_in, std::vector<double>& value_out)
{
    size_t nx = value_in.size();
    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector
    Eigen::VectorXd rhs(nx);                // RHS vector

    size_t i = 0;
    A.coeffRef(i, i) = -alpha;
    A.coeffRef(i, i + 1) = alpha;
    rhs[i] = value_in[i];
    for (size_t i = 1; i < nx-1; ++i)
    {
        A.coeffRef(i, i - 1) = alpha;
        A.coeffRef(i, i) = 1.0 - 2.0 *alpha;
        A.coeffRef(i, i + 1) = alpha;
        rhs[i] = value_in[i];
    }
    i = nx - 1;
    A.coeffRef(i, i) = alpha;
    A.coeffRef(i, i - 1) = -alpha;
    rhs[i] = value_in[i];

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solution = solver.solveWithGuess(rhs, solution);

    for (size_t i = 0; i < nx; ++i)
    {
        value_out[i] = solution[i]; // h, continuity-eq
    }
    return 0;
}

int read_bed_level(std::string filename, std::vector<double> & value)
{
    int nrow, ncol;
    double dummy;
    std::string record;

    size_t nx = (size_t) value.size();
    std::ifstream input_file;

    input_file.open(filename);
    if (!input_file.is_open())
    {
        std::cout << "Cannot open input stream: " << filename << std::endl;
        exit(0);
    }
    for (std::string record; std::getline(input_file, record); )
    {
        if (record[0] != '*') 
        {
            break;
        }
    }
    input_file >> nrow >> ncol;
    if (nrow != (int)nx)
    {
        std::cout << "Dimension of bed_filename does not match;  " << nrow << " != " << nx << std::endl;
        return 1;
    }
    for (size_t k = 0; k < nx; ++k)
    {
        input_file >> dummy >> value[k];
    }
    input_file.close();
    return 0;
}

double scv(double h1, double h2)
{
    // value at subcontrol volume
    return 0.25 * (3. * h1 + h2);
}
double Fabs(double q, double eps_abs)
{
    // absolute value as continue function
    return pow( pow(q, 4.) + pow(eps_abs, 4.), 0.25);
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
