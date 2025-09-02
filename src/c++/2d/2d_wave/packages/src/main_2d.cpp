//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
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

#define _USE_MATH_DEFINES
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <string>
#include <thread>

#include <vector>

#include <toml.h>
// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "bed_level.h"
#include "bed_shear_stress.h"
#include "boundary_condition.h"
#include "build_matrix_pattern.h"
#include "cfts.h"
#include "compile_date_and_time.h"
#include "convection.h"
#include "grid.h"
#include "initial_conditions.h"
#include "perf_timer.h"
#include "matrix_assembly_boundaries.h"
#include "matrix_assembly_corners.h"
#include "matrix_assembly_interior.h"
#include "regularization.h"
#include "ugrid2d.h"

void GetArguments(long argc, char** argv, std::filesystem::path & file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

int idx(int i, int j, int ny);
double c_scv(double, double, double, double);
double scvf_n(double c0, double c1, double c2, double c3);
double dcdx_scv(double, double, double, double);
double dcdy_scv(double, double, double, double);
double dcdx_scvf_n(double, double, double, double);
double dcdx_scvf_t(double, double, double, double);
double dcdy_scvf_n(double, double, double, double);
double dcdy_scvf_n(double, double, double, double);
std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name);
std::string string_format_with_zeros(double value, int width);


// Solve the linear wave equation
// Continuity equation: d(h)/dt + d(q)/dx = 0
// Momentum equation: d(q)/dt + gh d(zeta)/dx = 0
// Momentum equation: d(r)/dt + gh d(zeta)/dy = 0

std::string compileDateTime()
{
    std::string str1(compileYear());
    std::string str2(compileMonth());
    std::string str3(compileDay());
    if (str3.size() == 1)
    {
        str3.resize(2);
        str3[1] = str3[0];
        str3[0] = '0';
    }
    return str1 + "-" + str2 + "-" + str3 + " " + __TIME__;
}
int main(int argc, char *argv[])
{
    bool stationary = false;
    double sign = 1.0;
    std::filesystem::path toml_file_name("---not-defined---");
    int status = -1;

    int BC_NORTH = 0;
    int BC_EAST = 1;
    int BC_SOUTH = 2;
    int BC_WEST = 3;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path start_dir;
    std::filesystem::path output_dir;
    std::filesystem::path input_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();
    start_dir = std::filesystem::current_path();

    toml::table tbl;
    toml::table tbl_chp;  // table for a chapter
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
        tbl = toml::parse_file(toml_file_name.c_str());
        // std::cout << tbl << "\n";
        output_dir = input_dir;
        output_dir += "/output/";
        std::filesystem::create_directory(output_dir);
    }
    else
    {
        std::cout << "======================================================" << std::endl;
        std::cout << "Executable compiled: " << compileDateTime() << std::endl;
        std::cout << std::endl;
        std::cout << "usage: 2d_wave.exe --toml <input_file>" << std::endl;
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
    std::cout << "Executable compiled : " << compileDateTime() << std::endl;
    std::cout << "Start time          : " << start_date_time << std::endl;
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

    std::string logging = tbl["Logging"].value_or("None");

    ss << "2d_wave";
    out_file = output_dir.string() + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Linear wave equation, BiCGstab");

    std::ofstream log_file;
    log_file.open(log_filename);

    std::cout << "=== Input file =======================================" << std::endl;
    std::cout << std::filesystem::absolute(toml_file_name) << std::endl;
    std::cout << "======================================================" << std::endl;
    log_file << std::endl;
    log_file << "======================================================" << std::endl;
    log_file << "Executable compiled: " << compileDateTime() << std::endl;
    log_file << "Start time         : " << start_date_time << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;
    log_file << "=== Copy of the input file ============================" << std::endl;
    log_file << tbl << "\n";  // Write input TOML file to log_file
    log_file << "=======================================================" << std::endl;
    log_file << std::endl;

    tbl_chp = *tbl["Numerics"].as_table();
    double dt = tbl_chp["dt"].value_or(double(0.0));  // dt == 0 => default stationary
    double dtpseu = tbl_chp["dtpseu"].value_or(double(0.0));  // dt == 0 => default stationary
    if (dt == 0.0) { stationary = true;  }
    // Time
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    double tstop = tbl_chp["tstop"].value_or(double(60.));
    
    // Initial
    std::vector<std::string> ini_vars;  // get the element as an array
    tbl_chp = *tbl["Initial"].as_table();
    status = get_toml_array(tbl_chp, "ini_vars", ini_vars);
    double gauss_amp  = tbl_chp["gauss_amp"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary
    double gauss_mu   = tbl_chp["gauss_mu"].value_or(double(-INFINITY));
    double gauss_mu_x = tbl_chp["gauss_mu_x"].value_or(double(0.0));
    double gauss_mu_y = tbl_chp["gauss_mu_y"].value_or(double(0.0));
    double gauss_sigma = tbl_chp["gauss_sigma"].value_or(double(-INFINITY));
    double gauss_sigma_x = tbl_chp["gauss_sigma_x"].value_or(double(1.0));
    double gauss_sigma_y = tbl_chp["gauss_sigma_y"].value_or(double(1.0));
    if (gauss_mu != -INFINITY)
    {
        gauss_mu_x = gauss_mu; 
        gauss_mu_y = 0.0;  // only shift on x-axis
    }
    if (gauss_sigma != -INFINITY)
    {
        gauss_sigma_x = gauss_sigma; 
        gauss_sigma_y = gauss_sigma; 
    }

    // Domain
    tbl_chp = *tbl["Domain"].as_table();
    std::string grid_filename = tbl_chp["mesh_file"].value_or("--none--");
    std::filesystem::path full_grid_filename = input_dir;
    full_grid_filename += grid_filename;
    std::string bed_level_filename = tbl_chp["bed_level_file"].value_or("--none--");
    std::filesystem::path full_bed_level_filename = input_dir;
    full_bed_level_filename += bed_level_filename;

    //Physics
    tbl_chp = *tbl["Physics"].as_table();
    double g = tbl_chp["g"].value_or(double(9.81));  // Gravitational acceleration
    bool do_continuity = tbl_chp["do_continuity"].value_or(bool(true));  // default, continuity
    bool do_q_equation = tbl_chp["do_q_equation"].value_or(bool(true));  // default, q_equation
    bool do_r_equation = tbl_chp["do_r_equation"].value_or(bool(true));  // default, r_equation
    bool do_convection = tbl_chp["do_convection"].value_or(bool(false));  // default, no convection
    
    bool do_viscosity = tbl_chp["do_viscosity"].value_or(bool(false));  // default, no viscosity
    double visc_const = tbl_chp["viscosity"].value_or(double(0.0001));  // default 1e-4

    bool do_bed_shear_stress = tbl_chp["do_bed_shear_stress"].value_or(bool(false));  // default, no bed shear stress
    double chezy_coefficient = tbl_chp["chezy_coefficient"].value_or(double(50.0));
    if (do_q_equation && do_r_equation)
    {
        model_title = "Linear wave equation, BiCGstab";
        if (do_bed_shear_stress)
        {
            model_title = "Linear wave equation + bed shear stress, BiCGstab";
        }
        if (do_convection && do_bed_shear_stress)
        {
            model_title = "Linear wave equation + convection + bed shear stress, BiCGstab";
        }
        if (do_convection)
        {
            model_title = "Linear wave equation + convection, BiCGstab";
        }
    }
    else if (do_q_equation)
    {
        model_title = "Linear wave equation, only q-equation, BiCGstab";
    }
    else if (do_r_equation)
    {
        model_title = "Linear wave equation, only r-equation, BiCGstab";
    }
    else
    {
        model_title = "No q- and no r-equation (=> no waves computed), BiCGstab";
    }

    // Boundary
    tbl_chp = *tbl["Boundary"].as_table();
    double eps_bc_corr = tbl_chp["eps_bc_corr"].value_or(double(0.0001));  // default 1e-4
    std::vector<std::string> bc_type;
    status = get_toml_array(tbl_chp, "bc_type", bc_type);
    double treg = tbl_chp["treg"].value_or(double(150.0));
    if (stationary) { treg = 0.0; }
    std::vector<std::string> bc_vars;  // get the element as an array
    status = get_toml_array(tbl_chp, "bc_vars", bc_vars);
    std::vector<double> bc_vals;
    status = get_toml_array(tbl_chp, "bc_vals", bc_vals);
    std::vector<bool> bc_absorbing;
    status = get_toml_array(tbl_chp, "bc_absorbing", bc_absorbing);

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
//    double dt = tbl["Numerics"]["dt"].value_or(double(1.0));  // default stationary
    if (dt == 0.0) { stationary = true; }
    double theta = tbl_chp["theta"].value_or(double(0.501));
    if (stationary) { theta = 1.0; }
    double c_psi = tbl_chp["c_psi"].value_or(double(4.));
    int iter_max = tbl_chp["iter_max"].value_or(int(50));
    double eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));
    double eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));
    double eps_abs = tbl_chp["eps_absolute"].value_or(double(1.0e-2));  // epsilon needed to approximate the abs-function by a continues function
    bool regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    bool regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    bool regularization_time = tbl_chp["regularization_time"].value_or(bool(false));

    // Output
    tbl_chp = *tbl["Output"].as_table();
    double dt_his = tbl_chp["dt_his"].value_or(double(1.0));  // write interval to his-file
    double dt_map = tbl_chp["dt_map"].value_or(double(0.0));  // write interval to his-file

                          
    struct _mesh2d * mesh2d;
    SGRID* sgrid = new SGRID();
    status = sgrid->open(full_grid_filename.string());
    if (status != 0) 
    {
        log_file << "Error: Failed to open grid file: " << full_grid_filename.string() << std::endl;
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    status = sgrid->read();  // reading mesh
    mesh2d = sgrid->get_mesh_2d();

    int nx = mesh2d->node[0]->dims[0]; // including virtual points
    int ny = mesh2d->node[0]->dims[1]; // including virtual points

    std::vector<double>& x = mesh2d->node[0]->x; // reference to original x-coordinate
    std::vector<double>& y = mesh2d->node[0]->y; // reference to original y-coordinate

    BED_LEVEL * bed = new BED_LEVEL();
    status = bed->open(full_bed_level_filename.string());
    if (status != 0)
    {
        std::cout << "Failed to open file: " << full_bed_level_filename.string() << std::endl;
        log_file << "Failed to open file: " << full_bed_level_filename.string() << std::endl;
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    status = bed->read(nx, ny);
    std::vector<double> zb_giv = bed->get_bed_level();

    double dx = x[ny] - x[0];
    double dy = y[1] - y[0];
    double Lx = dx * mesh2d->face[0]->dims[0];
    double Ly = dy * mesh2d->face[0]->dims[1];

    double dxinv = 1./dx;                               // invers grid size [m]
    double dyinv = 1./dy;                               // invers grid size [m]
    int nxny = nx * ny;                                   // total number of nodes
    double dxdy = dx * dy ;                               // area of control volume
    //if (viscosity == 0.0)
    //{
    //    viscosity = 0.2 * std::sqrt(dx*dx + dy*dy);
    //}

    int total_time_steps = int((tstop - tstart) / dt) + 1;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    double dtpseuinv = 0.0;                                     // Inverse of dtpseu
    int wrihis;                                           // write interval to his-file
    int wrimap;                                           // write interval to map-file
    if (stationary)
    {
        dt = 0.0;                                         // Time step size [s]
        dtinv = 0.0;                                      // stationary solution
        eps_bc_corr = 1.0;
        iter_max = 2 * iter_max;
        theta = 1.0;                                      // Stationary solution
        tstop = 1.;
        total_time_steps = 2;                             // initial step (step 1), stationary result (step 2)
        treg = 0.0;                                       // Thatcher-Harleman return time [s], when zero supply boundary value immediately
        wrihis = 1;                                       // write interval to his-file
        wrimap = 1;                                       // write interval to map-file
    }
    else
    {
        dtinv = 1. / dt;                                  // Inverse of dt [1/s]
        wrihis = std::max(int(dt * dtinv), int(dt_his * dtinv));      // write interval to his-file (every delta t)
        if (dt_map == 0.0)
        {
            wrimap = total_time_steps - 1;  // write only first and last time step
        }
        else
        {
            wrimap = std::max(int(dt * dtinv), int(dt_map * dtinv));     // write interval to map-file (every 1 sec , or every delta t)
        }
    }

    REGULARIZATION* regularization = new REGULARIZATION(iter_max, g);

    log_file << "=== Used input variables ==============================" << std::endl;
    log_file << "[Domain]" << std::endl;
    log_file << "mesh_file = " << grid_filename << std::endl;
    log_file << "bed_level_file = " << bed_level_filename << std::endl;

    log_file << std::endl << "[Time]" << std::endl;
    log_file << "tstart = " << tstart << std::endl;
    log_file << "tstop = " << tstop << std::endl;

    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "ini_vars = ";
    for (int i = 0; i < ini_vars.size(); ++i)
    {
        log_file << ini_vars[i];
        if (i < ini_vars.size() - 1) { log_file << ", "; }
    }
    log_file << std::endl;
    log_file << "Gauss_amp = " << gauss_amp << std::endl;
    log_file << "Gauss_mu = " << gauss_mu << std::endl;
    log_file << "Gauss_mu_x = " << gauss_mu_x << std::endl;
    log_file << "Gauss_mu_y = " << gauss_mu_y << std::endl;
    log_file << "Gauss_sigma = " << gauss_sigma << std::endl;
    log_file << "Gauss_sigma_x = " << gauss_sigma_x << std::endl;
    log_file << "Gauss_sigma_y = " << gauss_sigma_y << std::endl;

    log_file << std::endl << "[Boundary]" << std::endl;
    log_file << "bc_type = ";
    for (int i = 0; i < bc_type.size() - 1; ++i) { log_file << bc_type[i] << ", "; }
    log_file << bc_type[bc_type.size() - 1] << std::endl;
    log_file << "bc_vars = ";
    for (int i = 0; i < bc_vars.size() - 1; ++i) { log_file << bc_vars[i] << ", "; }
    log_file << bc_vars[bc_vars.size() - 1] << std::endl;
    log_file << "bc_vals = ";
    for (int i = 0; i < bc_vals.size() - 1; ++i) { log_file << bc_vals[i] << ", "; }
    log_file << bc_vals[bc_vals.size() - 1] << std::endl;
    log_file << "bc_absorbing = ";
    for (int i = 0; i < bc_absorbing.size() - 1; ++i) { log_file << bc_absorbing[i] << ", "; }
    log_file << bc_absorbing[bc_absorbing.size() - 1] << std::endl;
    log_file << "treg = " << treg << std::endl;
    log_file << "eps_bc_corr = " << eps_bc_corr << std::endl;

    log_file << std::endl << "[Physics]" << std::endl;
    log_file << "do_continuity = " << do_continuity << std::endl;
    log_file << "do_q_equation = " << do_q_equation << std::endl;
    log_file << "do_r_equation = " << do_r_equation << std::endl;
    log_file << "do_convection = " << do_convection << std::endl;
    log_file << "do_bed_shear_stress = " << do_bed_shear_stress << std::endl;
    log_file << "chezy_coefficient = " << chezy_coefficient << std::endl;
    log_file << "do_viscosity = " << do_viscosity << std::endl;
    log_file << "viscosity = " << visc_const << std::endl;

    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "dt = " << dt << std::endl;
    log_file << "dx = " << dx << std::endl;
    log_file << "dy = " << dy << std::endl;
    log_file << "theta = " << theta << std::endl;
    log_file << "c_psi = " << c_psi << std::endl;
    log_file << "iter_max = " << iter_max << std::endl;
    log_file << "eps_newton = " << eps_newton << std::endl;
    log_file << "eps_bicgstab = " << eps_bicgstab << std::endl;
    log_file << "eps_abs_function = " << eps_abs << std::endl;
    log_file << "regularization_init = " << regularization_init << std::endl;
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    log_file << "regularization_time = " << regularization_time << std::endl;

    log_file << std::endl << "[Output]" << std::endl;
    log_file << "dt_his = " << dt_his << std::endl;
    log_file << "dt_map = " << dt_map << std::endl;

    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix in x-direction

    std::vector<double> zb(nxny, 0.);                     // regularized bed level
    std::vector<double> hn(nxny, 0.);                     // water depth at (n)
    std::vector<double> qn(nxny, 0.);                     // q-flux at (n)
    std::vector<double> rn(nxny, 0.);                     // r-flux at (n)
    std::vector<double> hp(nxny, 0.);                     // total depth at (n+1,p), previous iteration
    std::vector<double> qp(nxny, 0.);                     // x-flow at (n+1,p), previous iteration
    std::vector<double> rp(nxny, 0.);                     // y-flow at (n+1,p), previous iteration
    std::vector<double> dh(nxny, 0.);                     // delta for water depth
    std::vector<double> dq(nxny, 0.);                     // delta for q-flux
    std::vector<double> dr(nxny, 0.);                     // delta for r-flux
    std::vector<double> s_giv(nxny, 0.);                     // water level, given
    std::vector<double> u_giv(nxny, 0.);                     // u-velocity, given
    std::vector<double> v_giv(nxny, 0.);                     // v-velocity, given
    std::vector<double> s(nxny, 0.);                     // water level, needed for post-processing
    std::vector<double> u(nxny, 0.);                     // u-velocity, needed for post-processing
    std::vector<double> v(nxny, 0.);                     // v-velocity, needed for post-processing
    std::vector<double> delta_h(nxny, 0.);
    std::vector<double> delta_q(nxny, 0.);
    std::vector<double> delta_r(nxny, 0.);
    std::vector<double> post_q;              // needed for postprocessing; bed shear stress; convection;
    std::vector<double> post_r;              // needed for postprocessing; bed shear stress; convection;
    if (do_bed_shear_stress || do_convection)
    {
        post_q.resize(nxny, 0.);              // needed for postprocessing; bed shear stress; convection;
        post_r.resize(nxny, 0.);              // needed for postprocessing; bed shear stress; convection;
    }
    std::vector<double> psi_11;  // needed when regularization of given initial function (ie bed level);
    std::vector<double> psi_22;  // needed when regularization of given initial function (ie bed level);
    std::vector<double> eq8;     // needed when regularization of given initial function (ie bed level);
    if (regularization_init)
    {
        psi_11.resize(nxny, 0.);  // needed when regularization of given initial function (ie bed level);
        psi_22.resize(nxny, 0.);  // needed when regularization of given initial function (ie bed level);
        eq8.resize(nxny, 0.);     // needed when regularization of given initial function (ie bed level);
    }
    std::vector<double> htheta(nxny);
    std::vector<double> qtheta(nxny);
    std::vector<double> rtheta(nxny);

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(3 * nxny, 3 * nxny);
    Eigen::VectorXd solution(3 * nxny);                     // solution vector [h, q, r]^{n}
    Eigen::VectorXd rhs(3 * nxny);                          // RHS vector [h, q, r]^{n}
    std::vector< Eigen::Triplet<double> > triplets; 
    triplets.reserve(9 * 3 * nxny); // 9-points per row, nrow = 3 * nxny; large enough to avoid reallocation
    
    status = build_matrix_pattern(triplets, nx, ny);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed(); // Very important for valuePtr() access
    
    solution.setZero(); // Delta h, Delta q and Delta r
    rhs.setZero();

    double cf = g / (chezy_coefficient * chezy_coefficient);  // bed friction coefficient
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
    w_ess[0] = w_nat[0];
    w_ess[1] = w_nat[1];
    w_ess[2] = w_nat[2];

    //initialize water level
    std::cout << "Initialisation" << std::endl;
    status = 0;
    double min_zb = *std::min_element(zb_giv.begin(), zb_giv.end());

    log_file << "-------------------------------------------------------" << std::endl;
    log_file << "Nodes   : " << nx << "x" << ny << "=" << nxny << std::endl;
    log_file << "Elements: " << (nx - 1) << "x" << (ny - 1) << "=" << (nx - 1) * (ny - 1) << std::endl;
    log_file << "Volumes : " << (nx - 2) << "x" << (ny - 2) << "=" << (nx - 2) * (ny - 2) << std::endl;
    log_file << "CFL (2D): " << std::sqrt(g * std::abs(min_zb)) * dt * std::sqrt(( 1./(dx*dx) + 1./(dy*dy))) << std::endl;
    log_file << "CFL (1D): " << std::sqrt(g * std::abs(min_zb)) * dt /dx << std::endl;
    log_file << "LxLy    : " << Lx - 2. * dx << "x" << Ly - 2. * dy << " without virtual cells" << std::endl;
    log_file << "dxdy    : " << dx << "x" << dy << "=" << dxdy << " [m2]" << std::endl;
    log_file << "nxny    : " << nx << "x" << ny << "=" << nxny << std::endl;
    log_file << "=======================================================" << std::endl;
    std::cout << "    LxLy: " << Lx - 2. * dx << "x" << Ly - 2. * dy << " without virtual cells" << std::endl;
    std::cout << "    dxdy: " << dx << "x" << dy << "=" << dxdy << " [m2]" << std::endl;
    std::cout << "    nxny: " << nx << "x" << ny << "=" << nxny << std::endl;
    std::cout << "======================================================" << std::endl;

    STOP_TIMER(Writing log-file);  // but two write statements are not timed
    if (status != 0) {
        log_file.close();
        exit(1);
    }

    (void) initial_conditions(x, y, nx, ny,
        s_giv, u_giv, v_giv, 
        ini_vars, gauss_amp, gauss_mu_x, gauss_mu_y, gauss_sigma_x, gauss_sigma_y);

    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        regularization->given_function(zb, psi_11, psi_22, eq8, zb_giv, nx, ny, dx, dy, c_psi, log_file);
        //regularization->given_function(visc_reg, psi, visc_giv, nx, ny, dx, dy, c_psi);
        STOP_TIMER(Regularization_init);
    }
    else
    {
        for (int i = 0; i < zb_giv.size(); ++i)
        {
            zb[i] = zb_giv[i];
        }
    }
    for (int k = 0; k < zb_giv.size(); ++k)
    {
        hn[k] = s_giv[k] - zb[k];  // Initial water depth
        qn[k] = hn[k] * u_giv[k];  // Initial q=hu -velocity
        rn[k] = hn[k] * v_giv[k];  // Initial r=hv -velocity
        // 
        hp[k] = hn[k]; 
        qp[k] = qn[k]; 
        rp[k] = rn[k]; 
    }

    for (int i = 0; i < nx*ny; ++i)
    {
        s[i] = hn[i] + zb[i];
        u[i] = qn[i] / hn[i];
        v[i] = rn[i] / hn[i];
    }

    double time = double(0) * dt;
    ////////////////////////////////////////////////////////////////////////////
    // Create map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    UGRID2D* map_file = new UGRID2D();
    status = map_file->open(nc_mapfilename, model_title);
    if (status != NC_NOERR)
    {
        std::cout << "Failed to open file (probably in use): \'" << map_filename << "\'" << std::endl;
        log_file  << "Failed to open file (probably in use): \'" << map_filename << "\'" << std::endl;
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    status = map_file->mesh2d();  // initialize the mesh2D variable

    int nr_nodes = nx * ny;
    int nr_edges = (nx - 1) * ny + nx * (ny - 1);
    int nr_faces = (nx - 1) * (ny - 1);
    int mesh2d_nmax_face_nodes = 4;  // all elements are quads
    status = map_file->def_dimensions(nr_nodes, nr_edges, nr_faces, mesh2d_nmax_face_nodes);

    std::vector<std::string> dim_names;
    dim_names.clear();
    dim_names.push_back("mesh2d_nEdges");
    dim_names.push_back("Two");
    status = map_file->add_mesh2d_edge_nodes("mesh2d_edge_nodes", dim_names, "Each edge connects two nodes");
    dim_names.clear();
    dim_names.push_back("mesh2d_nFaces");
    dim_names.push_back("mesh2d_nMax_face_nodes");
    status = map_file->add_mesh2d_edge_nodes("mesh2d_face_nodes", dim_names, "Each face contains four nodes");
    std::vector<int> mesh2d_edge_nodes;
    std::vector<int> mesh2d_face_nodes;
    std::vector<int> nodes_per_face;
    std::vector<int> node_mask(nr_nodes, 0);

    int p0 = 0;
    int p1 = 0;
    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j, ny);
            p1 = idx(i + 1, j, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i + 1, j    , ny);
            p1 = idx(i + 1, j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i + 1, j + 1, ny);
            p1 = idx(i    , j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i    , j + 1, ny);
            p1 = idx(i    , j    , ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                //node_mask[p1] += 1;
            }
        }
    }
    status = map_file->put_variable_2("mesh2d_edge_nodes", mesh2d_edge_nodes);
    
    int p2;
    int p3;
    double fill_value = mesh2d->node[0]->fill_value;
    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            if ((x[p0] != fill_value || y[p0] != fill_value) &&
                (x[p1] != fill_value || y[p1] != fill_value) &&
                (x[p2] != fill_value || y[p2] != fill_value) &&
                (x[p3] != fill_value || y[p3] != fill_value) )
            {
                mesh2d_face_nodes.push_back(p0);
                mesh2d_face_nodes.push_back(p1);
                mesh2d_face_nodes.push_back(p2);
                mesh2d_face_nodes.push_back(p3);
            }
        }
    }
    status = map_file->put_variable_4("mesh2d_face_nodes", mesh2d_face_nodes);

    fill_value = mesh2d->node[0]->fill_value;
    dim_names.clear();
    dim_names.push_back("mesh2d_nNodes");
    status = map_file->add_variable("mesh2d_node_x", dim_names, "projection_x_coordinate", "x", "m");
    status = map_file->add_attribute("mesh2d_node_x", "_FillValue", fill_value);
    status = map_file->add_variable("mesh2d_node_y", dim_names, "projection_y_coordinate", "y", "m");
    status = map_file->add_attribute("mesh2d_node_y", "_FillValue", fill_value);
    status = map_file->put_variable("mesh2d_node_x", x);
    status = map_file->put_variable("mesh2d_node_y", y);

    // Compute edges centres
    std::vector<double> edge_x;
    std::vector<double> edge_y;

    for (int i = 0; i < mesh2d_edge_nodes.size(); i += 2)
    {
        p0 = mesh2d_edge_nodes[i];
        p1 = mesh2d_edge_nodes[i+1];
        edge_x.push_back(0.5 * (x[p0] + x[p1]));
        edge_y.push_back(0.5 * (y[p0] + y[p1]));
    }
    dim_names.clear();
    dim_names.push_back("mesh2d_nEdges");
    status = map_file->add_variable("mesh2d_edge_x", dim_names, "projection_x_coordinate", "x", "m");
    status = map_file->add_attribute("mesh2d_edge_x", "_FillValue", fill_value);
    status = map_file->add_variable("mesh2d_edge_y", dim_names, "projection_y_coordinate", "y", "m");
    status = map_file->add_attribute("mesh2d_edge_y", "_FillValue", fill_value);
    status = map_file->put_variable("mesh2d_edge_x", edge_x);
    status = map_file->put_variable("mesh2d_edge_y", edge_y);

    // Compute mass centres of faces
    std::vector<double> xmc;
    std::vector<double> ymc;
    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            xmc.push_back(0.25 * (x[p0] + x[p1] + x[p2] + x[p3]));
            ymc.push_back(0.25 * (y[p0] + y[p1] + y[p2] + y[p3]));
        }
    }
    fill_value = mesh2d->face[0]->fill_value;
    dim_names.clear();
    dim_names.push_back("mesh2d_nFaces");
    status = map_file->add_variable("mesh2d_face_x", dim_names, "projection_x_coordinate", "x", "m");
    status = map_file->add_attribute("mesh2d_face_x", "_FillValue", fill_value);
    status = map_file->add_variable("mesh2d_face_y", dim_names, "projection_y_coordinate", "y", "m");
    status = map_file->add_attribute("mesh2d_face_y", "_FillValue", fill_value);
    status = map_file->put_variable("mesh2d_face_x", xmc);
    status = map_file->put_variable("mesh2d_face_y", ymc);

    // Compute area of faces
    std::vector<double> cell_area;
    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            double area = std::abs(0.5 * ((x[p0] * y[p1] + x[p1] * y[p2] + x[p2] * y[p3] + x[p3] * y[p0]) - (y[p0] * x[p1] + y[p1] * x[p2] + y[p2] * x[p3] + y[p3] * x[p0])));
            cell_area.push_back(area);
        }
    }
    dim_names.clear();
    dim_names.push_back("mesh2d_nFaces");
    status = map_file->add_variable("cell_area", dim_names, "cell_area", "-", "m2", "mesh2D", "face");
    status = map_file->add_attribute("cell_area", "coordinates", "mesh2d_face_x, mesh2d_face_y");
    status = map_file->put_variable("cell_area", cell_area);

    status = map_file->add_time_series();

    dim_names.clear();
    dim_names.push_back("time");
    dim_names.push_back("mesh2d_nNodes");
    std::string map_h_name("hn_2d");
    std::string map_q_name("qn_2d");
    std::string map_r_name("rn_2d");
    std::string map_dh_name("delta_h_2d");
    std::string map_dq_name("delta_q_2d");
    std::string map_dr_name("delta_r_2d");
    std::string map_s_name("s_2d");
    std::string map_u_name("u_2d");
    std::string map_v_name("v_2d");
    std::string map_zb_name("zb_2d");
    std::string map_psi_11_name("psi_11");
    std::string map_psi_22_name("psi_22");
    std::string map_eq8_name("eq8");
    std::string map_beds_q_name("bed_stress_q");
    std::string map_beds_r_name("bed_stress_r");
    std::string map_conv_q_name("convection_q");
    std::string map_conv_r_name("convection_r");

    status = map_file->add_variable(map_h_name, dim_names, "sea_floor_depth_below_sea_surface", "Water depth", "m", "mesh2D", "node");
    status = map_file->add_variable(map_q_name, dim_names, "", "Water flux (x)", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_r_name, dim_names, "", "Water flux (y)", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dh_name, dim_names, "", "Delta h^{n+1,p+1}", "m", "mesh2D", "node");
    status = map_file->add_variable(map_dq_name, dim_names, "", "Delta q^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dr_name, dim_names, "", "Delta r^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_s_name, dim_names, "sea_surface_height_above_geoid", "WaterLevel", "m", "mesh2D", "node");
    status = map_file->add_variable(map_u_name, dim_names, "sea_water_x_velocity", "Velocity (x)", "m s-1", "mesh2D", "node");
    status = map_file->add_variable(map_v_name, dim_names, "sea_water_y_velocity", "Velocity (y)", "m s-1", "mesh2D", "node");
    status = map_file->add_variable(map_zb_name, dim_names, "", "BedLevel", "m", "mesh2D", "node");
    if (regularization_init)
    {
        status = map_file->add_variable(map_psi_11_name, dim_names, "", "Psi_11", "m2 s-1", "mesh2D", "node");
        status = map_file->add_variable(map_psi_22_name, dim_names, "", "Psi_22", "m2 s-1", "mesh2D", "node");
        status = map_file->add_variable(map_eq8_name, dim_names, "", "Eq8", "-", "mesh2D", "node");
    }
    if (do_bed_shear_stress)
    {
        status = map_file->add_variable(map_beds_q_name, dim_names, "", "Bed shear stress (x)", "m2 s-2", "mesh2D", "node");
        status = map_file->add_variable(map_beds_r_name, dim_names, "", "Bed shear stress (y)", "m2 s-2", "mesh2D", "node");
    }
    if (do_convection)
    {
        status = map_file->add_variable(map_conv_q_name, dim_names, "", "Convection (x)", "m2 s-2", "mesh2D", "node");
        status = map_file->add_variable(map_conv_r_name, dim_names, "", "Convection (y)", "m2 s-2", "mesh2D", "node");
    }

    // Put data on map file
    START_TIMER(Writing map-file);
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    map_file->put_time_variable(map_h_name, nst_map, hn);
    map_file->put_time_variable(map_q_name, nst_map, qn);
    map_file->put_time_variable(map_r_name, nst_map, rn);
    map_file->put_time_variable(map_dh_name, nst_map, dh);
    map_file->put_time_variable(map_dq_name, nst_map, dq);
    map_file->put_time_variable(map_dr_name, nst_map, dr);
    map_file->put_time_variable(map_s_name, nst_map, s);
    map_file->put_time_variable(map_u_name, nst_map, u);
    map_file->put_time_variable(map_v_name, nst_map, v);
    map_file->put_time_variable(map_zb_name, nst_map, zb_giv);
    if (regularization_init)
    {
        map_file->put_time_variable(map_psi_11_name, nst_map, psi_11);
        map_file->put_time_variable(map_psi_22_name, nst_map, psi_22);
        map_file->put_time_variable(map_eq8_name, nst_map, eq8);
    }
    if (do_bed_shear_stress)
    {
        bed_shear_stress_rhs(post_q, post_r, hn, qn, rn, cf, nx, ny);
        map_file->put_time_variable(map_beds_q_name, nst_map, post_q);
        map_file->put_time_variable(map_beds_r_name, nst_map, post_r);
    }
    if (do_convection)
    {
        convection_rhs(post_q, post_r, hn, qn, rn, dx, dy, nx, ny);
        map_file->put_time_variable(map_conv_q_name, nst_map, post_q);
        map_file->put_time_variable(map_conv_r_name, nst_map, post_r);
    }
    STOP_TIMER(Writing map-file);

    // End define map file
    ////////////////////////////////////////////////////////////////////////////
    // Create time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations
    int centre = idx(nx / 2, ny / 2, ny);
    int w_bnd  = idx(1     , ny / 2, ny);  // at boundary, not at virtual point
    int e_bnd  = idx(nx - 2, ny / 2, ny);  // at boundary, not at virtual point
    int s_bnd  = idx(nx / 2, 1     , ny);
    int n_bnd  = idx(nx / 2, ny - 2, ny);

    int sw_bnd = idx(1     , 1     , ny);
    int ne_bnd = idx(nx - 2, ny - 2, ny);
    int nw_bnd = idx(1     , ny - 2, ny);
    int se_bnd = idx(nx - 2, 1     , ny);

    //Station halfway to boundary

    int hw_bnd  = idx(    nx / 4 + 1,     ny / 2    , ny);
    int he_bnd  = idx(3 * nx / 4 - 1,     ny / 2    , ny);
    int hs_bnd  = idx(    nx / 2    ,     ny / 4 + 1, ny);
    int hn_bnd  = idx(    nx / 2    , 3 * ny / 4 - 1, ny);

    int hsw_bnd = idx(    nx / 4 + 1,     ny / 4 + 1, ny);
    int hne_bnd = idx(3 * nx / 4 - 1, 3 * ny / 4 - 1, ny);
    int hnw_bnd = idx(    nx / 4 + 1, 3 * ny / 4 - 1, ny);
    int hse_bnd = idx(3 * nx / 4 - 1,     ny / 4 + 1, ny);

    int p_a;
    int p_b;
    int p_c;
    int p_d;

    if (int(2500. / dx) * dx != 2500. || int(2500. / dy) * dy != 2500.)
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "dx=" << dx << " or dy=" << dy << " is not a divider of 2500 [m]" << std::endl;
        std::cout << "Press Enter to finish";
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    double x_a = std::min(Lx / 2., 2500.);
    double x_b = std::min(Lx / 2., 2000.);
    double x_c = std::min(Lx / 2., 1500.);
    double x_d = std::min(Lx / 2., 0.);
    double y_a = std::min(Ly / 2., 0.);
    double y_b = std::min(Ly / 2., 1500.);
    double y_c = std::min(Ly / 2., 2000.);
    double y_d = std::min(Ly / 2., 2500.);
    p_a = idx(int((y_a / dy) + nx / 2), int((x_a / dx) + ny / 2), ny);
    p_b = idx(int((y_b / dy) + nx / 2), int((x_b / dx) + ny / 2), ny);
    p_c = idx(int((y_c / dy) + nx / 2), int((x_c / dx) + ny / 2), ny);
    p_d = idx(int((y_d / dy) + nx / 2), int((x_d / dx) + ny / 2), ny);

    std::vector<double> x_obs = { x[p_a], x[p_b], x[p_c], x[p_d], x[centre], 
        x[n_bnd], x[ne_bnd], x[e_bnd], x[se_bnd], x[s_bnd], x[sw_bnd], x[w_bnd], x[nw_bnd],
        x[hn_bnd], x[hne_bnd], x[he_bnd], x[hse_bnd], x[hs_bnd], x[hsw_bnd], x[hw_bnd], x[hnw_bnd]
    };
    std::vector<double> y_obs = { y[p_a], y[p_b], y[p_c], y[p_d], y[centre], 
        y[n_bnd], y[ne_bnd], y[e_bnd], y[se_bnd], y[s_bnd], y[sw_bnd], y[w_bnd], y[nw_bnd],
        y[hn_bnd], y[hne_bnd], y[he_bnd], y[hse_bnd], y[hs_bnd], y[hsw_bnd], y[hw_bnd], y[hnw_bnd],
    };
    int nsig = 0;
    for (int i = 0; i < x_obs.size(); ++i)
    {
        nsig = std::max(nsig, (int)std::log10(x_obs[i]));
    }
    //nsig += 1;

    std::vector<std::string> obs_stations;
    obs_stations.push_back(setup_obs_name(x[p_a], y[p_a], nsig, "A"));
    obs_stations.push_back(setup_obs_name(x[p_b], y[p_b], nsig, "B"));
    obs_stations.push_back(setup_obs_name(x[p_c], y[p_c], nsig, "C"));
    obs_stations.push_back(setup_obs_name(x[p_d], y[p_d], nsig, "D"));
    obs_stations.push_back(setup_obs_name(x[centre], y[centre], nsig, "Centre"));
    obs_stations.push_back(setup_obs_name(x[n_bnd] , y[n_bnd] , nsig, "N"));
    obs_stations.push_back(setup_obs_name(x[ne_bnd], y[ne_bnd], nsig, "NE"));
    obs_stations.push_back(setup_obs_name(x[e_bnd] , y[e_bnd] , nsig, "E"));
    obs_stations.push_back(setup_obs_name(x[se_bnd], y[se_bnd], nsig, "SE"));
    obs_stations.push_back(setup_obs_name(x[s_bnd] , y[s_bnd] , nsig, "S"));
    obs_stations.push_back(setup_obs_name(x[sw_bnd], y[sw_bnd], nsig, "SW"));
    obs_stations.push_back(setup_obs_name(x[w_bnd] , y[w_bnd] , nsig, "W"));
    obs_stations.push_back(setup_obs_name(x[nw_bnd], y[nw_bnd], nsig, "NW"));

    obs_stations.push_back(setup_obs_name(x[hn_bnd] , y[hn_bnd] , nsig, "Halfway N"));
    obs_stations.push_back(setup_obs_name(x[hne_bnd], y[hne_bnd], nsig, "Halfway NE"));
    obs_stations.push_back(setup_obs_name(x[he_bnd] , y[he_bnd] , nsig, "Halfway E"));
    obs_stations.push_back(setup_obs_name(x[hse_bnd], y[hse_bnd], nsig, "Halfway SE"));
    obs_stations.push_back(setup_obs_name(x[hs_bnd] , y[hs_bnd] , nsig, "Halfway S"));
    obs_stations.push_back(setup_obs_name(x[hsw_bnd], y[hsw_bnd], nsig, "Halfway SW"));
    obs_stations.push_back(setup_obs_name(x[hw_bnd] , y[hw_bnd] , nsig, "Halfway W"));
    obs_stations.push_back(setup_obs_name(x[hnw_bnd], y[hnw_bnd], nsig, "Halfway NW"));

    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_h_name("hn_2d");
    std::string his_q_name("qn_2d");
    std::string his_r_name("rn_2d");
    std::string his_s_name("water_level");
    std::string his_u_name("u_velocity");
    std::string his_v_name("v_velocity");

    his_file->add_variable(his_h_name, "sea_floor_depth_below_sea_surface", "Water depth", "m");
    his_file->add_variable(his_q_name, "", "Water flux (x)", "m2 s-1");
    his_file->add_variable(his_r_name, "", "Water flux (y)", "m2 s-1");
    his_file->add_variable(his_s_name, "sea_surface_height", "Water level", "m");
    his_file->add_variable(his_u_name, "sea_water_x_velocity", "Velocity xdir", "m s-1");
    his_file->add_variable(his_v_name, "sea_water_y_velocity", "Velocity ydir", "m s-1");

    // Put data on time history file
    START_TIMER(Writing his-file);
    int nst_his = 0;
    his_file->put_time(nst_his, time);
    std::vector<double> his_values = { hn[p_a], hn[p_b], hn[p_c], hn[p_d], hn[centre], 
        hn[n_bnd], hn[ne_bnd], hn[e_bnd], hn[se_bnd], hn[s_bnd], hn[sw_bnd], hn[w_bnd], hn[nw_bnd],
        hn[hn_bnd], hn[hne_bnd], hn[he_bnd], hn[hse_bnd], hn[hs_bnd], hn[hsw_bnd], hn[hw_bnd], hn[hnw_bnd]
    };
    his_file->put_variable(his_h_name, nst_his, his_values);

    his_values = { qn[p_a], qn[p_b], qn[p_c], qn[p_d], qn[centre], 
        qn[n_bnd], qn[ne_bnd], qn[e_bnd], qn[se_bnd], qn[s_bnd], qn[sw_bnd], qn[w_bnd], qn[nw_bnd],
        qn[hn_bnd], qn[hne_bnd], qn[he_bnd], qn[hse_bnd], qn[hs_bnd], qn[hsw_bnd], qn[hw_bnd], qn[hnw_bnd]
    };
    his_file->put_variable(his_q_name, nst_his, his_values);
 
    his_values = { rn[p_a], rn[p_b], rn[p_c], rn[p_d], rn[centre], 
        rn[n_bnd], rn[ne_bnd], rn[e_bnd], rn[se_bnd], rn[s_bnd], rn[sw_bnd], rn[w_bnd], rn[nw_bnd],
        rn[hn_bnd], rn[hne_bnd], rn[he_bnd], rn[hse_bnd], rn[hs_bnd], rn[hsw_bnd], rn[hw_bnd], rn[hnw_bnd]
    };
    his_file->put_variable(his_r_name, nst_his, his_values);

    his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], 
        s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd],
        s[hn_bnd], s[hne_bnd], s[he_bnd], s[hse_bnd], s[hs_bnd], s[hsw_bnd], s[hw_bnd], s[hnw_bnd]
    };
    his_file->put_variable(his_s_name, nst_his, his_values);

    his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], 
        u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd],
        u[hn_bnd], u[hne_bnd], u[he_bnd], u[hse_bnd], u[hs_bnd], u[hsw_bnd], u[hw_bnd], u[hnw_bnd]
    };
    his_file->put_variable(his_u_name, nst_his, his_values);

    his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], 
        v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd],
        v[hn_bnd], v[hne_bnd], v[he_bnd], v[hse_bnd], v[hs_bnd], v[hsw_bnd], v[hw_bnd], v[hnw_bnd]
    };
    his_file->put_variable(his_v_name, nst_his, his_values);

    std::string his_newton_iter_name("his_newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "iterations", "Newton iteration", "-");
    his_values = { 0.0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    std::string his_BiCGstab_iter_name("his_BiCGstab_iterations");
    his_file->add_variable_without_location(his_BiCGstab_iter_name, "iterations", "BiCGstab iteration", "-");
    his_values = { 0.0 };
    his_file->put_variable(his_BiCGstab_iter_name, nst_his, his_values);

    std::string his_BiCGstab_iter_error_name("BiCGstab_iteration_error");
    his_file->add_variable_without_location(his_BiCGstab_iter_error_name, "iteration_error", "BiCGstab iteration error", "-");
    his_values = { 0.0 };
    his_file->put_variable(his_BiCGstab_iter_error_name, nst_his, his_values);
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
    STOP_TIMER(Initialization);
    
    std::cout << "Start time-loop" << std::endl;
    if (stationary) { std::cout << "Stationary solution" << std::endl; }
    else { std::cout << "Time dependent simulation" << std::endl; }
    std::cout << std::fixed << std::setprecision(3) << "tstart= " << tstart + time << ";   tstop= " << tstart + tstop << ";   dt= " << dt << std::endl;

    // Start time loop
    START_TIMER(Time loop);
    double dh_max = 0.0;
    double dq_max = 0.0;
    double dr_max = 0.0;
    int dh_maxi = 0;
    int dq_maxi = 0;
    int dr_maxi = 0;

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
    Eigen::IncompleteLUT<double> ilu;
    //Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    //Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > solver;

    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        time = dt * double(nst);

        int select = 1;  // Constant boundary condition
        std::vector<double> bc(4, 0.0);
        boundary_condition(bc[BC_NORTH], bc_vals[BC_NORTH], time, treg, select);
        boundary_condition(bc[BC_EAST ], bc_vals[BC_EAST ], time, treg, select);
        boundary_condition(bc[BC_SOUTH], bc_vals[BC_SOUTH], time, treg, select);
        boundary_condition(bc[BC_WEST ], bc_vals[BC_WEST ], time, treg, select);

        if (logging == "pattern")
        {
            log_file << "=== Matrix build matrix pattern =======================" << std::endl;
            for (int i = 0; i < 3 * nxny; ++i)
            {
                for (int j = 0; j < 3*nxny; ++j)
                {
                    if (A.coeff(i, j) != 0.0)
                    {
                        log_file << "* ";
                    }
                    else
                    {
                        log_file << "- ";
                        //log_file << std::showpos << std::setprecision(3) << std::scientific << A.coeff(i, j) << " ";
                    }
                    if (std::fmod(j+1,3) == 0) { log_file << "| "; }
                }
                log_file << std::endl;
                if (std::fmod(i+1,3) == 0) { log_file << std::endl; }
            }
        }
        int used_newton_iter = 0;
        int used_lin_solv_iter = 0;
        START_TIMER(Newton iteration);
        for (int iter = 0; iter < iter_max; ++iter)
        {
            used_newton_iter += 1;
            if (logging == "iteration")
            {
                log_file << "Iteration: " << used_newton_iter << std::endl;
            }
            //
            // interior nodes
            //
            START_TIMER(Linear wave);
            for (int k = 0; k < nxny; ++k)
            {
                htheta[k] = theta * hp[k] + (1.0 - theta) * hn[k];
                qtheta[k] = theta * qp[k] + (1.0 - theta) * qn[k];
                rtheta[k] = theta * rp[k] + (1.0 - theta) * rn[k];
            }

            double* values = A.valuePtr();         // pointer to all non-zero values
            const int* outer = A.outerIndexPtr();   // row start pointers

            //const int* inner = A.innerIndexPtr();   // column indices
            //int ncol = A.cols();  // ==A.outerSize()
            //int nrow = A.rows();  // ==A.outerSize()
            //int nr_nodes = A.cols()/3;

            // row 0, 1, 2; sw-corner
            // row 3,..., 3 * (ny - 1) - 1: west boundary
            // row 3 * (ny - 1), +1, +2; nw-corner
            // row std::fmod(col , 3 * ny), +1, +2; south boundary 
            // row std::fmod(col + 3, 3 * ny), +1, +2; north boundary
            // row 3 * (nx - 1) * ny, +1, +2; se-corner
            // row 3 * (nx - 1) * ny + 3, ..., 3 * nx * ny - 3 - 1: east boundary
            // row 3 * nx * ny - 3, +1, +2 : ne-corner

           // south-west corner
            for (int row = 0; row < 3; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status =  corner_south_west(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // west boundary
            for (int row = 3; row < 3 * (ny - 1); row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status =  boundary_west(values, row, c_eq, q_eq, r_eq, rhs, 
                                        dtinv, dxinv, theta, g, eps_bc_corr, 
                                        stationary, do_convection, nx, ny,
                                        hn, qn, rn,
                                        hp, qp, rp,
                                        htheta, qtheta, rtheta,
                                        zb, bc_type,  bc_vars, BC_WEST, bc,
                                        w_nat, w_ess);
            }
            // north-west corner
            for (int row = 3 * (ny - 1); row < 3 * ny; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                 status =  corner_north_west(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // interior with south and north boundary
            for (int row = 3 * ny; row < 3 * (nx - 1) * ny; row += 3) 
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status =  interior(values, row, c_eq, q_eq, r_eq, rhs, 
                                    dtinv, dxinv, theta, g, do_convection, nx, ny,
                                    hn, qn, rn,
                                    hp, qp, rp,
                                    htheta, qtheta, rtheta,
                                    zb, dx, dy, dxdy, mass);

                if (std::fmod(row , 3 * ny) == 0) {
                    // south boundary, over write coefficients
                    status = boundary_south(values, row, c_eq, q_eq, r_eq, rhs, 
                                            dtinv, dxinv, theta, g, eps_bc_corr, 
                                            stationary, do_convection, nx, ny,
                                            hn, qn, rn,
                                            hp, qp, rp,
                                            htheta, qtheta, rtheta,
                                            zb, bc_type,  bc_vars, BC_SOUTH, bc,
                                            w_nat, w_ess);
                } 
                if (std::fmod(row + 3, 3 * ny) == 0) {
                    // north boundary, over write coefficients
                    status = boundary_north(values, row, c_eq, q_eq, r_eq, rhs, 
                                            dtinv, dxinv, theta, g, eps_bc_corr, 
                                            stationary, do_convection, nx, ny,
                                            hn, qn, rn,
                                            hp, qp, rp,
                                            htheta, qtheta, rtheta,
                                            zb, bc_type,  bc_vars, BC_NORTH, bc,
                                            w_nat, w_ess);
                }
            }
            // south-east corner
            for (int row = 3 * (nx - 1) * ny; row < 3 * (nx - 1) * ny + 3; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status = corner_south_east(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // east boundary
            for (int row = 3 * (nx - 1) * ny + 3; row < 3 * nx * ny - 3; row += 3) 
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status = boundary_east(values, row, c_eq, q_eq, r_eq, rhs, 
                                        dtinv, dxinv, theta, g, eps_bc_corr, 
                                        stationary, do_convection, nx, ny,
                                        hn, qn, rn,
                                        hp, qp, rp,
                                        htheta, qtheta, rtheta,
                                        zb, bc_type,  bc_vars, BC_EAST, bc,
                                        w_nat, w_ess);
            }
            // north-east corner
            for (int row = 3 * nx * ny - 3; row < 3 * nx * ny; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status =  corner_north_east(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            STOP_TIMER(Linear wave);

            if (do_bed_shear_stress)
            {
                // For the moment only interior nodes (2025-08-13)
                // Do not clear the rows, but add the bed shear stress terms

                START_TIMER(Bed shear stress);
                // corner_south_west
                // boundary_west
                // corner_north_west

                // interior with south and north boundary
                for (int row = 3 * ny; row < 3 * (nx - 1) * ny; row += 3) 
                {
                    int c_eq = outer[row    ];
                    int q_eq = outer[row + 1];
                    int r_eq = outer[row + 2];

                    status = bed_shear_stress_matrix_rhs(values, row, c_eq, q_eq, r_eq, rhs,
                                htheta, qtheta, rtheta, cf, theta, dx, dy, nx, ny);
                    // boundary_south
                    // boundray_north
                }
                // corner_south_east
                // boundary_east
                // corner_north_east

                STOP_TIMER(Bed shear stress);
            }
            if (do_convection)
            {
                // For the moment only interior nodes (2025-08-13)
                // Do not clear the rows, but add the bed shear stress terms

                START_TIMER(Convection);
                // corner_south_west
                // boundary_west
                // corner_north_west

                // interior with south and north boundary
                for (int row = 3 * ny; row < 3 * (nx - 1) * ny; row += 3) 
                {
                    int c_eq = outer[row    ];
                    int q_eq = outer[row + 1];
                    int r_eq = outer[row + 2];

                    status = convection_matrix_rhs(values, row, c_eq, q_eq, r_eq, rhs,
                                htheta, qtheta, rtheta, theta, dx, dy, nx, ny);
                    // boundary_south
                    // boundray_north
                }
                // corner_south_east
                // boundary_east
                // corner_north_east

                STOP_TIMER(Convection);
            }
            if (do_viscosity)
            {
                START_TIMER(Viscosity);
                //
                // viscosity
                //
                STOP_TIMER(Viscosity);
            }
            
            if (nst == 1 && iter == 0)
            {
                START_TIMER(BiCGStab_initialization);
            }
            else if (nst != 1)
            {
                START_TIMER(BiCGStab);
            }

            if (logging == "matrix" && (nst == 1 || nst == total_time_steps-1) && iter == 0)
            {
                log_file << "=== Matrix ============================================" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    for (int j = 0; j < 3*nxny; ++j)
                    {
                        log_file << std::showpos << std::setprecision(3) << std::scientific << A.coeff(i, j) << " ";
                        if (std::fmod(j+1,3) == 0) { log_file << "| "; }
                    }
                    log_file << std::endl;
                    if (std::fmod(i+1,3) == 0) { log_file << std::endl; }
                }
                log_file << "=== Diagonal dominant == diag, off_diag, -, + =========" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    double off_diag = 0.0;
                    double diag = 0.0;
                    for (int j = 0; j < 3*nxny; ++j)
                    {
                        if (i != j ) { off_diag += std::abs(A.coeff(i,j)); }
                        else { diag = std::abs(A.coeff(i,j)); }
                    }
                    log_file << std::showpos << std::setprecision(5) << std::scientific << diag << " " << off_diag << " " 
                        <<  diag - off_diag << " " <<  diag + off_diag << " ";
                    log_file << std::endl;
                    if (std::fmod(i+1,3) == 0) { log_file << std::endl; }
                }
                log_file << "=== RHS ===============================================" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
                    if (std::fmod(i+1,3) == 0) { log_file << std::endl; }
                }
            }

            if (nst == 1 && iter == 0)
            {
                ilu.setDroptol(1e-02);
                ilu.setFillfactor(9);
                //ilu.compute(A);
            }
            solver.compute(A);
            solver.setTolerance(eps_bicgstab);
            solution = solver.solve(rhs);
            //solution = solver.solveWithGuess(rhs, solution);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(BiCGStab_initialization);
            }
            else if (nst != 1)
            {
                STOP_TIMER(BiCGStab);
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(4) << std::scientific << time
                    << "    Newton iteration: " << used_newton_iter
                    << "    BiCGstab iterations: " << solver.iterations()
                    << "    estimated error: " << solver.error()
                    << std::endl;
            }
            // 
            // The new solution is the previous iterant plus the delta
            //
            dh_max = 0.0;
            dq_max = 0.0;
            dr_max = 0.0;
            dh_maxi = 0;
            dq_maxi = 0;
            dr_maxi = 0;
            int k = 0;
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    k = idx(i, j ,ny);
                    // needed for postprocessing
                    hp[k] += solution[3 * k];
                    qp[k] += solution[3 * k + 1];
                    rp[k] += solution[3 * k + 2];
                    delta_h[k] = solution[3 * k];  // h, continuity-eq
                    delta_q[k] = solution[3 * k + 1];  // q, momentum-eq
                    delta_r[k] = solution[3 * k + 2];  // r, momentum-eq
                    if (dh_max < std::abs(delta_h[k]))
                    {
                        dh_max = std::abs(delta_h[k]);
                        dh_maxi = k;
                    }
                    if (dq_max < std::abs(delta_q[k]))
                    {
                        dq_max = std::abs(delta_q[k]);
                        dq_maxi = k;
                    }
                    if (dr_max < std::abs(delta_r[k]))
                    {
                        dr_max = std::abs(delta_r[k]);
                        dr_maxi = k;
                    }
                }
            }
            if (logging == "matrix" && (nst == 1 || nst == total_time_steps-1) && iter == 0)
            {
                log_file << "=== Solution ==========================================" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << solution[i] << std::endl;
                }
                log_file << "=== hp, qp, rp, zeta ==================================" << std::endl;
                for (int i = 0; i < nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << hp[i] << ", " << qp[i] << ", " << rp[i] << ", " << hp[i] + zb[i] << std::endl;
                }
                log_file << "=======================================================" << std::endl;
            }
            used_lin_solv_iter = std::max(used_lin_solv_iter, (int)solver.iterations());
            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                //for (int i = 0; i < nr_nodes; ++i)
                //{
                //    u[i] = qp[i] / hp[i];
                //}
                //(void)regular->first_derivative(psi, visc_reg, u, dx);
                STOP_TIMER(Regularization_iter_loop);
                //for (int i = 0; i < nr_nodes; ++i)
                //{
                    //visc[i] = visc_reg[i] * std::abs(psi[i]);
                    //pe[i] = qp[i] / hp[i] * dx / visc[i];
                //}
            }
            if (dh_max < eps_newton && dq_max < eps_newton && dr_max < eps_newton)
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
                << "    Delta q^{n + 1,p + 1}: " << dq_max << " at index: " << dq_maxi << std::endl
                << "    Delta r^{n + 1,p + 1}: " << dr_max << " at index: " << dr_maxi << std::endl;
        }
        else
        {
            if (std::fmod(time, 60.) == 0)
            {
                std::cout << std::fixed << std::setprecision(2) << tstart + time << ";   " << tstart + tstop << std::endl;
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(4) << std::scientific << time
                    << std::setprecision(8) << std::scientific << std::endl
                    << "    Newton iterations  : " << used_newton_iter << std::endl
                    << "    Delta h^{n + 1,p + 1}: " << dh_max << " at index: " << dh_maxi << std::endl
                    << "    Delta q^{n + 1,p + 1}: " << dq_max << " at index: " << dq_maxi << std::endl
                    << "    Delta r^{n + 1,p + 1}: " << dr_max << " at index: " << dr_maxi << std::endl;

            }
        }
        if (used_newton_iter == iter_max)
        {
            if (dh_max > eps_newton || dq_max > eps_newton || dr_max > eps_newton)
            {
                log_file << "    ----    maximum number of iterations reached, probably not converged, at time: " <<  time << " [sec]" << std::endl;
            }
        }
        for (int k = 0; k < nxny; ++k)
        {
            hn[k] = hp[k];  // h, continuity-eq
            qn[k] = qp[k];  // q, momentum-eq
            rn[k] = rp[k];  // r, momentum-eq
            s[k] = hn[k] + zb[k];
            u[k] = qn[k] / hn[k];
            v[k] = rn[k] / hn[k];
        }

        //START_TIMER(Regularization_time_loop);
        //(void)regular->first_derivative(psi, visc_reg, u, dx);
        //STOP_TIMER(Regularization_time_loop);
        //if (regularization_time)
        //{
        //    for (int i = 0; i < nr_nodes; ++i)
        //    {
        //        visc[i] = visc_reg[i] * std::abs(psi[i]);
        //        pe[i] = u[i] * dx / visc[i];
        //    }
        //}


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
            map_file->put_time_variable(map_h_name, nst_map, hn);
            map_file->put_time_variable(map_q_name, nst_map, qn);
            map_file->put_time_variable(map_r_name, nst_map, rn);
            map_file->put_time_variable(map_dh_name, nst_map, delta_h);
            map_file->put_time_variable(map_dq_name, nst_map, delta_q);
            map_file->put_time_variable(map_dr_name, nst_map, delta_r);
            map_file->put_time_variable(map_s_name, nst_map, s);
            map_file->put_time_variable(map_u_name, nst_map, u);
            map_file->put_time_variable(map_v_name, nst_map, v);
            map_file->put_time_variable(map_zb_name, nst_map, zb);
            if (regularization_init)
            {
                map_file->put_time_variable(map_psi_11_name, nst_map, psi_11);
                map_file->put_time_variable(map_psi_22_name, nst_map, psi_22);
                map_file->put_time_variable(map_eq8_name, nst_map, eq8);
            }
            if (do_bed_shear_stress)
            {
                bed_shear_stress_rhs(post_q, post_r, hn, qn, rn, cf, nx, ny);
                map_file->put_time_variable(map_beds_q_name, nst_map, post_q);
                map_file->put_time_variable(map_beds_r_name, nst_map, post_r);
            }
            if (do_convection)
            {
                convection_rhs(post_q, post_r, hn, qn, rn, dx, dy, nx, ny);
                map_file->put_time_variable(map_conv_q_name, nst_map, post_q);
                map_file->put_time_variable(map_conv_r_name, nst_map, post_r);
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
            his_file->put_time(nst_his, time);

            his_values = { hn[p_a], hn[p_b], hn[p_c], hn[p_d], hn[centre], 
                hn[n_bnd], hn[ne_bnd], hn[e_bnd], hn[se_bnd], hn[s_bnd], hn[sw_bnd], hn[w_bnd], hn[nw_bnd],
                hn[hn_bnd], hn[hne_bnd], hn[he_bnd], hn[hse_bnd], hn[hs_bnd], hn[hsw_bnd], hn[hw_bnd], hn[hnw_bnd]
            };
            his_file->put_variable(his_h_name, nst_his, his_values);

            his_values = { qn[p_a], qn[p_b], qn[p_c], qn[p_d], qn[centre], 
                qn[n_bnd], qn[ne_bnd], qn[e_bnd], qn[se_bnd], qn[s_bnd], qn[sw_bnd], qn[w_bnd], qn[nw_bnd],
                qn[hn_bnd], qn[hne_bnd], qn[he_bnd], qn[hse_bnd], qn[hs_bnd], qn[hsw_bnd], qn[hw_bnd], qn[hnw_bnd]
            };
            his_file->put_variable(his_q_name, nst_his, his_values);
 
            his_values = { rn[p_a], rn[p_b], rn[p_c], rn[p_d], rn[centre], 
                rn[n_bnd], rn[ne_bnd], rn[e_bnd], rn[se_bnd], rn[s_bnd], rn[sw_bnd], rn[w_bnd], rn[nw_bnd],
                rn[hn_bnd], rn[hne_bnd], rn[he_bnd], rn[hse_bnd], rn[hs_bnd], rn[hsw_bnd], rn[hw_bnd], rn[hnw_bnd]
            };
            his_file->put_variable(his_r_name, nst_his, his_values);

            his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], 
                s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd],
                s[hn_bnd], s[hne_bnd], s[he_bnd], s[hse_bnd], s[hs_bnd], s[hsw_bnd], s[hw_bnd], s[hnw_bnd]
            };
            his_file->put_variable(his_s_name, nst_his, his_values);

            his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], 
                u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd],
                u[hn_bnd], u[hne_bnd], u[he_bnd], u[hse_bnd], u[hs_bnd], u[hsw_bnd], u[hw_bnd], u[hnw_bnd]
            };
            his_file->put_variable(his_u_name, nst_his, his_values);

            his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], 
                v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd],
                v[hn_bnd], v[hne_bnd], v[he_bnd], v[hse_bnd], v[hs_bnd], v[hsw_bnd], v[hw_bnd], v[hnw_bnd]
            };
            his_file->put_variable(his_v_name, nst_his, his_values);

            his_values = { double(used_newton_iter) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            his_values = { double(solver.iterations()) };
            his_file->put_variable(his_BiCGstab_iter_name, nst_his, his_values);
            his_values = { solver.error() };
            his_file->put_variable(his_BiCGstab_iter_error_name, nst_his, his_values);
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
inline int idx(int i, int j, int ny)
{
    return i * ny + j;
}
inline double c_scv(double c0, double c1, double c2, double c3)
{
    // value at subcontrol volume
    return 0.0625 * (9. * c0 + 3. * c1 +  3. * c2 + c3);
}
inline double dcdx_scv(double c0, double c1, double c2, double c3)
{
    // value quadrature point (i+1/4, j+1/4) at subcontrol volume
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
inline double dcdy_scv(double c0, double c1, double c2, double c3)
{
    // value quadrature point (i+1/4, j+1/4) at subcontrol volume
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
inline double scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdx normal at subcontrol volume edge
    return 0.125 * (3. * c0 + 3. * c1 + c2 + c3);
}
inline double dcdx_scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdx normal at subcontrol volume edge
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
inline double dcdx_scvf_t(double c0, double c1, double c2, double c3)
{
    // dcdx tangential at subcontrol volume edge
    return 0.5 * (c1 - c0 + c2 - c3);
}
inline double dcdy_scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdy normal at subcontrol volume edge
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
inline double dcdy_scvf_t(double c0, double c1, double c2, double c3)
{
    // dcdy tangential at subcontrol volume edge
    return 0.5 * (c3 - c0 + c2 - c1);
}
std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name)
{
    std::string ss_x;
    std::string ss_y;
    ss_x = string_format_with_zeros(x_obs, nsig + 3);
    ss_y = string_format_with_zeros(y_obs, nsig + 3);

    //ss_x << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << x_obs;
    //ss_y << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << y_obs;
    return (obs_name + " (" + ss_x + ", " + ss_y + ") ");
}
std::string string_format_with_zeros(double value, int width) 
{
    std::ostringstream oss;
    if (value < 0) {
        oss << '-';
        oss << std::setw(width - 1) << std::setfill('0') << std::fixed << std::setprecision(1) << -value;
    } else {
        oss << std::setw(width) << std::setfill('0') << std::fixed << std::setprecision(1) << value;
    }
    return oss.str();
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

