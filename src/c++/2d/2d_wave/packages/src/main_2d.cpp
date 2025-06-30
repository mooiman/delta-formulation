//
// Programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
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

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <filesystem>
#include <thread>
#include <toml.h>
#include <string>
#include <algorithm>


// for BiCGstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "bed_level.h"
#include "bed_shear_stress.h"
#include "boundary_condition.h"
#include "cfts.h"
#include "convection.h"
#include "grid.h"
#include "initial_conditions.h"
#include "perf_timer.h"
#include "regularization.h"
#include "ugrid2d.h"

void GetArguments(long argc, char** argv, std::filesystem::path & file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

int p_index(int i, int j, int ny);
double c_scv(double, double, double, double);
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
// Momentum equation: d(q)/dt + 1/2  d(h^2)/dx = 0
// Momentum equation: d(r)/dt + 1/2  d(h^2)/dy = 0

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
        std::cout << "Executable compiled: " << __DATE__ << ", " << __TIME__ << std::endl;
        std::cout << std::endl;
        std::cout << "usage: 2d_wave.exe --toml <input_file>" << std::endl;
        std::cout << "======================================================" << std::endl;
        //std::cin.ignore();
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        exit(1);
    }

    std::cout << "======================================================" << std::endl;
    std::cout << "Executable compiled: " << __DATE__ << ", " << __TIME__ << std::endl;
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
    log_file << "======================================================" << std::endl;
    log_file << "Executable compiled: " << __DATE__ << ", " << __TIME__ << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;
    log_file << "=== Copy of the input file ============================" << std::endl;
    log_file << tbl << "\n";  // Write input TOML file to log_file
    log_file << "=======================================================" << std::endl;
    log_file << std::endl;

    double dt = tbl["Numerics"]["dt"].value_or(double(0.0));  // dt == 0 => default stationary
    if (dt == 0.0) { stationary = true;  }
    // Time
    double tstart = tbl["Time"]["tstart"].value_or(double(0.0));
    double tstop = tbl["Time"]["tstop"].value_or(double(1800.));
    
    // Initial
    std::vector<std::string> ini_vars;  // get the element as an array
    auto vars = tbl["Initial"];
    status = get_toml_array(*vars.as_table(), "ini_vars", ini_vars);
    double gauss_amp = tbl["Initial"]["gauss_amp"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary
    double gauss_mu = tbl["Initial"]["gauss_mu"].value_or(double(-INFINITY));
    double gauss_mu_x = tbl["Initial"]["gauss_mu_x"].value_or(double(0.0));
    double gauss_mu_y = tbl["Initial"]["gauss_mu_y"].value_or(double(0.0));
    if (gauss_mu != -INFINITY)
    {
        gauss_mu_x = gauss_mu; 
        gauss_mu_y = 0.0;  // only shift on x-axis
    }
    double gauss_sigma = tbl["Initial"]["gauss_sigma"].value_or(double(-INFINITY));
    double gauss_sigma_x = tbl["Initial"]["gauss_sigma"].value_or(double(1.0));
    double gauss_sigma_y = tbl["Initial"]["gauss_sigma"].value_or(double(1.0));
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
    if (do_q_equation && do_r_equation)
    {
        model_title = "Linear wave equation, BiCGstab";
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
    bool do_q_convection = tbl_chp["do_q_convection"].value_or(bool(false));  // default, no convection
    bool do_r_convection = tbl_chp["do_r_convection"].value_or(bool(false));  // default, no convection
    
    bool do_viscosity = tbl_chp["do_viscosity"].value_or(bool(false));  // default, no viscosity
    double visc_const = tbl_chp["viscosity"].value_or(double(0.0001));  // default 1e-4

    bool do_bed_shear_stress = tbl_chp["do_bed_shear_stress"].value_or(bool(false));  // default, no bed shear stress
    double chezy_coefficient = tbl_chp["chezy_coefficient"].value_or(double(50.0));

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
    double eps_fabs = tbl_chp["eps_fabs"].value_or(double(1.0e-2));  // epsilon needed to approximate the abs-function by a continues function
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
        log_file << "Error: Failed to open file: " << full_grid_filename.string() << std::endl;
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
    CONVECTION* convection = new CONVECTION(theta, dx, dy, nx, ny);
    double cf = g / (chezy_coefficient * chezy_coefficient);
    BEDSHEARSTRESS* bed_stress = new BEDSHEARSTRESS(theta, dx, dy, cf, eps_fabs, nx, ny);

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
    log_file << "do_q_convection = " << do_q_convection << std::endl;
    log_file << "do_r_convection = " << do_r_convection << std::endl;
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
    log_file << "eps_abs_function = " << eps_fabs << std::endl;
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
    std::vector<double> s_giv(nxny, 0.);                     // water level, needed for post-processing
    std::vector<double> u_giv(nxny, 0.);                     // u-velocity, needed for post-processing
    std::vector<double> v_giv(nxny, 0.);                     // v-velocity, needed for post-processing
    std::vector<double> s(nxny, 0.);                     // water level, needed for post-processing
    std::vector<double> u(nxny, 0.);                     // u-velocity, needed for post-processing
    std::vector<double> v(nxny, 0.);                     // v-velocity, needed for post-processing
    std::vector<double> delta_h(nxny, 0.);
    std::vector<double> delta_q(nxny, 0.);
    std::vector<double> delta_r(nxny, 0.);
    std::vector<double> psi(nxny, 0.);

    Eigen::SparseMatrix<double> A(3 * nxny, 3 * nxny);
    Eigen::VectorXd solution(3 * nxny);                     // solution vector [h, q, r]^{n}
    Eigen::VectorXd rhs(3 * nxny);                          // RHS vector [h, q, r]^{n}
    for (int i = 0; i < 3 * nxny; ++i)
    {
        A.coeffRef(i, i) = 1.0;  // 8.8888888;
        solution[i] = 0.0;  // Delta h, Delta q and Delta r
        rhs[i] = 0.0; 
    }

    double alpha = 1. / 8.;                               // Linear (spatial) interpolation coefficient

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
    std::cout << "Initialisation" << std::endl;
    status = 0;
    double min_zb = *std::min_element(zb_giv.begin(), zb_giv.end());

    log_file << "=======================================================" << std::endl;
    log_file << "Nodes   : " << nx << "x" << ny << "=" << nxny << std::endl;
    log_file << "Elements: " << (nx - 1) << "x" << (ny - 1) << "=" << (nx - 1) * (ny - 1) << std::endl;
    log_file << "Volumes : " << (nx - 2) << "x" << (ny - 2) << "=" << (nx - 2) * (ny - 2) << std::endl;
    log_file << "CFL     : " << std::sqrt(g * std::abs(min_zb)) * dt * std::sqrt(( 1./(dx*dx) + 1./(dy*dy))) << std::endl;
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
        regularization->given_function(zb, psi, zb_giv, nx, ny, dx, dy, c_psi, log_file);
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
        std::cout << "Failed to open file: \'" << map_filename << "\' probably in use" << std::endl;
        log_file << "Failed to open file: \'" << map_filename << "\' probably in use" << std::endl;
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
            p0 = p_index(i    , j, ny);
            p1 = p_index(i + 1, j, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = p_index(i + 1, j    , ny);
            p1 = p_index(i + 1, j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = p_index(i + 1, j + 1, ny);
            p1 = p_index(i    , j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                mesh2d_edge_nodes.push_back(p0);
                mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = p_index(i    , j + 1, ny);
            p1 = p_index(i    , j    , ny);
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
            p0 = p_index(i    , j    , ny);
            p1 = p_index(i + 1, j    , ny);
            p2 = p_index(i + 1, j + 1, ny);
            p3 = p_index(i    , j + 1, ny);
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
            int p0 = p_index(i    , j    , ny);
            int p1 = p_index(i + 1, j    , ny);
            int p2 = p_index(i + 1, j + 1, ny);
            int p3 = p_index(i    , j + 1, ny);
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
            int p0 = p_index(i    , j    , ny);
            int p1 = p_index(i + 1, j    , ny);
            int p2 = p_index(i + 1, j + 1, ny);
            int p3 = p_index(i    , j + 1, ny);
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
    status = map_file->add_variable(map_h_name, dim_names, "sea_floor_depth_below_sea_surface", "Water depth", "m", "mesh2D", "node");
    status = map_file->add_variable(map_q_name, dim_names, "", "Water flux (x)", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_r_name, dim_names, "", "Water flux (y)", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dh_name, dim_names, "", "Delta h^{n+1,p+1}", "m", "mesh2D", "node");
    status = map_file->add_variable(map_dq_name, dim_names, "", "Delta q^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dr_name, dim_names, "", "Delta r^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_s_name, dim_names, "sea_surface_height_above_geoid", "WaterLevel", "m", "mesh2D", "node");
    status = map_file->add_variable(map_u_name, dim_names, "sea_water_x_velocity", "Velocity (x)", "m s-1", "mesh2D", "node");
    status = map_file->add_variable(map_v_name, dim_names, "sea_water_y_velocity", "Velocity (y)", "m s-1", "mesh2D", "node");
    status = map_file->add_variable(map_zb_name, dim_names, "", "Bed level", "m", "mesh2D", "node");

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
    STOP_TIMER(Writing map-file);

    // End define map file
    ////////////////////////////////////////////////////////////////////////////
    // Create time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations
    int centre = p_index(nx / 2, ny / 2, ny);
    int w_bnd  = p_index(1     , ny / 2, ny);  // at boundary, not at virtual point
    int e_bnd  = p_index(nx - 2, ny / 2, ny);  // at boundary, not at virtual point
    int s_bnd  = p_index(nx / 2, 1     , ny);
    int n_bnd  = p_index(nx / 2, ny - 2, ny);

    int sw_bnd = p_index(1     , 1     , ny);
    int ne_bnd = p_index(nx - 2, ny - 2, ny);
    int nw_bnd = p_index(1     , ny - 2, ny);
    int se_bnd = p_index(nx - 2, 1     , ny);
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
    p_a = p_index(int((y_a / dy) + nx / 2), int((x_a / dx) + ny / 2), ny);
    p_b = p_index(int((y_b / dy) + nx / 2), int((x_b / dx) + ny / 2), ny);
    p_c = p_index(int((y_c / dy) + nx / 2), int((x_c / dx) + ny / 2), ny);
    p_d = p_index(int((y_d / dy) + nx / 2), int((x_d / dx) + ny / 2), ny);

    std::vector<double> x_obs = { x[p_a], x[p_b], x[p_c], x[p_d], x[centre], x[n_bnd], x[ne_bnd], x[e_bnd], x[se_bnd], x[s_bnd], x[sw_bnd], x[w_bnd], x[nw_bnd] };
    std::vector<double> y_obs = { y[p_a], y[p_b], y[p_c], y[p_d], y[centre], y[n_bnd], y[ne_bnd], y[e_bnd], y[se_bnd], y[s_bnd], y[sw_bnd], y[w_bnd], y[nw_bnd] };
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
    obs_stations.push_back(setup_obs_name(x[n_bnd] , y[n_bnd] , nsig, "N station"));
    obs_stations.push_back(setup_obs_name(x[ne_bnd], y[ne_bnd], nsig, "NE station"));
    obs_stations.push_back(setup_obs_name(x[e_bnd] , y[e_bnd] , nsig, "E station"));
    obs_stations.push_back(setup_obs_name(x[se_bnd], y[se_bnd], nsig, "SE station"));
    obs_stations.push_back(setup_obs_name(x[s_bnd] , y[s_bnd] , nsig, "S station"));
    obs_stations.push_back(setup_obs_name(x[sw_bnd], y[sw_bnd], nsig, "SW station"));
    obs_stations.push_back(setup_obs_name(x[w_bnd] , y[w_bnd] , nsig, "W station"));
    obs_stations.push_back(setup_obs_name(x[nw_bnd], y[nw_bnd], nsig, "NW station"));

    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_h_name("hn_2d");
    std::string his_q_name("qn_2d");
    std::string his_r_name("rn_2d");
    std::string his_s_name("Water_level");
    std::string his_u_name("U-velocity");
    std::string his_v_name("V-velocity");

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
    std::vector<double> his_values = { hn[p_a], hn[p_b], hn[p_c], hn[p_d], hn[centre], hn[n_bnd], hn[ne_bnd], hn[e_bnd], hn[se_bnd], hn[s_bnd], hn[sw_bnd], hn[w_bnd], hn[nw_bnd] };
    his_file->put_variable(his_h_name, nst_his, his_values);

    his_values = { qn[p_a], qn[p_b], qn[p_c], qn[p_d], qn[centre], qn[n_bnd], qn[ne_bnd], qn[e_bnd], qn[se_bnd], qn[s_bnd], qn[sw_bnd], qn[w_bnd], qn[nw_bnd] };
    his_file->put_variable(his_q_name, nst_his, his_values);
 
    his_values = { rn[p_a], rn[p_b], rn[p_c], rn[p_d], rn[centre], rn[n_bnd], rn[ne_bnd], rn[e_bnd], rn[se_bnd], rn[s_bnd], rn[sw_bnd], rn[w_bnd], rn[nw_bnd] };
    his_file->put_variable(his_r_name, nst_his, his_values);

    his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd] };
    his_file->put_variable(his_s_name, nst_his, his_values);

    his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd] };
    his_file->put_variable(his_u_name, nst_his, his_values);

    his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd] };
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
    std::cout << std::fixed << std::setprecision(2) << "tstart= " << tstart + time << ";   tstop= " << tstart + tstop << ";   dt= " << dt << std::endl;

    // Start time loop
    START_TIMER(Time loop);
    double dh_max = 0.0;
    double dq_max = 0.0;
    double dr_max = 0.0;
    int dh_maxi = 0;
    int dq_maxi = 0;
    int dr_maxi = 0;
    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        time = dt * double(nst);

        double hn_im2, hn_im1, hn_0, hn_ip1, hn_ip2, hn_jm2, hn_jm1, hn_jp1, hn_jp2;  // previous timestep (n)
        double hn_im12, hn_ip12, hn_jm12, hn_jp12;
        double hp_im2, hp_im1, hp_0, hp_ip1, hp_ip2, hp_jm2, hp_jm1, hp_jp1, hp_jp2;  // previous iteration (p)
        double hp_im12, hp_ip12, hp_jm12, hp_jp12;

        double qn_im2, qn_im1, qn_0, qn_ip1, qn_ip2, qn_jm2, qn_jm1, qn_jp1, qn_jp2;  // timestep (n)
        double qn_im12, qn_ip12, qn_jm12, qn_jp12;
        double qp_im2, qp_im1, qp_0, qp_ip1, qp_ip2, qp_jm2, qp_jm1, qp_jp1, qp_jp2;  // iteration (p)
        double qp_im12, qp_ip12, qp_jm12, qp_jp12;

        double rn_im2, rn_im1, rn_0, rn_ip1, rn_ip2, rn_jm2, rn_jm1, rn_jp1, rn_jp2;  // timestep (n)
        double rn_im12, rn_ip12, rn_jm12, rn_jp12;
        double rp_im2, rp_im1, rp_0, rp_ip1, rp_ip2, rp_jm2, rp_jm1, rp_jp1, rp_jp2;  // iteration (p)
        double rp_im12, rp_ip12, rp_jm12, rp_jp12;

        double htheta_0;
        double htheta_im1, htheta_im2, htheta_ip1, htheta_ip2;
        double htheta_im12, htheta_ip12;
        double htheta_jm1, htheta_jm2, htheta_jp1, htheta_jp2;
        double htheta_jm12, htheta_jp12;
        double qtheta_0;
        double qtheta_im1, qtheta_im2, qtheta_ip1, qtheta_ip2, qtheta_jm2, qtheta_jm1, qtheta_jp1, qtheta_jp2;
        double qtheta_im12, qtheta_ip12, qtheta_jm12, qtheta_jp12;
        double rtheta_0;
        double rtheta_ip1, rtheta_ip2, rtheta_im1, rtheta_im2, rtheta_jm2, rtheta_jm1, rtheta_jp1, rtheta_jp2;
        double rtheta_im12, rtheta_ip12, rtheta_jm12, rtheta_jp12;

        double bc0;
        double bc1;
        std::vector<double> bc(4, 0.0);

        int select = 1;  // Constant boundary condition
        boundary_condition(bc[BC_NORTH], bc_vals[BC_NORTH], time, treg, select);
        boundary_condition(bc[BC_EAST ], bc_vals[BC_EAST ], time, treg, select);
        boundary_condition(bc[BC_SOUTH], bc_vals[BC_SOUTH], time, treg, select);
        boundary_condition(bc[BC_WEST ], bc_vals[BC_WEST ], time, treg, select);

        int used_newton_iter = 0;
        int used_lin_solv_iter = 0;
        START_TIMER(Newton iteration);
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
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
            for (int i = 1; i < nx - 1; i++)
            {
                for (int j = 1; j < ny - 1; j++)
                {
                    if (i == nx / 2 + 1 && j == ny / 2 + 1)
                    {
                        int c = 1;
                    }
                    //log_file << "Node: " << i << "  " << j << std::endl;
                    int ph_sw = p_index(i - 1, j - 1, ny);  
                    int ph_s  = p_index(i    , j - 1, ny);  
                    int ph_se = p_index(i + 1, j - 1, ny);  
                    int ph_w  = p_index(i - 1, j    , ny);  
                    int ph_0  = p_index(i    , j    , ny); // central point of control volume
                    int ph_e  = p_index(i + 1, j    , ny);  
                    int ph_nw = p_index(i - 1, j + 1, ny);  
                    int ph_n  = p_index(i    , j + 1, ny);  
                    int ph_ne = p_index(i + 1, j + 1, ny);  

                    // ph_0 + 1;  q-momentum equation
                    // ph_0 + 2;  r-momentum equation

                    double zb_0  = zb[ph_0];
                    double zb_n  = zb[ph_n];
                    double zb_ne = zb[ph_ne];
                    double zb_e  = zb[ph_e];
                    double zb_se = zb[ph_se];
                    double zb_s  = zb[ph_s];
                    double zb_sw = zb[ph_sw];
                    double zb_w  = zb[ph_w];
                    double zb_nw = zb[ph_nw];

                    double htheta_0  = theta * hp[ph_0 ] + (1.0 - theta) * hn[ph_0 ];
                    double htheta_n  = theta * hp[ph_n ] + (1.0 - theta) * hn[ph_n ];
                    double htheta_ne = theta * hp[ph_ne] + (1.0 - theta) * hn[ph_ne];
                    double htheta_e  = theta * hp[ph_e ] + (1.0 - theta) * hn[ph_e ];
                    double htheta_se = theta * hp[ph_se] + (1.0 - theta) * hn[ph_se];
                    double htheta_s  = theta * hp[ph_s ] + (1.0 - theta) * hn[ph_s ];
                    double htheta_sw = theta * hp[ph_sw] + (1.0 - theta) * hn[ph_sw];
                    double htheta_w  = theta * hp[ph_w ] + (1.0 - theta) * hn[ph_w ];
                    double htheta_nw = theta * hp[ph_nw] + (1.0 - theta) * hn[ph_nw];

                    double qtheta_0  = theta * qp[ph_0 ] + (1.0 - theta) * qn[ph_0 ];
                    double qtheta_n  = theta * qp[ph_n ] + (1.0 - theta) * qn[ph_n ];
                    double qtheta_ne = theta * qp[ph_ne] + (1.0 - theta) * qn[ph_ne];
                    double qtheta_e  = theta * qp[ph_e ] + (1.0 - theta) * qn[ph_e ];
                    double qtheta_se = theta * qp[ph_se] + (1.0 - theta) * qn[ph_se];
                    double qtheta_s  = theta * qp[ph_s ] + (1.0 - theta) * qn[ph_s ];
                    double qtheta_sw = theta * qp[ph_sw] + (1.0 - theta) * qn[ph_sw];
                    double qtheta_w  = theta * qp[ph_w ] + (1.0 - theta) * qn[ph_w ];
                    double qtheta_nw = theta * qp[ph_nw] + (1.0 - theta) * qn[ph_nw];

                    double rtheta_0  = theta * rp[ph_0 ] + (1.0 - theta) * rn[ph_0 ];
                    double rtheta_n  = theta * rp[ph_n ] + (1.0 - theta) * rn[ph_n ];
                    double rtheta_ne = theta * rp[ph_ne] + (1.0 - theta) * rn[ph_ne];
                    double rtheta_e  = theta * rp[ph_e ] + (1.0 - theta) * rn[ph_e ];
                    double rtheta_se = theta * rp[ph_se] + (1.0 - theta) * rn[ph_se];
                    double rtheta_s  = theta * rp[ph_s ] + (1.0 - theta) * rn[ph_s ];
                    double rtheta_sw = theta * rp[ph_sw] + (1.0 - theta) * rn[ph_sw];
                    double rtheta_w  = theta * rp[ph_w ] + (1.0 - theta) * rn[ph_w ];
                    double rtheta_nw = theta * rp[ph_nw] + (1.0 - theta) * rn[ph_nw];
                    //
                    // continuity equation (dh/dt ... = 0)
                    //
                    if (do_continuity)
                    {
                        int c_eq = 3 * ph_0;
                        if (i == nx / 2 && j == ny / 2)
                        {
                            int c = 1;
                        }
                        //
                        // time integration, continuity equation
                        //
                        A.coeffRef(c_eq, 3 * ph_sw) = dtinv * dxdy * mass[0] * mass[0];
                        A.coeffRef(c_eq, 3 * ph_s ) = dtinv * dxdy * mass[0] * mass[1];
                        A.coeffRef(c_eq, 3 * ph_se) = dtinv * dxdy * mass[0] * mass[2];
                        //
                        A.coeffRef(c_eq, 3 * ph_w) = dtinv * dxdy * mass[1] * mass[0];
                        A.coeffRef(c_eq, 3 * ph_0) = dtinv * dxdy * mass[1] * mass[1];
                        A.coeffRef(c_eq, 3 * ph_e) = dtinv * dxdy * mass[1] * mass[2];
                        //
                        A.coeffRef(c_eq, 3 * ph_nw) = dtinv * dxdy * mass[2] * mass[0];
                        A.coeffRef(c_eq, 3 * ph_n ) = dtinv * dxdy * mass[2] * mass[1];
                        A.coeffRef(c_eq, 3 * ph_ne) = dtinv * dxdy * mass[2] * mass[2];
                        //
                        rhs[c_eq] =
                            -dtinv * dxdy * mass[0] * mass[0] * (hp[ph_sw] - hn[ph_sw]) +
                            -dtinv * dxdy * mass[0] * mass[1] * (hp[ph_s ] - hn[ph_s ]) +
                            -dtinv * dxdy * mass[0] * mass[2] * (hp[ph_se] - hn[ph_se]) +
                            //
                            -dtinv * dxdy * mass[1] * mass[0] * (hp[ph_w] - hn[ph_w]) +
                            -dtinv * dxdy * mass[1] * mass[1] * (hp[ph_0] - hn[ph_0]) +
                            -dtinv * dxdy * mass[1] * mass[2] * (hp[ph_e] - hn[ph_e]) +
                            //
                            -dtinv * dxdy * mass[2] * mass[0] * (hp[ph_nw] - hn[ph_nw]) +
                            -dtinv * dxdy * mass[2] * mass[1] * (hp[ph_n ] - hn[ph_n ]) +
                            -dtinv * dxdy * mass[2] * mass[2] * (hp[ph_ne] - hn[ph_ne]);
                        //
                        // mass flux, continuity equation
                        //
                        if (false)
                        {
                            if (do_q_equation)
                            {
                                A.coeffRef(c_eq, 3 * ph_e + 1) = 0.5 * theta * dy;
                                A.coeffRef(c_eq, 3 * ph_0 + 1) = 0.0;  //0.5 * theta * dy - 0.5 * theta * dy;
                                A.coeffRef(c_eq, 3 * ph_w + 1) = -0.5 * theta * dy;
                            }
                            //
                            if (do_r_equation)
                            {

                                A.coeffRef(c_eq, 3 * ph_n + 2) = 0.5 * theta * dx;
                                A.coeffRef(c_eq, 3 * ph_0 + 2) = 0.0;  // 0.5 * theta * dx - 0.5 * theta * dx;
                                A.coeffRef(c_eq, 3 * ph_s + 2) = -0.5 * theta * dx;
                            }
                        }
                        else
                        {
                            A.coeffRef(c_eq, 3 * ph_0 + 1) = 0.0;  //
                            if (do_q_equation)
                            {
                                A.coeffRef(c_eq, 3 * ph_sw + 1) = theta * 0.125 * 0.5 * (-dy);  // flux_7
                                A.coeffRef(c_eq, 3 * ph_s  + 1) = theta * 0.125 * 0.5 * (dy - dy); // flux_2 and flux_7
                                A.coeffRef(c_eq, 3 * ph_se + 1) = theta * 0.125 * 0.5 * (dy);  // flux_2
                                A.coeffRef(c_eq, 3 * ph_e  + 1) = theta * 0.125 * 0.5 * (dy) * (3. + 3.);  // flux_2 and flux_3
                                A.coeffRef(c_eq, 3 * ph_ne + 1) = theta * 0.125 * 0.5 * (dy);  // flux_3
                                A.coeffRef(c_eq, 3 * ph_n  + 1) = theta * 0.125 * 0.5 * (dy - dy);  // flux_3 and flux_6
                                A.coeffRef(c_eq, 3 * ph_nw + 1) = theta * 0.125 * 0.5 * (-dy);  // flux_6
                                A.coeffRef(c_eq, 3 * ph_w  + 1) = theta * 0.125 * 0.5 * (-dy) * (3. + 3.);  // flux_6 and flux_7
                            }
                            //
                            if (do_r_equation)
                            {
                                A.coeffRef(c_eq, 3 * ph_sw + 2) = theta * 0.125 * 0.5 * (-dx);  // flux_0
                                A.coeffRef(c_eq, 3 * ph_s  + 2) = theta * 0.125 * 0.5 * (-dx) * (3. + 3.);  // flux_0 and flux_1
                                A.coeffRef(c_eq, 3 * ph_se + 2) = theta * 0.125 * 0.5 * (-dx);  // flux 1
                                A.coeffRef(c_eq, 3 * ph_e  + 2) = theta * 0.125 * 0.5 * (-dx + dx);  // flux_1 and flux_4
                                A.coeffRef(c_eq, 3 * ph_ne + 2) = theta * 0.125 * 0.5 * (dx);  // flux_4
                                A.coeffRef(c_eq, 3 * ph_n  + 2) = theta * 0.125 * 0.5 * (dx) * (3. + 3.);  // flux_4 and flux_5
                                A.coeffRef(c_eq, 3 * ph_nw + 2) = theta * 0.125 * 0.5 * (dx);  // flux_5
                                A.coeffRef(c_eq, 3 * ph_w  + 2) = theta * 0.125 * 0.5 * (-dx + dx);  // flux_0 and flux_5
                            }
                        }
                        //
                        // RHS continuity equation
                        //
                        // eight fluxes
                        //

                        double fluxy_0 = 0.25 * dx * dcdy_scvf_n(rtheta_0, rtheta_s, rtheta_w, rtheta_sw);
                        double fluxy_1 = 0.25 * dx * dcdy_scvf_n(rtheta_0, rtheta_s, rtheta_e, rtheta_se);
                        double fluxx_2 = 0.25 * dy * dcdx_scvf_n(qtheta_e, qtheta_0, qtheta_se, qtheta_s);
                        double fluxx_3 = 0.25 * dy * dcdx_scvf_n(qtheta_e, qtheta_0, qtheta_ne, qtheta_n);
                        double fluxy_4 = 0.25 * dx * dcdy_scvf_n(rtheta_n, rtheta_0, rtheta_ne, rtheta_e);
                        double fluxy_5 = 0.25 * dx * dcdy_scvf_n(rtheta_n, rtheta_0, rtheta_nw, rtheta_w);
                        double fluxx_6 = 0.25 * dy * dcdx_scvf_n(qtheta_0, qtheta_w, qtheta_n, qtheta_nw);
                        double fluxx_7 = 0.25 * dy * dcdx_scvf_n(qtheta_0, qtheta_w, qtheta_s, qtheta_sw);
                        //
                        double a = (fluxy_0 + fluxy_1 + fluxx_2 + fluxx_3 + fluxy_4 + fluxy_5 + fluxx_6 + fluxx_7);
                        rhs[c_eq] += -a;
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
                        //log_file << std::scientific << std::setprecision(8) << "p(" << i << "," << j << ") " << "c_eq. " << c_eq << "  RHS: " << dqdx << "  b:  " << b << std::endl;
                    }
                    //
                    // q-momentum(dq/dt ... = 0) 
                    //
                    if (do_q_equation)
                    {
                        int q_eq = 3 * ph_0 + 1;
                        if (i == nx / 2 && j == ny / 2)
                        {
                            int c = 1;
                        }
                        //
                        // time integration, q-momentum equation
                        //
                        A.coeffRef(q_eq, 3 * ph_sw + 1) = dtinv * dxdy * mass[0] * mass[0];
                        A.coeffRef(q_eq, 3 * ph_s  + 1) = dtinv * dxdy * mass[0] * mass[1];
                        A.coeffRef(q_eq, 3 * ph_se + 1) = dtinv * dxdy * mass[0] * mass[2];
                        //              
                        A.coeffRef(q_eq, 3 * ph_w + 1) = dtinv * dxdy * mass[1] * mass[0];
                        A.coeffRef(q_eq, 3 * ph_0 + 1) = dtinv * dxdy * mass[1] * mass[1];
                        A.coeffRef(q_eq, 3 * ph_e + 1) = dtinv * dxdy * mass[1] * mass[2];
                        //            
                        A.coeffRef(q_eq, 3 * ph_nw + 1) = dtinv * dxdy * mass[2] * mass[0];
                        A.coeffRef(q_eq, 3 * ph_n  + 1) = dtinv * dxdy * mass[2] * mass[1];
                        A.coeffRef(q_eq, 3 * ph_ne + 1) = dtinv * dxdy * mass[2] * mass[2];
                        //
                        rhs[q_eq] =
                            -dtinv * dxdy * mass[0] * mass[0] * (qp[ph_sw] - qn[ph_sw]) +
                            -dtinv * dxdy * mass[0] * mass[1] * (qp[ph_s ] - qn[ph_s ]) +
                            -dtinv * dxdy * mass[0] * mass[2] * (qp[ph_se] - qn[ph_se]) +
                            //
                            -dtinv * dxdy * mass[1] * mass[0] * (qp[ph_w] - qn[ph_w]) +
                            -dtinv * dxdy * mass[1] * mass[1] * (qp[ph_0] - qn[ph_0]) +
                            -dtinv * dxdy * mass[1] * mass[2] * (qp[ph_e] - qn[ph_e]) +
                            //
                            -dtinv * dxdy * mass[2] * mass[0] * (qp[ph_nw] - qn[ph_nw]) +
                            -dtinv * dxdy * mass[2] * mass[1] * (qp[ph_n ] - qn[ph_n ]) +
                            -dtinv * dxdy * mass[2] * mass[2] * (qp[ph_ne] - qn[ph_ne]);
                        //
                        double depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
                        double depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
                        double depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
                        double depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

                        double scv_area = 0.25 * dxdy;
                        double dzetadx_0 = scv_area / dx * dcdx_scv(htheta_0 + zb_0, htheta_w + zb_w, htheta_s + zb_s , htheta_sw + zb_sw );
                        double dzetadx_1 = scv_area / dx * dcdx_scv(htheta_e + zb_e, htheta_0 + zb_0, htheta_se + zb_se, htheta_s  + zb_s );
                        double dzetadx_2 = scv_area / dx * dcdx_scv(htheta_e + zb_e, htheta_0 + zb_0, htheta_ne + zb_ne, htheta_n  + zb_n );
                        double dzetadx_3 = scv_area / dx * dcdx_scv(htheta_0 + zb_0, htheta_w + zb_w, htheta_n  + zb_n , htheta_nw + zb_nw);
                        //
                        if (false)
                        {
                            A.coeffRef(q_eq, 3 * ph_e) =  0.5 * theta * g * dy * htheta_ip12;
                            A.coeffRef(q_eq, 3 * ph_0) = -0.5 * theta * g * dy * htheta_im12 + 0.5 * theta * g * dy * htheta_ip12;
                            A.coeffRef(q_eq, 3 * ph_w) = -0.5 * theta * g * dy * htheta_im12;
                        }
                        else
                        {   // theta * dzeta/dx * Delta h
                            double fac = theta * scv_area * g * 0.0625;
                            A.coeffRef(q_eq, 3 * ph_0 ) = fac * (9. * dzetadx_0 + 9. * dzetadx_1 + 9. * dzetadx_2 + 9. * dzetadx_3) ;   
                            A.coeffRef(q_eq, 3 * ph_sw) = fac * (1. * dzetadx_0);
                            A.coeffRef(q_eq, 3 * ph_s ) = fac * (3. * dzetadx_1 + 3. * dzetadx_0); 
                            A.coeffRef(q_eq, 3 * ph_se) = fac * (1. * dzetadx_1);
                            A.coeffRef(q_eq, 3 * ph_e ) = fac * (3. * dzetadx_1 + 3. * dzetadx_2);
                            A.coeffRef(q_eq, 3 * ph_ne) = fac * (1. * dzetadx_2);
                            A.coeffRef(q_eq, 3 * ph_n ) = fac * (3. * dzetadx_2 + 3. * dzetadx_3);   
                            A.coeffRef(q_eq, 3 * ph_nw) = fac * (1. * dzetadx_3);
                            A.coeffRef(q_eq, 3 * ph_w)  = fac * (3. * dzetadx_0 + 3. * dzetadx_3);
                            // theta * h * d(Delta zeta)/dx 
                            fac = theta * scv_area * g * 0.25/dx;
                            A.coeffRef(q_eq, 3 * ph_0 ) += fac * (+3. * depth_0 - 3. * depth_1 - 3. * depth_2 + 3. * depth_3);
                            A.coeffRef(q_eq, 3 * ph_sw) += fac * (-1. * depth_0);
                            A.coeffRef(q_eq, 3 * ph_s ) += fac * (+1. * depth_0 - 1. * depth_1);
                            A.coeffRef(q_eq, 3 * ph_se) += fac * (+1. * depth_1 );
                            A.coeffRef(q_eq, 3 * ph_e ) += fac * (+3. * depth_1 + 3. * depth_3);
                            A.coeffRef(q_eq, 3 * ph_ne) += fac * (+1. * depth_2);
                            A.coeffRef(q_eq, 3 * ph_n ) += fac * (-1. * depth_2 + 1. * depth_3);
                            A.coeffRef(q_eq, 3 * ph_nw) += fac * (-1. * depth_3);
                            A.coeffRef(q_eq, 3 * ph_w ) += fac * (-3. * depth_0 - 3. * depth_3);

                        }
                        //
                        // RHS q-momentum equation
                        //
                        double rhs_q = g * (depth_0 * dzetadx_0 + depth_1 * dzetadx_1 + depth_2 * dzetadx_2 + depth_3 * dzetadx_3);
                        rhs[q_eq] += -rhs_q;
                    }
                    // 
                    // r-momentum (dr/dt ... = 0)
                    //
                    if (do_r_equation)
                    {
                        int r_eq = 3 * ph_0 + 2;
                        //
                        // time integration, r-momentum equation
                        //
                        A.coeffRef(r_eq, 3 * ph_sw + 2) = dtinv * dxdy * mass[0] * mass[0];
                        A.coeffRef(r_eq, 3 * ph_s  + 2) = dtinv * dxdy * mass[0] * mass[1];
                        A.coeffRef(r_eq, 3 * ph_se + 2) = dtinv * dxdy * mass[0] * mass[2];
                        //
                        A.coeffRef(r_eq, 3 * ph_w + 2) = dtinv * dxdy * mass[1] * mass[0];
                        A.coeffRef(r_eq, 3 * ph_0 + 2) = dtinv * dxdy * mass[1] * mass[1];
                        A.coeffRef(r_eq, 3 * ph_e + 2) = dtinv * dxdy * mass[1] * mass[2];
                        //
                        A.coeffRef(r_eq, 3 * ph_nw + 2) = dtinv * dxdy * mass[2] * mass[0];
                        A.coeffRef(r_eq, 3 * ph_n  + 2) = dtinv * dxdy * mass[2] * mass[1];
                        A.coeffRef(r_eq, 3 * ph_ne + 2) = dtinv * dxdy * mass[2] * mass[2];
                        //
                        rhs[r_eq] =
                            -dtinv * dxdy * mass[0] * mass[0] * (rp[ph_sw] - rn[ph_sw]) +
                            -dtinv * dxdy * mass[0] * mass[1] * (rp[ph_s ] - rn[ph_s ]) +
                            -dtinv * dxdy * mass[0] * mass[2] * (rp[ph_se] - rn[ph_se]) +
                            //
                            -dtinv * dxdy * mass[1] * mass[0] * (rp[ph_w] - rn[ph_w]) +
                            -dtinv * dxdy * mass[1] * mass[1] * (rp[ph_0] - rn[ph_0]) +
                            -dtinv * dxdy * mass[1] * mass[2] * (rp[ph_e] - rn[ph_e]) +
                            //
                            -dtinv * dxdy * mass[2] * mass[0] * (rp[ph_nw] - rn[ph_nw]) +
                            -dtinv * dxdy * mass[2] * mass[1] * (rp[ph_n ] - rn[ph_n ]) +
                            -dtinv * dxdy * mass[2] * mass[2] * (rp[ph_ne] - rn[ph_ne]);
                        //
                        double scv_area = 0.25 * dxdy;

                        double depth_0 = c_scv(htheta_0, htheta_w, htheta_s, htheta_sw);
                        double depth_1 = c_scv(htheta_0, htheta_s, htheta_e, htheta_se);
                        double depth_2 = c_scv(htheta_0, htheta_e, htheta_n, htheta_ne);
                        double depth_3 = c_scv(htheta_0, htheta_n, htheta_w, htheta_nw);

                        double dzetady_0 = scv_area / dy * dcdy_scv(htheta_0 + zb_0, htheta_s + zb_s, htheta_w  + zb_w , htheta_sw + zb_sw);
                        double dzetady_1 = scv_area / dy * dcdy_scv(htheta_0 + zb_0, htheta_s + zb_s, htheta_e  + zb_e , htheta_se + zb_se);
                        double dzetady_2 = scv_area / dy * dcdy_scv(htheta_n + zb_n, htheta_0 + zb_0, htheta_ne + zb_ne, htheta_e  + zb_e );
                        double dzetady_3 = scv_area / dy * dcdy_scv(htheta_n + zb_n, htheta_0 + zb_0, htheta_nw + zb_nw, htheta_w  + zb_w );

                        if (false)
                        {
                            A.coeffRef(r_eq, 3 * ph_n) =  0.5 * theta * g * dx * htheta_jp12;
                            A.coeffRef(r_eq, 3 * ph_0) = -0.5 * theta * g * dx * htheta_jm12 + 0.5 * theta * g * dx * htheta_jp12;
                            A.coeffRef(r_eq, 3 * ph_s) = -0.5 * theta * g * dx * htheta_jm12;
                        }
                        else
                        {   // theta * dzeta/dy * Delta h
                            double fac = theta * scv_area * g * 0.0625;
                            A.coeffRef(r_eq, 3 * ph_0 ) = (9. * dzetady_0 + 9. * dzetady_1 + 9. * dzetady_2 + 9. * dzetady_3);
                            A.coeffRef(r_eq, 3 * ph_sw) = fac * (1. * dzetady_0);
                            A.coeffRef(r_eq, 3 * ph_s ) = fac * (3. * dzetady_1 + 3. * dzetady_0);
                            A.coeffRef(r_eq, 3 * ph_se) = fac * (1. * dzetady_1);
                            A.coeffRef(r_eq, 3 * ph_e ) = fac * (3. * dzetady_1 + 3. * dzetady_2);
                            A.coeffRef(r_eq, 3 * ph_ne) = fac * (1. * dzetady_2);
                            A.coeffRef(r_eq, 3 * ph_n ) = fac * (3. * dzetady_2 + 3. * dzetady_3);
                            A.coeffRef(r_eq, 3 * ph_nw) = fac * (1. * dzetady_3);
                            A.coeffRef(r_eq, 3 * ph_w ) = fac * (3. * dzetady_0 + 3. * dzetady_3);
                            // theta * h * d(Delta zeta)/dy
                            fac = theta * scv_area * g * 0.25 / dy;
                            A.coeffRef(r_eq, 3 * ph_0 ) += fac * (+3. * depth_0 + 3. * depth_1 - 3. * depth_2 - 3* depth_3);
                            A.coeffRef(r_eq, 3 * ph_sw) += fac * (-1. * depth_0);
                            A.coeffRef(r_eq, 3 * ph_s ) += fac * (-3. * depth_0 - 3. * depth_1);
                            A.coeffRef(r_eq, 3 * ph_se) += fac * (-1. * depth_1);
                            A.coeffRef(r_eq, 3 * ph_e ) += fac * (+1. * depth_1 - 1. * depth_2);
                            A.coeffRef(r_eq, 3 * ph_ne) += fac * (+1. * depth_2);
                            A.coeffRef(r_eq, 3 * ph_n ) += fac * (+3. * depth_2 + 3. * depth_3);
                            A.coeffRef(r_eq, 3 * ph_nw) += fac * (+1. * depth_3);
                            A.coeffRef(r_eq, 3 * ph_w ) += fac * (+1. * depth_0 - 1. * depth_3);
                        }
                        // 
                        // RHS r-momentum equation
                        //
                        double rhs_r = g * (depth_0 * dzetady_0 + depth_1 * dzetady_1 + depth_2 * dzetady_2 + depth_3 * dzetady_3);
                        rhs[r_eq] += -rhs_r;
                    }
                }  // end interior nodes
            }
            STOP_TIMER(Linear wave);
            if (do_q_convection)
            {
                START_TIMER(Convection q);
                //
                // convection
                //
                STOP_TIMER(Convection q);
            }
            if (do_r_convection)
            {
                START_TIMER(Convection r);
                //
                // convection
                //
                STOP_TIMER(Convection r);
            }
            if (do_bed_shear_stress)
            {
                START_TIMER(Bed shear stress);
                bed_stress->matrix_2d_qr_eq(A, rhs, hp, qp, rp);
                STOP_TIMER(Bed shear stress);
            }
            if (do_viscosity)
            {
                START_TIMER(Viscosity);
                //
                // viscosity
                //
                STOP_TIMER(Viscosity);
            }

            //
            //boundary nodes
            //
            if (true)  // nnorth boundary
            {
                int j = ny - 1;
                for (int i = 1; i < nx - 1; ++i)
                {
                    int ph_0  = p_index(i, j    , ny);  // central point of control volume
                    int ph_s  = p_index(i, j - 1, ny);
                    int ph_ss = p_index(i, j - 2, ny); 

                    int c_eq = 3 * ph_0;
                    int q_eq = c_eq + 1;
                    int r_eq = c_eq + 2;

                    hn_0 = hn[ph_0];
                    hn_jm1 = hn[ph_s];
                    hn_jm2 = hn[ph_ss];
                    hp_0 = hp[ph_0];
                    hp_jm1 = hp[ph_s];
                    hp_jm2 = hp[ph_ss];
                    htheta_0 = theta * hp_0 + (1.0 - theta) * hn_0;
                    htheta_jm1 = theta * hp_jm1 + (1.0 - theta) * hn_jm1;
                    htheta_jm2 = theta * hp_jm2 + (1.0 - theta) * hn_jm2;

                    hn_jm12 = w_nat[0] * hn_0 + w_nat[1] * hn_jm1 + w_nat[2] * hn_jm2;
                    hp_jm12 = w_nat[0] * hp_0 + w_nat[1] * hp_jm1 + w_nat[2] * hp_jm2;
                    htheta_jm12 = w_nat[0] * htheta_0 + w_nat[1] * htheta_jm1 + w_nat[2] * htheta_jm2;

                    qn_0 = qn[ph_0];
                    qn_jm1 = qn[ph_s];
                    qn_jm2 = qn[ph_ss];
                    qp_0 = qp[ph_0];
                    qp_jm1 = qp[ph_s];
                    qp_jm2 = qp[ph_ss];
                    qtheta_0 = theta * qp_0 + (1.0 - theta) * qn_0;
                    qtheta_jm1 = theta * qp_jm1 + (1.0 - theta) * qn_jm1;
                    qtheta_jm2 = theta * qp_jm2 + (1.0 - theta) * qn_jm2;

                    qn_jm12 = w_nat[0] * qn_0 + w_nat[1] * qn_jm1 + w_nat[2] * qn_jm2;
                    qp_jm12 = w_nat[0] * qp_0 + w_nat[1] * qp_jm1 + w_nat[2] * qp_jm2;
                    qtheta_jm12 = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_jm1 + w_nat[2] * qtheta_jm2;

                    rn_0 = rn[ph_0];
                    rn_jm1 = rn[ph_s];
                    rn_jm2 = rn[ph_ss];
                    rp_0 = rp[ph_0];
                    rp_jm1 = rp[ph_s];
                    rp_jm2 = rp[ph_ss];
                    rtheta_0 = theta * rp_0 + (1.0 - theta) * rn_0;
                    rtheta_jm1 = theta * rp_jm1 + (1.0 - theta) * rn_jm1;
                    rtheta_jm2 = theta * rp_jm2 + (1.0 - theta) * rn_jm2;

                    rn_jm12 = w_nat[0] * rn_0 + w_nat[1] * rn_jm1 + w_nat[2] * rn_jm2;
                    rp_jm12 = w_nat[0] * rp_0 + w_nat[1] * rp_jm1 + w_nat[2] * rp_jm2;
                    rtheta_jm12 = w_nat[0] * rtheta_0 + w_nat[1] * rtheta_jm1 + w_nat[2] * rtheta_jm2;

                    double zb_jm12 = w_nat[0] * zb[ph_0] + w_nat[1] * zb[ph_s] + w_nat[2] * zb[ph_ss];
                    double h_infty = -zb_jm12;
                    double c_wave = std::sqrt(g * h_infty);

                    if (bc_type[BC_NORTH] == "dirichlet" || bc_absorbing[BC_NORTH] == false)
                    {
                        //
                        // continuity equation
                        //
                        A.coeffRef(c_eq, 3 * ph_0) = 1.0;
                        A.coeffRef(c_eq, 3 * ph_s) = -1.0;
                        A.coeffRef(c_eq, 3 * ph_ss) = 0.0;
                        rhs[c_eq] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(q_eq, 3 * ph_0 + 1) = 1.0;
                        A.coeffRef(q_eq, 3 * ph_s + 1) = -1.0;
                        A.coeffRef(q_eq, 3 * ph_ss + 1) = 0.0;
                        rhs[q_eq] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(r_eq, 3 * ph_0 + 2) = 0.5;
                        A.coeffRef(r_eq, 3 * ph_s + 2) = 0.5;
                        A.coeffRef(r_eq, 3 * ph_ss + 2) = 0.0;
                        rhs[r_eq] = 0.0;
                    }
                    else
                    {
                        if (bc_absorbing[BC_NORTH])
                        {
                            double con_fac = c_wave;
                            if (bc_type[BC_NORTH] == "mooiman")
                            {
                                //
                                // continuity
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = w_nat[0] * -c_wave;
                                A.coeffRef(c_eq, 3 * ph_s ) = w_nat[1] * -c_wave;
                                A.coeffRef(c_eq, 3 * ph_ss) = w_nat[2] * -c_wave;
                                // flow_x 
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_s  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ss + 1) = 0.0;
                                // flow y
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = w_nat[0];
                                A.coeffRef(c_eq, 3 * ph_s  + 2) = w_nat[1];
                                A.coeffRef(c_eq, 3 * ph_ss + 2) = w_nat[2];
                                //
                                rhs[c_eq] = -(rp_jm12 - c_wave * (hp_jm12 - h_infty));
                            }
                            if (bc_type[BC_NORTH] == "borsboom")
                            {
                                //
                                // Essential boundary condition
                                // ----------------------------
                                //
                                if (do_r_convection) { con_fac = c_wave - rp_jm12 / hp_jm12; }
                                A.coeffRef(c_eq, 3 * ph_0)  = dtinv * -con_fac * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_s)  = dtinv * -con_fac * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ss) = dtinv * -con_fac * w_ess[2];
                                //
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_s  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ss + 1) = 0.0;
                                //
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = dtinv * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_s  + 2) = dtinv * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ss + 2) = dtinv * w_ess[2];
                                //
                                double dhdt = dtinv * (hp_jm12 - hn_jm12);
                                double drdt = dtinv * (rp_jm12 - rn_jm12);
                                rhs[c_eq] = -(drdt - con_fac * dhdt);

                                double corr_term = 0.0;
                                if (bc_vars[BC_NORTH] == "zeta")
                                {
                                    if (stationary) { sign = -1.0; }
                                    corr_term = + sign * eps_bc_corr * (hp_im12 - (bc[BC_NORTH] - zb_jm12));
                                    rhs[c_eq] += corr_term;
                                    sign = 1.0;
                                }
                                if (bc_vars[BC_NORTH] == "q")
                                {
                                    if (stationary) { sign = -1.0; }
                                    corr_term = - sign * eps_bc_corr * (qp_im12 - bc[BC_NORTH]);
                                    rhs[c_eq] += corr_term;
                                    sign = 1.0;
                                }
                            }
                            //
                            // natural boundary condition (third equation)
                            // --------------------------
                            //
                            double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                                + dtinv * (hp_jm1 - hn_jm1) * w_nat[1]
                                + dtinv * (hp_jm2 - hn_jm2) * w_nat[2];
                            double drdy = dyinv * (rtheta_0 - rtheta_jm1);
                            double drdt = dtinv * (rp_0 - rn_0) * w_nat[0]
                                + dtinv * (rp_jm1 - rn_jm1) * w_nat[1]
                                + dtinv * (rp_jm2 - rn_jm2) * w_nat[2];
                            double dzetady = dyinv * (htheta_0 + zb[ph_0] - htheta_jm1 - zb[ph_s]);
                            //
                            // q-momentum (tangential equation, q == 0)
                            //
                            A.coeffRef(q_eq, 3 * ph_0 ) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_s ) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ss) = 0.0;
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 1) = 1.0;
                            A.coeffRef(q_eq, 3 * ph_s  + 1) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ss + 1) = 0.0;
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_s  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ss + 2) = 0.0;
                            //
                            rhs[q_eq] = 0.0;
                            //
                            // third equation (mom_eq + c_wave * con_eq)
                            //
                            A.coeffRef(r_eq, 3 * ph_0 ) = w_nat[0] * theta * g * dzetady + dyinv * theta * g * htheta_jm12;
                            A.coeffRef(r_eq, 3 * ph_s ) = w_nat[1] * theta * g * dzetady - dyinv * theta * g * htheta_jm12;
                            A.coeffRef(r_eq, 3 * ph_ss) = w_nat[2] * theta * g * dzetady;
                            //
                            A.coeffRef(r_eq, 3 * ph_0  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_s  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ss + 1) = 0.0;
                            //
                            A.coeffRef(r_eq, 3 * ph_0  + 2) = dtinv * w_nat[0];
                            A.coeffRef(r_eq, 3 * ph_s  + 2) = dtinv * w_nat[1];
                            A.coeffRef(r_eq, 3 * ph_ss + 2) = dtinv * w_nat[2];
                            //
                            rhs[r_eq] = -(drdt + g * htheta_jm12 * dzetady);
                            //
                            // continuity part (added and multiplied by +c_wave)
                            //
                             if (do_r_convection) { con_fac = c_wave + rp_jm12 / hp_jm12; }
                            // Contribution Delta h
                            A.coeffRef(r_eq, 3 * ph_0 ) += con_fac * dtinv * w_nat[0];
                            A.coeffRef(r_eq, 3 * ph_s ) += con_fac * dtinv * w_nat[1];
                            A.coeffRef(r_eq, 3 * ph_ss) += con_fac * dtinv * w_nat[2];
                            // Contribution Delta q
                            A.coeffRef(r_eq, 3 * ph_0  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_s  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ss + 1) = 0.0;
                            // Contribution Delta r
                            A.coeffRef(r_eq, 3 * ph_0  + 2) += con_fac * dyinv * theta;
                            A.coeffRef(r_eq, 3 * ph_s  + 2) += con_fac * -dyinv * theta;
                            A.coeffRef(r_eq, 3 * ph_ss + 2) += 0.0;
                            //
                            rhs[r_eq] += -con_fac * (dhdt + drdy);
                       }
                    }
                }
            }
            if (true)  // eeast boundary
            {
                int i = nx - 1;
                for (int j = 1; j < ny - 1; ++j)
                {
                    int ph_0  = p_index(i    , j, ny);  // central point of control volume
                    int ph_w  = p_index(i - 1, j, ny);
                    int ph_ww = p_index(i - 2, j, ny);

                    int c_eq = 3 * ph_0;
                    int q_eq = c_eq + 1;
                    int r_eq = c_eq + 2;

                    hn_0 = hn[ph_0];
                    hn_im1 = hn[ph_w];
                    hn_im2 = hn[ph_ww];
                    hp_0 = hp[ph_0];
                    hp_im1 = hp[ph_w];
                    hp_im2 = hp[ph_ww];
                    htheta_0 = theta * hp_0 + (1.0 - theta) * hn_0;
                    htheta_im1 = theta * hp_im1 + (1.0 - theta) * hn_im1;
                    htheta_im2 = theta * hp_im2 + (1.0 - theta) * hn_im2;

                    hn_im12 = w_nat[0] * hn_0 + w_nat[1] * hn_im1 + w_nat[2] * hn_im2;
                    hp_im12 = w_nat[0] * hp_0 + w_nat[1] * hp_im1 + w_nat[2] * hp_im2;
                    htheta_im12 = w_nat[0] * htheta_0 + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;

                    qn_0 = qn[ph_0];
                    qn_im1 = qn[ph_w];
                    qn_im2 = qn[ph_ww];
                    qp_0 = qp[ph_0];
                    qp_im1 = qp[ph_w];
                    qp_im2 = qp[ph_ww];
                    qtheta_0 = theta * qp_0 + (1.0 - theta) * qn_0;
                    qtheta_im1 = theta * qp_im1 + (1.0 - theta) * qn_im1;
                    qtheta_im2 = theta * qp_im2 + (1.0 - theta) * qn_im2;

                    qn_im12 = w_nat[0] * qn_0 + w_nat[1] * qn_im1 + w_nat[2] * qn_im2;
                    qp_im12 = w_nat[0] * qp_0 + w_nat[1] * qp_im1 + w_nat[2] * qp_im2;
                    qtheta_im12 = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;

                    rn_0 = rn[ph_0];
                    rn_im1 = rn[ph_w];
                    rn_im2 = rn[ph_ww];
                    rp_0 = rp[ph_0];
                    rp_im1 = rp[ph_w];
                    rp_im2 = rp[ph_ww];
                    rtheta_0 = theta * rp_0 + (1.0 - theta) * rn_0;
                    rtheta_im1 = theta * rp_im1 + (1.0 - theta) * rn_im1;
                    rtheta_im2 = theta * rp_im2 + (1.0 - theta) * rn_im2;

                    rn_im12 = w_nat[0] * rn_0 + w_nat[1] * rn_im1 + w_nat[2] * rn_im2;
                    rp_im12 = w_nat[0] * rp_0 + w_nat[1] * rp_im1 + w_nat[2] * rp_im2;
                    rtheta_im12 = w_nat[0] * rtheta_0 + w_nat[1] * rtheta_im1 + w_nat[2] * rtheta_im2;

                    double zb_im12 = w_nat[0] * zb[ph_0] + w_nat[1] * zb[ph_w] + w_nat[2] * zb[ph_ww];
                    double h_infty = -zb_im12;
                    double c_wave = std::sqrt(g * h_infty);

                    if (bc_type[BC_EAST] == "dirichlet" || bc_absorbing[BC_EAST] == false)
                    {
                        //
                        // continuity equation
                        //
                        A.coeffRef(c_eq, 3 * ph_0) = 1.0;
                        A.coeffRef(c_eq, 3 * ph_w) = -1.0;
                        A.coeffRef(c_eq, 3 * ph_ww ) = 0.0;
                        rhs[c_eq] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(q_eq, 3 * ph_0 + 1) = 0.5;
                        A.coeffRef(q_eq, 3 * ph_w + 1) = 0.5;
                        A.coeffRef(q_eq, 3 * ph_ww + 1) = 0.0;
                        rhs[q_eq] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(r_eq, 3 * ph_0 + 2) = 1.0;
                        A.coeffRef(r_eq, 3 * ph_w + 2) = -1.0;
                        A.coeffRef(r_eq, 3 * ph_ww + 2) = 0.0;
                        rhs[r_eq] = 0.0;
                    }
                    else
                    {
                        if (bc_absorbing[BC_EAST])
                        {
                            double con_fac = c_wave;
                            if (bc_type[BC_EAST] == "mooiman")
                            {
                                //
                                // continuity
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = w_nat[0] * -c_wave;
                                A.coeffRef(c_eq, 3 * ph_w ) = w_nat[1] * -c_wave;
                                A.coeffRef(c_eq, 3 * ph_ww) = w_nat[2] * -c_wave;
                                // flow x
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = w_nat[0];
                                A.coeffRef(c_eq, 3 * ph_w  + 1) = w_nat[1];
                                A.coeffRef(c_eq, 3 * ph_ww + 1) = w_nat[2];
                                // flow y
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_w  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ww + 2) = 0.0;
                                //
                                rhs[c_eq] = -(qp_im12 - c_wave * (hp_im12 - h_infty));
                            }
                            if (bc_type[BC_EAST] == "borsboom")
                            {
                                //
                                // Essential boundary condition
                                // ----------------------------
                                if (do_q_convection) { con_fac = c_wave - qp_im12 / hp_im12; }
                                //
                                A.coeffRef(c_eq, 3 * ph_0) = dtinv * -con_fac * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_w) = dtinv * -con_fac * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ww) = dtinv * -con_fac * w_ess[2];
                                // flow x
                                A.coeffRef(c_eq, 3 * ph_0 + 1) = dtinv * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_w + 1) = dtinv * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ww + 1) = dtinv * w_ess[2];
                                // flow y
                                A.coeffRef(c_eq, 3 * ph_0 + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_w + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ww + 2) = 0.0;
                                //
                                double dhdt = dtinv * (hp_im12 - hn_im12);
                                double dqdt = dtinv * (qp_im12 - qn_im12);
                                rhs[c_eq] = -(dqdt - con_fac * dhdt);

                                double corr_term = 0.0;
                                if (bc_vars[BC_EAST] == "zeta")
                                {
                                    if (stationary) { sign = -1.0; }
                                    corr_term = + sign * eps_bc_corr * (hp_im12 - (bc[BC_EAST] - zb_im12));
                                    rhs[c_eq] += corr_term;
                                    sign = 1.0;
                                }
                                if (bc_vars[BC_EAST] == "q")
                                {
                                    if (stationary) { sign = -1.0; }
                                    corr_term = - sign * eps_bc_corr * (qp_im12 - bc[BC_EAST]);
                                    rhs[c_eq] += corr_term;
                                    sign = 1.0;
                                }
                            }
                            //
                            // natural boundary condition (second equation)
                            // --------------------------
                            //
                            double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                                + dtinv * (hp_im1 - hn_im1) * w_nat[1]
                                + dtinv * (hp_im2 - hn_im2) * w_nat[2];
                            double dqdx = dxinv * (qtheta_0 - qtheta_im1);
                            double dqdt = dtinv * (qp_0 - qn_0) * w_nat[0]
                                + dtinv * (qp_im1 - qn_im1) * w_nat[1]
                                + dtinv * (qp_im2 - qn_im2) * w_nat[2];
                            double dzetadx = dxinv * (htheta_0 + zb[ph_0] - htheta_im1 - zb[ph_w]);
                            //
                            // q-momentum
                            //
                            A.coeffRef(q_eq, 3 * ph_0 ) = w_nat[0] * theta * g * dzetadx + dxinv * theta * g * htheta_im12;
                            A.coeffRef(q_eq, 3 * ph_w ) = w_nat[1] * theta * g * dzetadx - dxinv * theta * g * htheta_im12;
                            A.coeffRef(q_eq, 3 * ph_ww) = w_nat[2] * theta * g * dzetadx;

                            A.coeffRef(q_eq, 3 * ph_0  + 1) = dtinv * w_nat[0];
                            A.coeffRef(q_eq, 3 * ph_w  + 1) = dtinv * w_nat[1];
                            A.coeffRef(q_eq, 3 * ph_ww + 1) = dtinv * w_nat[2];
                            // flow y
                            A.coeffRef(q_eq, 3 * ph_0  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_w  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ww + 2) = 0.0;
                            //
                            rhs[q_eq] = -(dqdt + g * htheta_im12 * dzetadx);
                            //
                            // continuity part (added and multiplied by +c_wave)
                            //
                            if (do_q_convection) { con_fac = c_wave + qp_im12 / hp_im12; }
                            // Contribution Delta h
                            A.coeffRef(q_eq, 3 * ph_0 ) += con_fac * dtinv * w_nat[0];
                            A.coeffRef(q_eq, 3 * ph_w ) += con_fac * dtinv * w_nat[1];
                            A.coeffRef(q_eq, 3 * ph_ww) += con_fac * dtinv * w_nat[2];
                            // Contribution Delta q
                            A.coeffRef(q_eq, 3 * ph_0  + 1) += con_fac * dxinv * theta;
                            A.coeffRef(q_eq, 3 * ph_w  + 1) += con_fac * -dxinv * theta;
                            A.coeffRef(q_eq, 3 * ph_ww + 1) += 0.0;
                            // Contribution Delta r
                            A.coeffRef(q_eq, 3 * ph_0  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_w  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ww + 2) = 0.0;
                            //
                            rhs[q_eq] += -con_fac * (dhdt + dqdx);
                            //
                            // r-momentum (tangential equation)
                            //
                            A.coeffRef(r_eq, 3 * ph_0 ) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_w ) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ww) = 0.0;
                            // flow x
                            A.coeffRef(r_eq, 3 * ph_0  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_w  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ww + 1) = 0.0;
                            // flow y
                            A.coeffRef(r_eq, 3 * ph_0  + 2) = 1.0;
                            A.coeffRef(r_eq, 3 * ph_w  + 2) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ww + 2) = 0.0;
                            //
                            rhs[r_eq] = 0.0;
                        }
                    }
                }
            }
            if (true)  // ssouth boundary
            {
                int j = 0;
                for (int i = 1; i < nx - 1; ++i)
                {
                    int ph_0  = p_index(i, j    , ny);  // central point of control volume
                    int ph_n  = p_index(i, j + 1, ny);  //
                    int ph_nn = p_index(i, j + 2, ny);  //

                    int c_eq = 3 * ph_0;
                    int q_eq = c_eq + 1;
                    int r_eq = c_eq + 2;

                    hn_0 = hn[ph_0];
                    hn_jp1 = hn[ph_n];
                    hn_jp2 = hn[ph_nn];
                    hp_0 = hp[ph_0];
                    hp_jp1 = hp[ph_n];
                    hp_jp2 = hp[ph_nn];
                    htheta_0 = theta * hp_0 + (1.0 - theta) * hn_0;
                    htheta_jp1 = theta * hp_jp1 + (1.0 - theta) * hn_jp1;
                    htheta_jp2 = theta * hp_jp2 + (1.0 - theta) * hn_jp2;

                    hn_jp12 = w_nat[0] * hn_0 + w_nat[1] * hn_jp1 + w_nat[2] * hn_jp2;
                    hp_jp12 = w_nat[0] * hp_0 + w_nat[1] * hp_jp1 + w_nat[2] * hp_jp2;
                    htheta_jp12 = w_nat[0] * htheta_0 + w_nat[1] * htheta_jp1 + w_nat[2] * htheta_jp2;

                    qn_0 = qn[ph_0];
                    qn_jp1 = qn[ph_n];
                    qn_jp2 = qn[ph_nn];
                    qp_0 = qp[ph_0];
                    qp_jp1 = qp[ph_n];
                    qp_jp2 = qp[ph_nn];
                    qtheta_0 = theta * qp_0 + (1.0 - theta) * qn_0;
                    qtheta_jp1 = theta * qp_jp1 + (1.0 - theta) * qn_jp1;
                    qtheta_jp2 = theta * qp_jp2 + (1.0 - theta) * qn_jp2;

                    qn_jp12 = w_nat[0] * qn_0 + w_nat[1] * qn_jp1 + w_nat[2] * qn_jp2;
                    qp_jp12 = w_nat[0] * qp_0 + w_nat[1] * qp_jp1 + w_nat[2] * qp_jp2;
                    qtheta_jp12 = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_jp1 + w_nat[2] * qtheta_jp2;

                    rn_0 = rn[ph_0];
                    rn_jp1 = rn[ph_n];
                    rn_jp2 = rn[ph_nn];
                    rp_0 = rp[ph_0];
                    rp_jp1 = rp[ph_n];
                    rp_jp2 = rp[ph_nn];
                    rtheta_0 = theta * rp_0 + (1.0 - theta) * rn_0;
                    rtheta_jp1 = theta * rp_jp1 + (1.0 - theta) * rn_jp1;
                    rtheta_jp2 = theta * rp_jp2 + (1.0 - theta) * rn_jp2;

                    rn_jp12 = w_nat[0] * rn_0 + w_nat[1] * rn_jp1 + w_nat[2] * rn_jp2;
                    rp_jp12 = w_nat[0] * rp_0 + w_nat[1] * rp_jp1 + w_nat[2] * rp_jp2;
                    rtheta_jp12 = w_nat[0] * rtheta_0 + w_nat[1] * rtheta_jp1 + w_nat[2] * rtheta_jp2;

                    double zb_jp12 = w_nat[0] * zb[ph_0] + w_nat[1] * zb[ph_n] + w_nat[2] * zb[ph_nn];
                    double h_infty = -zb_jp12;
                    double c_wave = std::sqrt(g * h_infty);

                    if (bc_type[BC_SOUTH] == "dirichlet" || bc_absorbing[BC_SOUTH] == false)
                    {
                        //
                        // continuity equation
                        //
                        A.coeffRef(c_eq, 3 * ph_0 ) = -1.0;
                        A.coeffRef(c_eq, 3 * ph_n ) = 1.0;
                        A.coeffRef(c_eq, 3 * ph_nn) = 0.0;
                        rhs[c_eq] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(q_eq, 3 * ph_0  + 1) = -1.0;
                        A.coeffRef(q_eq, 3 * ph_n  + 1) = 1.0;
                        A.coeffRef(q_eq, 3 * ph_nn + 1) = 0.0;
                        rhs[q_eq] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(r_eq, 3 * ph_0  + 2) = 0.5;
                        A.coeffRef(r_eq, 3 * ph_n  + 2) = 0.5;
                        A.coeffRef(r_eq, 3 * ph_nn + 2) = 0.0;
                        rhs[r_eq] = 0.0;
                    }
                    else
                    {
                        if (bc_absorbing[BC_SOUTH])
                        {
                            double con_fac = c_wave;
                            if (do_r_convection) { con_fac = c_wave + rp_jp12 / hp_jp12; }
                            if (bc_type[BC_SOUTH] == "mooiman")
                            {
                                //
                                // continuity
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = w_nat[0] * c_wave;
                                A.coeffRef(c_eq, 3 * ph_n ) = w_nat[1] * c_wave;
                                A.coeffRef(c_eq, 3 * ph_nn) = w_nat[2] * c_wave;
                                // Contribution Delta q
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_n  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_nn + 1) = 0.0;
                                // Contribution Delta r
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = w_nat[0];
                                A.coeffRef(c_eq, 3 * ph_n  + 2) = w_nat[1];
                                A.coeffRef(c_eq, 3 * ph_nn + 2) = w_nat[2];
                                //
                                rhs[c_eq] = -(rp_jp12 + c_wave * (hp_jp12 - h_infty));
                            }
                            if (bc_type[BC_SOUTH] == "borsboom")
                            {
                                //
                                // Essential boundary condition
                                // ----------------------------
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = dtinv * con_fac * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_n ) = dtinv * con_fac * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_nn) = dtinv * con_fac * w_ess[2];
                                // Contribution Delta q
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_n  + 1) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_nn + 1) = 0.0;
                                // Contribution Delta r
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = dtinv * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_n  + 2) = dtinv * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_nn + 2) = dtinv * w_ess[2];
                                //
                                double dhdt = dtinv * (hp_jp12 - hn_jp12);
                                double drdt = dtinv * (rp_jp12 - rn_jp12);
                                rhs[c_eq] = -(drdt + con_fac * dhdt);

                                double corr_term = 0.0;
                                if (bc_vars[BC_SOUTH] == "zeta")
                                {
                                    corr_term = - eps_bc_corr * (hp_ip12 - (bc[BC_SOUTH] - zb_jp12));
                                    rhs[c_eq] += corr_term;
                                }
                                if (bc_vars[BC_WEST] == "q")
                                {
                                    corr_term =  - eps_bc_corr * (qp_ip12 - bc[BC_SOUTH]);
                                    rhs[c_eq] += corr_term;
                                }
                            }
                            //
                            // natural boundary condition (second equation)
                            // --------------------------
                            //
                            if (do_r_convection) { con_fac = c_wave - rp_jp12 / hp_jp12; }
                            double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                                + dtinv * (hp_jp1 - hn_jp1) * w_nat[1]
                                + dtinv * (hp_jp2 - hn_jp2) * w_nat[2];
                            double drdy = dyinv * (rtheta_jp1 - rtheta_0);
                            double drdt = dtinv * (rp_0 - rn_0) * w_nat[0]
                                + dtinv * (rp_jp1 - rn_jp1) * w_nat[1]
                                + dtinv * (rp_jp2 - rn_jp2) * w_nat[2];
                            double dzetady = dyinv * (htheta_jp1 + zb[ph_n] - htheta_0 - zb[ph_0]);
                            //
                            // second equation (q == 0)
                            //
                            A.coeffRef(q_eq, 3 * ph_0) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_n) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_nn) = 0.0;
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 1) = 1.0;
                            A.coeffRef(q_eq, 3 * ph_n  + 1) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_nn + 1) = 0.0;
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_n  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_nn + 2) = 0.0;
                            //
                            rhs[q_eq] = 0.0;
                            //
                            // third equation (mom_eq - c_wave con_eq)
                            //
                            // Contribution Delta h
                            A.coeffRef(r_eq, 3 * ph_0 ) = w_nat[0] * theta * g * dzetady - dyinv * theta * g * htheta_jp12;
                            A.coeffRef(r_eq, 3 * ph_n ) = w_nat[1] * theta * g * dzetady + dyinv * theta * g * htheta_jp12;
                            A.coeffRef(r_eq, 3 * ph_nn) = w_nat[2] * theta * g * dzetady;
                            //
                            A.coeffRef(r_eq, 3 * ph_0  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_n  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_nn + 1) = 0.0;
                            //
                            A.coeffRef(r_eq, 3 * ph_0  + 2) = dtinv * w_nat[0];
                            A.coeffRef(r_eq, 3 * ph_n  + 2) = dtinv * w_nat[1];
                            A.coeffRef(r_eq, 3 * ph_nn + 2) = dtinv * w_nat[2];
                            //
                            rhs[r_eq] = -(drdt + g * htheta_jp12 * dzetady);

                            // Contribution Delta h
                            A.coeffRef(r_eq, 3 * ph_0 ) += -con_fac * dtinv * w_nat[0];
                            A.coeffRef(r_eq, 3 * ph_n ) += -con_fac * dtinv * w_nat[1];
                            A.coeffRef(r_eq, 3 * ph_nn) += -con_fac * dtinv * w_nat[2];
                            // Contribution Delta q
                            A.coeffRef(r_eq, 3 * ph_0  + 1) += 0.0;
                            A.coeffRef(r_eq, 3 * ph_n  + 1) += 0.0;
                            A.coeffRef(r_eq, 3 * ph_nn + 1) += 0.0;
                            // Contribution Delta q
                            A.coeffRef(r_eq, 3 * ph_0  + 2) += -con_fac * -dyinv * theta;
                            A.coeffRef(r_eq, 3 * ph_n  + 2) += -con_fac * dyinv * theta;
                            A.coeffRef(r_eq, 3 * ph_nn + 2) += 0.0;
                            //
                            rhs[r_eq] += con_fac * (dhdt + drdy);
                        }
                    }
                }
            }
            if (true)  // wwest boundary(2D)
            {
                int i = 0;
                for (int j = 1; j < ny - 1; ++j)
                {
                    int ph_0  = p_index(i    , j, ny);  // central point of control volume
                    int ph_e  = p_index(i + 1, j, ny); 
                    int ph_ee = p_index(i + 2, j, ny); 
                    
                    int c_eq = 3 * ph_0;   // continuity equation
                    int q_eq = c_eq + 1;   // q-momentum equation
                    int r_eq = c_eq + 2;   // r-momentum equation

                    hn_0 = hn[ph_0];
                    hn_ip1 = hn[ph_e];
                    hn_ip2 = hn[ph_ee];
                    hp_0 = hp[ph_0];
                    hp_ip1 = hp[ph_e];
                    hp_ip2 = hp[ph_ee];
                    htheta_0 = theta * hp_0 + (1.0 - theta) * hn_0;
                    htheta_ip1 = theta * hp_ip1 + (1.0 - theta) * hn_ip1;
                    htheta_ip2 = theta * hp_ip2 + (1.0 - theta) * hn_ip2;

                    hn_ip12 = w_nat[0] * hn_0 + w_nat[1] * hn_ip1 + w_nat[2] * hn_ip2;
                    hp_ip12 = w_nat[0] * hp_0 + w_nat[1] * hp_ip1 + w_nat[2] * hp_ip2;
                    htheta_ip12 = w_nat[0] * htheta_0 + w_nat[1] * htheta_ip1 + w_nat[2] * htheta_ip2;

                    qn_0 = qn[ph_0];
                    qn_ip1 = qn[ph_e];
                    qn_ip2 = qn[ph_ee];
                    qp_0 = qp[ph_0];
                    qp_ip1 = qp[ph_e];
                    qp_ip2 = qp[ph_ee];
                    qtheta_0 = theta * qp_0 + (1.0 - theta) * qn_0;
                    qtheta_ip1 = theta * qp_ip1 + (1.0 - theta) * qn_ip1;
                    qtheta_ip2 = theta * qp_ip2 + (1.0 - theta) * qn_ip2;

                    qn_ip12 = w_nat[0] * qn_0 + w_nat[1] * qn_ip1 + w_nat[2] * qn_ip2;
                    qp_ip12 = w_nat[0] * qp_0 + w_nat[1] * qp_ip1 + w_nat[2] * qp_ip2;
                    qtheta_ip12 = w_nat[0] * qtheta_0 + w_nat[1] * qtheta_ip1 + w_nat[2] * qtheta_ip2;

                    rn_0 = rn[ph_0];
                    rn_ip1 = rn[ph_e];
                    rn_ip2 = rn[ph_ee];
                    rp_0 = rp[ph_0];
                    rp_ip1 = rp[ph_e];
                    rp_ip2 = rp[ph_ee];
                    rtheta_0 = theta * rp_0 + (1.0 - theta) * rn_0;
                    rtheta_ip1 = theta * rp_ip1 + (1.0 - theta) * rn_ip1;
                    rtheta_ip2 = theta * rp_ip2 + (1.0 - theta) * rn_ip2;

                    rn_ip12 = w_nat[0] * rn_0 + w_nat[1] * rn_ip1 + w_nat[2] * rn_ip2;
                    rp_ip12 = w_nat[0] * rp_0 + w_nat[1] * rp_ip1 + w_nat[2] * rp_ip2;
                    rtheta_ip12 = w_nat[0] * rtheta_0 + w_nat[1] * rtheta_ip1 + w_nat[2] * rtheta_ip2;

                    double zb_ip12 = w_nat[0] * zb[ph_0] + w_nat[1] * zb[ph_e] + w_nat[2] * zb[ph_ee];
                    double h_infty = -zb_ip12;
                    double c_wave = std::sqrt(g * h_infty);

                    if (bc_type[BC_WEST] == "dirichlet" || bc_absorbing[BC_WEST] == false)
                    {
                        //
                        // continuity equation
                        //
                        A.coeffRef(c_eq, 3 * ph_0) = -1.0;
                        A.coeffRef(c_eq, 3 * ph_e) = 1.0;
                        A.coeffRef(c_eq, 3 * ph_ee) = 0.0;
                        rhs[c_eq] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(q_eq, 3 * ph_0 + 1) = 0.5;
                        A.coeffRef(q_eq, 3 * ph_e + 1) = 0.5;
                        A.coeffRef(q_eq, 3 * ph_ee + 1) = 0.0;
                        rhs[q_eq] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(r_eq, 3 * ph_0 + 2) = -1.0;
                        A.coeffRef(r_eq, 3 * ph_e + 2) = 1.0;
                        A.coeffRef(r_eq, 3 * ph_ee + 2) = 0.0;
                        rhs[r_eq] = 0.0;
                    }
                    else
                    {
                        if (bc_absorbing[BC_WEST])
                        {
                            double con_fac = c_wave;
                            if (do_q_convection) { con_fac = c_wave + qp_ip12 / hp_ip12; }
                            if (bc_type[BC_WEST] == "mooiman")
                            {
                                //
                                // first equation
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = w_nat[0] * c_wave;
                                A.coeffRef(c_eq, 3 * ph_e ) = w_nat[1] * c_wave;
                                A.coeffRef(c_eq, 3 * ph_ee) = w_nat[2] * c_wave;
                                //
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = w_nat[0];
                                A.coeffRef(c_eq, 3 * ph_e  + 1) = w_nat[1];
                                A.coeffRef(c_eq, 3 * ph_ee + 1) = w_nat[2];
                                //
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_e  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ee + 2) = 0.0;
                                //
                                rhs[c_eq] = -(qp_ip12 + c_wave * (hp_ip12 - h_infty));
                            }
                            if (bc_type[BC_WEST] == "borsboom")
                            {
                                //
                                // Essential boundary condition
                                // ----------------------------
                                // momentum + c_wave * continuity
                                //
                                A.coeffRef(c_eq, 3 * ph_0 ) = dtinv * con_fac * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_e ) = dtinv * con_fac * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ee) = dtinv * con_fac * w_ess[2];
                                // Contribution Delta q
                                A.coeffRef(c_eq, 3 * ph_0  + 1) = dtinv * w_ess[0];
                                A.coeffRef(c_eq, 3 * ph_e  + 1) = dtinv * w_ess[1];
                                A.coeffRef(c_eq, 3 * ph_ee + 1) = dtinv * w_ess[2];
                                // No contribution for Delta r
                                A.coeffRef(c_eq, 3 * ph_0  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_e  + 2) = 0.0;
                                A.coeffRef(c_eq, 3 * ph_ee + 2) = 0.0;
                                //
                                double dhdt = dtinv * (hp_ip12 - hn_ip12);
                                double dqdt = dtinv * (qp_ip12 - qn_ip12);
                                rhs[c_eq] = -(dqdt + con_fac * dhdt);

                                double corr_term = 0.0;
                                if (bc_vars[BC_WEST] == "zeta")
                                {
                                    //dhdt = dtinv * (bc[BC_WEST] - hn_ip12);
                                    //A.coeffRef(c_eq, ph_0 ) += dtinv * w_ess[0] + eps_bc_corr * w_ess[0];
                                    //A.coeffRef(c_eq, ph_e ) += dtinv * w_ess[1] + eps_bc_corr * w_ess[1];
                                    //A.coeffRef(c_eq, ph_ee) += dtinv * w_ess[2] + eps_bc_corr * w_ess[2];
                                    //corr_term = -dhdt - eps_bc_corr * (hp_ip12 - (bc[BC_WEST] - zb_ip12));
                                    corr_term = - eps_bc_corr * (hp_ip12 - (bc[BC_WEST] - zb_ip12));
                                    rhs[c_eq] += corr_term;
                                }
                                if (bc_vars[BC_WEST] == "q")
                                {
                                    //A.coeffRef(c_eq, ph_0  + 1) += 0.0 * dtinv * w_ess[0] + eps_bc_corr * theta * w_ess[0];
                                    //A.coeffRef(c_eq, ph_e  + 1) += 0.0 * dtinv * w_ess[1] + eps_bc_corr * theta * w_ess[1];
                                    //A.coeffRef(c_eq, ph_ee + 1) += 0.0 * dtinv * w_ess[2] + eps_bc_corr * theta * w_ess[2];
                                    //corr_term = -2. * dqdt - eps_bc_corr * (qp_ip12 - bc[BC_WEST]);
                                    corr_term =  - eps_bc_corr * (qp_ip12 - bc[BC_WEST]);
                                    rhs[c_eq] += corr_term;
                                }
                            }
                            //
                            // natural boundary condition (second equation)
                            // --------------------------
                            // momentum - c_wave * continuity
                            //
                            if (do_q_convection) { con_fac = c_wave - qp_ip12 / hp_ip12; }
                            double dhdt = dtinv * (hp_0 - hn_0) * w_nat[0]
                                + dtinv * (hp_ip1 - hn_ip1) * w_nat[1]
                                + dtinv * (hp_ip2 - hn_ip2) * w_nat[2];
                            double dqdx = dxinv * (qtheta_ip1 - qtheta_0);
                            double dqdt = dtinv * (qp_0 - qn_0) * w_nat[0]
                                + dtinv * (qp_ip1 - qn_ip1) * w_nat[1]
                                + dtinv * (qp_ip2 - qn_ip2) * w_nat[2];
                            double dzetadx = dxinv * (htheta_ip1 + zb[ph_e] - htheta_0 - zb[ph_0]);
                            //
                            // momentum part dq/dt + gh d(zeta)/dx
                            //
                            // Contribution Delta h
                            A.coeffRef(q_eq, 3 * ph_0 ) = w_nat[0] * theta * g * dzetadx - dxinv * theta * g * htheta_ip12;
                            A.coeffRef(q_eq, 3 * ph_e ) = w_nat[1] * theta * g * dzetadx + dxinv * theta * g * htheta_ip12;
                            A.coeffRef(q_eq, 3 * ph_ee) = w_nat[2] * theta * g * dzetadx;
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 1) = dtinv * w_nat[0];
                            A.coeffRef(q_eq, 3 * ph_e  + 1) = dtinv * w_nat[1];
                            A.coeffRef(q_eq, 3 * ph_ee + 1) = dtinv * w_nat[2];
                            //
                            A.coeffRef(q_eq, 3 * ph_0  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_e  + 2) = 0.0;
                            A.coeffRef(q_eq, 3 * ph_ee + 2) = 0.0;
                            //
                            rhs[q_eq] = -(dqdt + g * htheta_ip12 * dzetadx);
                            //
                            // continuity part (added and multiplied by -c_wave)
                            //
                            // Contribution Delta h
                            A.coeffRef(q_eq, 3 * ph_0 ) += -con_fac * dtinv * w_nat[0];
                            A.coeffRef(q_eq, 3 * ph_e ) += -con_fac * dtinv * w_nat[1];
                            A.coeffRef(q_eq, 3 * ph_ee) += -con_fac * dtinv * w_nat[2];
                            // Contribution Delta q
                            A.coeffRef(q_eq, 3 * ph_0  + 1) += -con_fac * -dxinv * theta;
                            A.coeffRef(q_eq, 3 * ph_e  + 1) += -con_fac * dxinv * theta;
                            A.coeffRef(q_eq, 3 * ph_ee + 1) += 0.0;
                            // Contribution Delta r
                            A.coeffRef(q_eq, 3 * ph_0  + 2) += 0.0;
                            A.coeffRef(q_eq, 3 * ph_e  + 2) += 0.0;
                            A.coeffRef(q_eq, 3 * ph_ee + 2) += 0.0;
                            
                            rhs[q_eq] += con_fac * (dhdt + dqdx);
                            //
                            // r-momentum (r == 0)
                            //
                            A.coeffRef(r_eq, 3 * ph_0 ) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_e ) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ee) = 0.0;
                            // flow x
                            A.coeffRef(r_eq, 3 * ph_0  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_e  + 1) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ee + 1) = 0.0;
                            // flow y
                            A.coeffRef(r_eq, 3 * ph_0  + 2) = 1.0;
                            A.coeffRef(r_eq, 3 * ph_e  + 2) = 0.0;
                            A.coeffRef(r_eq, 3 * ph_ee + 2) = 0.0;
                            //
                            rhs[r_eq] = 0.0;
                        }
                    }
                }
            }
            //
            //corner nodes
            //
            if (true)  // NE-corner
            {
                int i = nx - 1;
                int j = ny - 1;
                int ph = 3 * p_index(i, j, ny);  // continuity equation
                int pq = ph + 1;  // q-momentum equation
                int pr = ph + 2;  // r-momentum equation
                int p1 = 3 * p_index(i - 1, j , ny);
                int p2 = 3 * p_index(i, j - 1, ny);

                A.coeffRef(ph, p1) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, p2) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, p1 + 1) = 0.5;
                A.coeffRef(pq, pq    ) = -1.0;
                A.coeffRef(pq, p2 + 1) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, p1 + 2) = 0.5;
                A.coeffRef(pr, pr    ) = -1.0;
                A.coeffRef(pr, p2 + 2) = 0.5;
                rhs[pr] = 0.0;
            }
            if (true)  // SE-corner
            {
                int i = nx - 1;
                int j = 0;
                int ph = 3 * p_index(i, j, ny);  // continuity equation
                int pq = ph + 1;  // q-momentum equation
                int pr = ph + 2;  // r-momentum equation
                int p1 = 3 * p_index(i - 1, j, ny);
                int p2 = 3 * p_index(i, j + 1, ny);

                A.coeffRef(ph, p1) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, p2) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, p1 + 1) = 0.5;
                A.coeffRef(pq, pq    ) = -1.0;
                A.coeffRef(pq, p2 + 1) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, p1 + 2) = 0.5;
                A.coeffRef(pr, pr    ) = -1.0;
                A.coeffRef(pr, p2 + 2) = 0.5;
                rhs[pr] = 0.0;
            }
            if (true)  // SW-corner
            {
                int i = 0;
                int j = 0;
                int ph = 3 * p_index(i, j, ny);  // continuity equation
                int pq = ph + 1;  // q-momentum equation
                int pr = ph + 2;  // r-momentum equation
                int p1 = 3 * p_index(i + 1, j , ny);
                int p2 = 3 * p_index(i, j + 1, ny);

                A.coeffRef(ph, p1) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, p2) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, p1 + 1) = 0.5;
                A.coeffRef(pq, pq    ) = -1.0;
                A.coeffRef(pq, p2 + 1) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, p1 + 2) = 0.5;
                A.coeffRef(pr, pr    ) = -1.0;
                A.coeffRef(pr, p2 + 2) = 0.5;
                rhs[pr] = 0.0;
            }
            if (true)  // NW-corner
            {
                int i = 0;
                int j = ny - 1;
                int ph = 3 * p_index(i, j, ny);  // continuity equation
                int pq = ph + 1;  // q-momentum equation
                int pr = ph + 2;  // r-momentum equation
                int p1 = 3 * p_index(i + 1, j , ny);
                int p2 = 3 * p_index(i, j - 1, ny);

                A.coeffRef(ph, p1) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, p2) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, p1 + 1) = 0.5;
                A.coeffRef(pq, pq    ) = -1.0;
                A.coeffRef(pq, p2 + 1) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, p1 + 2) = 0.5;
                A.coeffRef(pr, pr    ) = -1.0;
                A.coeffRef(pr, p2 + 2) = 0.5;
                rhs[pr] = 0.0;
            }
            if (nst == 1 && iter == 0)
            {
                START_TIMER(BiCGStab_initialization);
            }
            else if (nst != 1)
            {
                START_TIMER(BiCGStab);
            }

            if (logging == "matrix")
            {
                log_file << "=== Matrix ============================================" << std::endl;
                log_file << std::showpos << std::setprecision(4) << std::scientific << A << std::endl;
                log_file << "=== Matrix diagonal====================================" << std::endl;
                Eigen::VectorXd diag = A.diagonal();
                for (int i = 0; i < diag.size(); ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << "Index " << i << ": " << diag[i] << std::endl;
                }
                //log_file << std::setprecision(8) << std::scientific << A.diagonal() << std::endl;
                log_file << "=== RHS ===============================================" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
                }
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
                    k = p_index(i, j ,ny);
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
            if (logging == "matrix")
            {
                log_file << "=== Solution ==========================================" << std::endl;
                for (int i = 0; i < 3 * nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << solution[i] << std::endl;
                }
                log_file << "=== hp, qp, rp ========================================" << std::endl;
                for (int i = 0; i < nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << hp[i] << ", " << qp[i] << ", " << rp[i] << std::endl;
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
            {
                int c = 1;
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
                log_file << "    ----    maximum number of iterations reached, probably not converged, at time: " <<  time << std::endl;
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
            his_values = { hn[p_a], hn[p_b], hn[p_c], hn[p_d], hn[centre], hn[n_bnd], hn[ne_bnd], hn[e_bnd], hn[se_bnd], hn[s_bnd], hn[sw_bnd], hn[w_bnd], hn[nw_bnd] };
            his_file->put_variable(his_h_name, nst_his, his_values);

            his_values = { qn[p_a], qn[p_b], qn[p_c], qn[p_d], qn[centre], qn[n_bnd], qn[ne_bnd], qn[e_bnd], qn[se_bnd], qn[s_bnd], qn[sw_bnd], qn[w_bnd], qn[nw_bnd] };
            his_file->put_variable(his_q_name, nst_his, his_values);
 
            his_values = { rn[p_a], rn[p_b], rn[p_c], rn[p_d], rn[centre], rn[n_bnd], rn[ne_bnd], rn[e_bnd], rn[se_bnd], rn[s_bnd], rn[sw_bnd], rn[w_bnd], rn[nw_bnd] };
            his_file->put_variable(his_r_name, nst_his, his_values);

            his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd] };
            his_file->put_variable(his_s_name, nst_his, his_values);

            his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd] };
            his_file->put_variable(his_u_name, nst_his, his_values);

            his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd] };
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
inline int p_index(int i, int j, int ny)
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
inline double dcdx_scvf_n(double c0, double c1, double c2, double c3)
{
    // dcdx normal at subcontrol volume edge
    return 0.25 * (3. * c0 - 3. * c1 + c2 - c3);
}
inline double dcdx_scv_t(double c0, double c1, double c2, double c3)
{
    // dcdy tangential at subcontrol volume edge
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

