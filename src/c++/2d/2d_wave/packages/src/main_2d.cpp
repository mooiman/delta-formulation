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

#define NOMINMAX 
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
#include <include/KDtree.hpp>
// for BiCGstab  solver
//#include <Eigen/Dense>
//#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/Core>

// for Algebraic Multigrid solver
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/profiler.hpp>

AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND()

#include "bed_level.h"
#include "bed_shear_stress.h"
#include "boundary_condition.h"
#include "build_matrix_pattern.h"
#include "cfts.h"
#include "compile_date_and_time.h"
#include "convection.h"
#include "data_input_struct.h"
#include "grid.h"
#include "initial_conditions.h"
#include "interpolations.h"
#include "perf_timer.h"
#include "matrix_assembly_boundaries.h"
#include "matrix_assembly_corners.h"
#include "matrix_assembly_interior.h"
#include "observation_stations.h"
#include "read_input_toml_file.h"
#include "regularization.h"
#include "ugrid2d.h"
#include "viscosity.h"

void GetArguments(long argc, char** argv, std::filesystem::path & file_name);

int write_used_input(struct _data_input data, std::ofstream & log_file);
inline size_t main_idx(size_t i, size_t j, size_t ny);

// Solve the linear wave equation
// Continuity equation: d(h)/dt + d(q)/dx = 0
// Momentum equation  : d(q)/dt + gh d(zeta)/dx + convection + bed_shear_stress - diffusivity = 0
// Momentum equation  : d(r)/dt + gh d(zeta)/dy + convection + bed_shear_stress - diffusivity = 0

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
    int solver_iterations = -1;
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

    struct _data_input input_data;
    input_data = read_toml_file(input_dir, toml_file_name);
    if (input_data.numerics.dt == 0.0) { stationary = true;  }
    if (stationary) { 
        input_data.boundary.treg = 0.0; 
        input_data.numerics.theta = 1.0;
    }
    if (input_data.initial.gauss_mu != -INFINITY)
    {
        input_data.initial.gauss_mu_x = input_data.initial.gauss_mu; 
        input_data.initial.gauss_mu_y = 0.0;  // only shift on x-axis
    }
    if (input_data.initial.gauss_sigma != -INFINITY)
    {
        input_data.initial.gauss_sigma_x = input_data.initial.gauss_sigma; 
        input_data.initial.gauss_sigma_y = input_data.initial.gauss_sigma; 
    }

    ss << "2d_wave";
    out_file = output_dir.string() + ss.str();
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
    log_file << "Executable compiled: " << compileDateTime() << std::endl;
    log_file << "Start time         : " << start_date_time << std::endl;
    log_file << "=== Input file =======================================" << std::endl;
    log_file << toml_file_name << std::endl;

    struct _mesh2d * mesh2d;
    SGRID* sgrid = new SGRID();
    status = sgrid->open(input_data.domain.full_grid_filename.string());
    if (status != 0) 
    {
        log_file << "Error: Failed to open grid file: " << input_data.domain.full_grid_filename.string() << std::endl;
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    status = sgrid->read();  // reading mesh
    mesh2d = sgrid->get_mesh_2d();

    size_t nx = mesh2d->node[0]->dims[0]; // including virtual points
    size_t ny = mesh2d->node[0]->dims[1]; // including virtual points

    std::vector<double>& x = mesh2d->node[0]->x; // reference to original x-coordinate
    std::vector<double>& y = mesh2d->node[0]->y; // reference to original y-coordinate

    BED_LEVEL * bed = new BED_LEVEL();
    status = bed->open(input_data.domain.full_bed_level_filename.string());
    if (status != 0)
    {
        std::cout << "Failed to open file: " << input_data.domain.full_bed_level_filename.string() << std::endl;
        log_file << "Failed to open file: " << input_data.domain.full_bed_level_filename.string() << std::endl;
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        //std::cin.ignore();
        exit(1);
    }
    status = bed->read(nx, ny);
    std::vector<double> zb_giv = bed->get_bed_level();

    //  Create kdtree, needed to locate the obsservation points
    std::vector<std::vector<double>> xy_points;
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::vector<double> point = {x[i], y[i]};
        xy_points.push_back(point);
    }
    KDTree xy_tree(xy_points);

    double dx = x[ny] - x[0];
    double dy = y[1] - y[0];
    double Lx = dx * mesh2d->face[0]->dims[0];
    double Ly = dy * mesh2d->face[0]->dims[1];

    double dxinv = 1./dx;                               // invers grid size [m]
    double dyinv = 1./dy;                               // invers grid size [m]
    size_t nxny = nx * ny;                                   // total number of nodes
    double dxdy = dx * dy ;                               // area of control volume
    //if (viscosity == 0.0)
    //{
    //    viscosity = 0.2 * std::sqrt(dx*dx + dy*dy);
    //}

    int total_time_steps = int((input_data.time.tstop - input_data.time.tstart) / input_data.numerics.dt) + 1;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    double dtpseuinv = 0.0;                                     // Inverse of dtpseu
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
    std::string solver_name("--- undefined ---");
    if (input_data.numerics.linear_solver == "bicgstab") { solver_name = "BiCGstab"; }
    if (input_data.numerics.linear_solver == "multigrid") { solver_name = "MultiGrid"; }

    std::string model_title("Linear wave equation");

    if (input_data.physics.do_q_equation && input_data.physics.do_r_equation)
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
    else if (input_data.physics.do_r_equation)
    {
        model_title = "Linear wave equation, only r-equation";
    }
    else
    {
        model_title = "No q- and no r-equation (=> no waves computed)";
    }
    model_title += ", " + solver_name;

    REGULARIZATION* regularization = new REGULARIZATION(input_data.numerics.iter_max, input_data.physics.g);

    log_file << "=== Used input variables ==============================" << std::endl;
    status = write_used_input(input_data, log_file);
    log_file << "=======================================================" << std::endl;

    // Copy input data to loccal data
    std::string logging = input_data.log.logging;

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
    double c_psi = input_data.numerics.c_psi;
    double eps_bicgstab = input_data.numerics.eps_bicgstab;
    double eps_newton = input_data.numerics.eps_newton;
    double theta = input_data.numerics.theta;
    double iter_max = input_data.numerics.iter_max;
    std::string linear_solver = input_data.numerics.linear_solver;
    bool regularization_iter = input_data.numerics.regularization_iter;
    bool regularization_init = input_data.numerics.regularization_init;

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
    std::vector<double> s_giv(nxny, 0.);                  // water level, given
    std::vector<double> u_giv(nxny, 0.);                  // u-velocity, given
    std::vector<double> v_giv(nxny, 0.);                  // v-velocity, given
    std::vector<double> s(nxny, 0.);                      // water level, needed for post-processing
    std::vector<double> u(nxny, 0.);                      // u-velocity, needed for post-processing
    std::vector<double> v(nxny, 0.);                      // v-velocity, needed for post-processing
    std::vector<double> visc_given(nxny, visc_const);     // Initialize viscosity array with given value
    std::vector<double> visc_reg(nxny, visc_const);       // Initialize given viscosity array with regularized value
    std::vector<double> visc(nxny, visc_const);           // Viscosity array used for computation, adjusted for cell peclet number
    std::vector<double> delta_h(nxny, 0.);
    std::vector<double> delta_q(nxny, 0.);
    std::vector<double> delta_r(nxny, 0.);
    std::vector<double> post_q;              // needed for postprocessing; bed shear stress; convection;
    std::vector<double> post_r;              // needed for postprocessing; bed shear stress; convection;
    if (do_bed_shear_stress || do_convection || do_viscosity)
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
    //w_ess[0] = w_nat[0];
    //w_ess[1] = w_nat[1];
    //w_ess[2] = w_nat[2];

    //initialize water level
    std::cout << "Initialisation" << std::endl;
    status = 0;
    double min_zb = *std::min_element(zb_giv.begin(), zb_giv.end());

    log_file << "Nodes   : " << nx << "x" << ny << "=" << nxny << std::endl;
    log_file << "Elements: " << (nx - 1) << "x" << (ny - 1) << "=" << (nx - 1) * (ny - 1) << std::endl;
    log_file << "Volumes : " << (nx - 2) << "x" << (ny - 2) << "=" << (nx - 2) * (ny - 2) << std::endl;
    log_file << "CFL (2D): " << std::sqrt(g * std::abs(min_zb)) * dt * std::sqrt(( 1./(dx*dx) + 1./(dy*dy))) << std::endl;
    log_file << "CFL (1D): " << std::sqrt(g * std::abs(min_zb)) * dt /dx << std::endl;
    if (do_viscosity) { log_file << "Ddt/d^2 : " << visc_const * dt * std::sqrt(( 1./(dx*dx) + 1./(dy*dy))) << std::endl; }
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
        regularization->given_function(visc_reg, psi_11, psi_22, eq8, visc_given, nx, ny, dx, dy, c_psi, log_file);
        for (size_t i = 0; i < visc_reg.size(); ++i)
        {
            visc[i] = visc_reg[i] + std::sqrt(psi_11[i] * psi_11[i] + psi_22[i] * psi_22[i]);
        }

        STOP_TIMER(Regularization_init);
    }
    else
    {
        for (size_t i = 0; i < zb_giv.size(); ++i)
        {
            zb[i] = zb_giv[i];
            visc_reg[i] = visc_given[i];
        }
    }
    for (size_t k = 0; k < zb_giv.size(); ++k)
    {
        hn[k] = s_giv[k] - zb[k];  // Initial water depth
        qn[k] = hn[k] * u_giv[k];  // Initial q=hu -velocity
        rn[k] = hn[k] * v_giv[k];  // Initial r=hv -velocity
        // 
        hp[k] = hn[k]; 
        qp[k] = qn[k]; 
        rp[k] = rn[k]; 
    }

    for (size_t i = 0; i < nx*ny; ++i)
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
    status = map_file->add_edge_nodes(nx, ny);
    status = map_file->add_face_nodes(x, y, mesh2d->node[0]->fill_value, nx, ny);
    status = map_file->add_node_coords(x, y, mesh2d->node[0]->fill_value);
    status = map_file->add_edge_coords(x, y, mesh2d->node[0]->fill_value);
    status = map_file->add_face_mass_centres(x, y, mesh2d->node[0]->fill_value, nx, ny);
    status = map_file->add_face_area(x, y, mesh2d->node[0]->fill_value, nx, ny);

    status = map_file->add_time_series();

    std::vector<std::string> dim_names;
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
    std::string map_visc_q_name("viscosity_q");
    std::string map_visc_r_name("viscosity_r");

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
    if (do_viscosity)
    {
        status = map_file->add_variable(map_visc_q_name, dim_names, "", "Viscosity (x)", "m2 s-1", "mesh2D", "node");
        status = map_file->add_variable(map_visc_r_name, dim_names, "", "Viscosity (y)", "m2 s-1", "mesh2D", "node");
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
        bed_shear_stress_post_rhs(post_q, post_r, hn, qn, rn, cf, nx, ny);
        map_file->put_time_variable(map_beds_q_name, nst_map, post_q);
        map_file->put_time_variable(map_beds_r_name, nst_map, post_r);
    }
    if (do_convection)
    {
        convection_post_rhs(post_q, post_r, hn, qn, rn, dx, dy, nx, ny);
        map_file->put_time_variable(map_conv_q_name, nst_map, post_q);
        map_file->put_time_variable(map_conv_r_name, nst_map, post_r);
    }
    if (do_viscosity)
    {
        viscosity_post_rhs(post_q, post_r, hn, qn, rn, visc, dx, dy, nx, ny);
        map_file->put_time_variable(map_visc_q_name, nst_map, post_q);
        map_file->put_time_variable(map_visc_r_name, nst_map, post_r);
    }
    STOP_TIMER(Writing map-file);

    // End define map file
    ////////////////////////////////////////////////////////////////////////////
    // Create time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

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

    std::vector<std::string> obs_station_names;
    status = def_observation_stations(obs_station_names, input_data.obs_points, xy_tree, x, y, Lx, Ly, dx, dy, nx, ny);
    {
        std::vector<double> x_obs;
        std::vector<double> y_obs;
        for (size_t i = 0; i < input_data.obs_points.size(); ++i)
        {
            auto obs = input_data.obs_points[i];
            x_obs.push_back(obs.x);
            y_obs.push_back(obs.y);
        }

        his_file->add_stations(obs_station_names, x_obs, y_obs);
    }
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
    std::vector<double> his_values;
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(hn[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_h_name, nst_his, his_values);

    his_values.clear();
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(qn[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_q_name, nst_his, his_values);

    his_values.clear();
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(rn[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_r_name, nst_his, his_values);

    his_values.clear();
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(s[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_s_name, nst_his, his_values);

    his_values.clear();
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(u[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_u_name, nst_his, his_values);

    his_values.clear();
    for (size_t i = 0; i < input_data.obs_points.size(); ++i)
    {
        his_values.push_back(v[input_data.obs_points[i].idx]);
    }
    his_file->put_variable(his_v_name, nst_his, his_values);

    std::string his_newton_iter_name("his_newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "iterations", "Newton iteration", "-");
    his_values.clear();
    his_values = { 0.0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    std::string his_LinSolver_iter_name("his_LinSolver_iterations");
    his_file->add_variable_without_location(his_LinSolver_iter_name, "iterations", "LinSolver iteration", "-");
    his_values.clear();
    his_values = { 0.0 };
    his_file->put_variable(his_LinSolver_iter_name, nst_his, his_values);

    std::string his_LinSolver_iter_error_name("LinSolver_iteration_error");
    his_file->add_variable_without_location(his_LinSolver_iter_error_name, "iteration_error", "LinSolver iteration error", "-");
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
    if (logging == "pattern")
    {
        log_file << "=== Matrix build matrix pattern =======================" << std::endl;
        for (size_t i = 0; i < 3 * nxny; ++i)
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
    // 
    // AMGCL setup
    using T = double;
    using Backend = amgcl::backend::builtin<T>;
    using AMG = amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    >;
    using Solver = amgcl::make_solver<AMG, amgcl::solver::bicgstab<Backend>>;
    Solver::params prm;
    prm.solver.tol = eps_bicgstab;
    prm.solver.maxiter = 1000;
               

    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        time = dt * double(nst);

        int select = 1;  // Constant boundary condition
        std::vector<double> bc(4, 0.0);
        boundary_condition(bc[BC_NORTH], bc_vals[BC_NORTH], time, treg, select);
        boundary_condition(bc[BC_EAST ], bc_vals[BC_EAST ], time, treg, select);
        boundary_condition(bc[BC_SOUTH], bc_vals[BC_SOUTH], time, treg, select);
        boundary_condition(bc[BC_WEST ], bc_vals[BC_WEST ], time, treg, select);

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
            for (size_t k = 0; k < nxny; ++k)
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

            if (do_linear_waves)
            {
            // row 0, 1, 2; sw-corner
            // row 3,..., 3 * (ny - 1) - 1: west boundary
            // row 3 * (ny - 1), +1, +2; nw-corner
            // row std::fmod(col , 3 * ny), +1, +2; south boundary 
            // row std::fmod(col + 3, 3 * ny), +1, +2; north boundary
            // row 3 * (nx - 1) * ny, +1, +2; se-corner
            // row 3 * (nx - 1) * ny + 3, ..., 3 * nx * ny - 3 - 1: east boundary
            // row 3 * nx * ny - 3, +1, +2 : ne-corner

           // south-west corner
            for (size_t row = 0; row < 3; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status =  corner_south_west(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // west boundary
            for (size_t row = 3; row < 3 * (ny - 1); row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status =  boundary_west(values, row, c_eq, q_eq, r_eq, rhs, 
                                        dtinv, dxinv, theta, g, eps_bc_corr, 
                                        stationary, do_convection, do_bed_shear_stress, do_viscosity,
                                        dx, dy, nx, ny,
                                        hn, qn, rn,
                                        hp, qp, rp,
                                        htheta, qtheta, rtheta,
                                        zb, cf, 
                                        bc_type,  bc_vars, BC_WEST, bc,
                                        w_nat, w_ess);
            }
            // north-west corner
            for (size_t row = 3 * (ny - 1); row < 3 * ny; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                 status =  corner_north_west(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // interior with south and north boundary
            for (size_t row = 3 * ny; row < 3 * (nx - 1) * ny; row += 3) 
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status =  interior(values, row, c_eq, q_eq, r_eq, rhs, 
                                    dtinv, dxinv, theta, g, do_convection, nx, ny,
                                    x, y,
                                    hn, qn, rn,
                                    hp, qp, rp,
                                    htheta, qtheta, rtheta,
                                    zb, dx, dy, dxdy, mass);

                if (std::fmod(row , 3 * ny) == 0) {
                    // south boundary, over write coefficients
                    status = boundary_south(values, row, c_eq, q_eq, r_eq, rhs, 
                                            dtinv, dxinv, theta, g, eps_bc_corr, 
                                            stationary, do_convection, do_bed_shear_stress, do_viscosity,
                                            dx, dy, nx, ny,
                                            hn, qn, rn,
                                            hp, qp, rp,
                                            htheta, qtheta, rtheta,
                                            zb, cf, 
                                            bc_type,  bc_vars, BC_SOUTH, bc,
                                            w_nat, w_ess);
                } 
                if (std::fmod(row + 3, 3 * ny) == 0) {
                    // north boundary, over write coefficients
                    status = boundary_north(values, row, c_eq, q_eq, r_eq, rhs, 
                                            dtinv, dxinv, theta, g, eps_bc_corr, 
                                            stationary, do_convection, do_bed_shear_stress, do_viscosity,
                                            dx, dy, nx, ny,
                                            hn, qn, rn,
                                            hp, qp, rp,
                                            htheta, qtheta, rtheta,
                                            zb, cf,
                                            bc_type,  bc_vars, BC_NORTH, bc,
                                            w_nat, w_ess);
                }
            }
            // south-east corner
            for (size_t row = 3 * (nx - 1) * ny; row < 3 * (nx - 1) * ny + 3; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status = corner_south_east(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            // east boundary
            for (size_t row = 3 * (nx - 1) * ny + 3; row < 3 * nx * ny - 3; row += 3) 
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];

                status = boundary_east(values, row, c_eq, q_eq, r_eq, rhs, 
                                       dtinv, dxinv, theta, g, eps_bc_corr, 
                                       stationary, do_convection, do_bed_shear_stress, do_viscosity,
                                       dx, dy, nx, ny,
                                       hn, qn, rn,
                                       hp, qp, rp,
                                       htheta, qtheta, rtheta,
                                       zb, cf,
                                       bc_type,  bc_vars, BC_EAST, bc,
                                       w_nat, w_ess);
            }
            // north-east corner
            for (size_t row = 3 * nx * ny - 3; row < 3 * nx * ny; row += 3)
            {
                int c_eq = outer[row    ];
                int q_eq = outer[row + 1];
                int r_eq = outer[row + 2];
                status =  corner_north_east(values, row, c_eq, q_eq, r_eq, rhs, 
                    theta, nx, ny, htheta, qtheta, rtheta);
            }
            } // end do_linear_wave
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

                    status = bed_shear_stress_matrix_and_rhs(values, row, c_eq, q_eq, r_eq, rhs,
                                htheta, qtheta, rtheta, cf, theta, dx, dy, nx, ny);
                    // boundary_south
                    // boundary_north
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

                    status = convection_matrix_and_rhs(values, row, c_eq, q_eq, r_eq, rhs,
                                htheta, qtheta, rtheta, theta, dx, dy, nx, ny);
                    // boundary_south
                    // boundary_north
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
                // For the moment only interior nodes (2025-08-13)
                // Do not clear the rows, but add the viscosity terms

                // corner_south_west
                // boundary_west
                // corner_north_west

                // interior with south and north boundary
                for (size_t row = 3 * ny; row < 3 * (nx - 1) * ny; row += 3) 
                {
                    int c_eq = outer[row    ];
                    int q_eq = outer[row + 1];
                    int r_eq = outer[row + 2];

                    status = viscosity_matrix_and_rhs(values, row, c_eq, q_eq, r_eq, rhs,
                                htheta, qtheta, rtheta, visc, theta, dx, dy, nx, ny);
                    // boundary_south
                    // boundary_north
                }
                // corner_south_east
                // boundary_east
                // corner_north_east

                STOP_TIMER(Viscosity);
            }

            ////////////////////////////////////////////////////////////////////
            if (linear_solver == "bicgstab")
            {

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
                    for (size_t i = 0; i < 3 * nxny; ++i)
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
                    for (size_t i = 0; i < 3 * nxny; ++i)
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
                    for (size_t i = 0; i < 3 * nxny; ++i)
                    {
                        log_file << std::showpos << std::setprecision(3) << i << "   ";
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
                solver_iterations = solver.iterations();
                if (logging == "iterations" || logging == "matrix")
                {
                    log_file << "time [sec]: " << std::setprecision(4) << std::scientific << time
                        << "    Newton iteration: " << used_newton_iter
                        << "    BiCGstab iterations: " << solver.iterations()
                        << "    estimated error: " << solver.error()
                        << std::endl;
                }
            }
            ////////////////////////////////////////////////////////////////////
            else if (linear_solver == "multigrid")
            {
                if (nst == 1 && iter == 0)
                {
                    START_TIMER(MULTIGRID_initialization);
                }
                else if (nst > 1)
                {
                    START_TIMER(MULTIGRID);
                }
               
                Solver solve(A, prm);
                auto [iters, err] = solve(rhs, solution);
               
                 if (nst == 1 && iter == 0)
                {
                    STOP_TIMER(MULTIGRID_initialization);
                }
                else if (nst > 1)
                {
                    STOP_TIMER(MULTIGRID);
                }
               solver_iterations = iters;
                if (logging == "iterations" || logging == "matrix")
                {
                    log_file << "time [sec]: " << std::setprecision(4) << std::scientific << time
                        << "    Newton iteration: " << used_newton_iter
                        << "    Multigrid iterations: " << iters
                        << "    estimated error: " << err
                        << std::endl;
                }
            }
            ////////////////////////////////////////////////////////////////////

            // 
            // The new solution is the previous iterant plus the delta
            //
            dh_max = 0.0;
            dq_max = 0.0;
            dr_max = 0.0;
            dh_maxi = 0;
            dq_maxi = 0;
            dr_maxi = 0;
            size_t k = 0;
            for (size_t i = 0; i < nx; i++)
            {
                for (size_t j = 0; j < ny; j++)
                {
                    k = main_idx(i, j ,ny);
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
                for (size_t i = 0; i < 3 * nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << solution[i] << std::endl;
                }
                log_file << "=== hp, qp, rp, zeta ==================================" << std::endl;
                for (size_t i = 0; i < nxny; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << hp[i] << ", " << qp[i] << ", " << rp[i] << ", " << hp[i] + zb[i] << std::endl;
                }
                log_file << "=======================================================" << std::endl;
            }
            used_lin_solv_iter = std::max(used_lin_solv_iter, (int) solver_iterations);
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
        //    for (size_t i = 0; i < nr_nodes; ++i)
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
            if (do_convection)
            {
                convection_post_rhs(post_q, post_r, hn, qn, rn, dx, dy, nx, ny);
                map_file->put_time_variable(map_conv_q_name, nst_map, post_q);
                map_file->put_time_variable(map_conv_r_name, nst_map, post_r);
            }
            if (do_bed_shear_stress)
            {
                bed_shear_stress_post_rhs(post_q, post_r, hn, qn, rn, cf, nx, ny);
                map_file->put_time_variable(map_beds_q_name, nst_map, post_q);
                map_file->put_time_variable(map_beds_r_name, nst_map, post_r);
            }
            if (do_viscosity)
            {
                viscosity_post_rhs(post_q, post_r, hn, qn, rn, visc, dx, dy, nx, ny);
                map_file->put_time_variable(map_visc_q_name, nst_map, post_q);
                map_file->put_time_variable(map_visc_r_name, nst_map, post_r);
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

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(hn[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_h_name, nst_his, his_values);

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(qn[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_q_name, nst_his, his_values);

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(rn[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_r_name, nst_his, his_values);

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(s[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_s_name, nst_his, his_values);

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(u[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_u_name, nst_his, his_values);

            his_values.clear();
            for (size_t i = 0; i < input_data.obs_points.size(); ++i)
            {
                his_values.push_back(v[input_data.obs_points[i].idx]);
            }
            his_file->put_variable(his_v_name, nst_his, his_values);

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
inline size_t main_idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
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
