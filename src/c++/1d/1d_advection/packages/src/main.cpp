//
// Programmer: J. Mooiman
// Date      : 2024-12-30
// email     : jan.mooiman@outlook.com
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

int main(int argc, char* argv[])
{
    bool stationary = false;
    std::string toml_file_name("--- not-defined ---");
    int status = -1;
    int nst_print_matrix = -1;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path current_dir;
    std::filesystem::path output_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();

    toml::table tbl;
    toml::table tbl_chp;
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
    }
    else
    {
        std::cout << "No \'toml\' file is read." << std::endl;
        current_dir = ".";
        output_dir = ".";
    }

    bool regularization_init = tbl["regularization_init"].value_or(bool(false));
    bool regularization_iter = tbl["regularization_iter"].value_or(bool(false));
    bool regularization_time = tbl["regularization_time"].value_or(bool(false));
    bool use_eq8 = tbl["use_eq8"].value_or(bool(false));
    if (!use_eq8 && (regularization_init || regularization_iter || regularization_time))
    {
        output_dir = current_dir.string() + "/output_eq8_not_used/";
    }
    std::filesystem::create_directory(output_dir);

    START_TIMERN(main);
    START_TIMER(Simulation initialization);

    std::string out_file;
    std::stringstream ss;
    ss << "advection";
    out_file = output_dir.string() + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Advection, BiCGSTab");

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
    log_file << "=== Used input variables ==============================" << std::endl;

    // Domain
    log_file << "[Domain]" << std::endl;
    tbl_chp = *tbl["Domain"].as_table();
    double Lx = tbl_chp["Lx"].value_or(double(1.0));
    log_file << "Lx = " << Lx << std::endl;
    double Interface = tbl_chp["Interface"].value_or(double(0.5));  // Location of interface[-], fraction of Lx
    log_file << "Interface = " << Interface << std::endl;

    // Time
    log_file << std::endl << "[Time]" << std::endl;
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    log_file << "tstart = " << tstart << std::endl;
    double tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.));
    log_file << "tstop = " << tstop << std::endl;
    double dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary
    if (dt == 0.0) { stationary = true;  }
    log_file << "dt = " << dt << std::endl;

    // Boundary
    log_file << std::endl << "[Boundary]" << std::endl;
    tbl_chp = *tbl["Boundary"].as_table();
    std::vector<double> bc_vals;
    status = get_toml_array(tbl_chp, "bc_vals", bc_vals);
    double treg = tbl_chp["treg"].value_or(double(300.0));
    if (stationary)
    {
        log_file << "treg = " << 0.0 << std::endl;
    }
    else
    {
        log_file << "treg = " << treg << std::endl;
    }

    //Physics
    tbl_chp = *tbl["Physics"].as_table();
    log_file << std::endl << "[Physics]" << std::endl;
    double u_const = tbl_chp["u_const"].value_or(double(0.0));  // default, no velocitty
    log_file << "u_const = " << u_const << std::endl;

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    tbl_chp = *tbl["Numerics"].as_table();
    double dx = tbl_chp["dx"].value_or(double(0.01));
    log_file << "dx = " << dx << std::endl;
    double c_psi = tbl_chp["c_psi"].value_or(double(4.0));
    log_file << "c_psi = " << c_psi << std::endl;
    double theta = tbl_chp["theta"].value_or(double(0.501));
    if (stationary)
    {
        log_file << "theta = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "theta = " << theta << std::endl;
    }
    double alpha = tbl_chp["alpha"].value_or(double(0.125));  // 0.125 (Mass matrix)
    log_file << "alpha = " << alpha << std::endl;
    int iter_max = tbl_chp["iter_max"].value_or(int(150));
    log_file << "iter_max = " << iter_max << std::endl;
    double eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));
    log_file << "eps_newton = " << eps_newton << std::endl;
    double eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));
    log_file << "eps_bicgstab = " << eps_bicgstab << std::endl;
    // already read: bool use_eq8 = tbl_chp["use_eq8"].value_or(bool(false));
    log_file << "use_eq8  = " << use_eq8 << std::endl;
    // bool regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    log_file << "regularization_int  = " << regularization_init << std::endl;
    // bool regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    // bool regularization_time = tbl_chp["regularization_time"].value_or(bool(false));
    log_file << "regularization_time = " << regularization_time << std::endl;

    //Output
    log_file << std::endl << "[Output]" << std::endl;
    tbl_chp = *tbl["Output"].as_table();
    double dt_his = tbl_chp["dt_his"].value_or(double(1.0));
    if (stationary)
    {
        log_file << "dt_his = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_his = " << dt_his << std::endl;
    }
    double dt_map = tbl_chp["dt_map"].value_or(double(10.0));
    if (stationary)
    {
        log_file << "dt_map = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_map = " << dt_map << std::endl;
    }

    //  select 
    //      1 uniform u>0
    //      2 uniform u>0, sine function at west boundary
    int select = 0;
    select = 2;


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

    double bc0;
    double bc1;

    std::vector<double> x(nx, 0.);  // x-coordinate
    std::vector<double> y(nx, 0.);  // y-coordinate
    std::vector<double> cn(nx, 0.0);  // constituent [-]
    std::vector<double> cp(nx, 0.);  // constituent [-]
    std::vector<double> q(nx, 0.);  // source
    std::vector<double> u(nx, u_const);  // advection velocity [m s-1], adv_diff_init_velocity
    std::vector<double> psi(nx, 0.);  //
    std::vector<double> eq8(nx, 0.);  //
    std::vector<double> delta_c(nx, 0.);
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix
    std::vector<double> w_bc(3, 0.);  // weighting coefficients on boundary
    std::vector<double> w_gp(3, 0.);  // weighting coefficients on grid point
    std::vector<double> w_cc(3, 0.);  // weighting coefficients on cell center
    std::vector<double> tmp(nx, 0.);  // dummy vector

    Eigen::VectorXd solution(nx);    // solution vector [c]^{n+1}
    Eigen::VectorXd rhs(nx);        // RHS vector [c]^n

    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    double alpha_bc = (2. * alpha - 1. / 2.);
    w_bc[0] = 1. + alpha_bc;  // 0.75
    w_bc[1] = 1. - 2. * alpha_bc; // 1.5
    w_bc[2] = alpha_bc;  // -0.25

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i-1) * dx;
    }

    // initial concentration
    (void)adv_diff_init_concentration(cn, select);
    (void)adv_diff_init_velocity(u, u_const, x, select);
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
    map_file->put_time_variable(map_psi_name, nst_map, psi);
    map_file->put_time_variable(map_eq8_name, nst_map, eq8);
    // End definition of map file
    ////////////////////////////////////////////////////////////////////////////
    // Define time history file
    std::cout << "    Create his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations (corrected for virtual points)
    int i_left = 1;
    int i_mid_left = nx / 4 + 1;
    int i_mid = nx / 2;
    int i_mid_right = 3 * nx / 4 - 1;
    int i_right = nx - 2;
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

    std::string his_newton_iter_name("newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "", "Newton iterations", "-");
    std::string his_lin_solv_iter_name("lin_solver_iterations");
    his_file->add_variable_without_location(his_lin_solv_iter_name, "", "Lin. solver iterations", "-");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, double(0)* dt);

    std::vector<double> his_values = { cn[i_left], cn[i_mid_left], cn[i_mid], cn[i_mid_right], cn[i_right] };
    his_file->put_variable(his_cn_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);

    ////////////////////////////////////////////////////////////////////////////

    // source term q
    (void)adv_diff_source(q, select);

    double time = tstart + dt * double(0);
    //(void)adv_diff_boundary_condition(cn[0], cn[nx - 1], time, treg, select);
    // compute Peclet number
    for (int i = 0; i < nx; i++)
    {
        cp[i] = cn[i];  // place the old solution in the new one, needed for the iteration procedure
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
        START_TIMER(Newton iteration);
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
            for (int i = 1; i < nx - 1; ++i)
            {
                double cn_im12 = 0.5 * (cn[i - 1] + cn[i]);
                double cn_ip12 = 0.5 * (cn[i] + cn[i + 1]);

                double cp_im12 = 0.5 * (cp[i - 1] + cp[i]);
                double cp_ip12 = 0.5 * (cp[i] + cp[i + 1]);

                double ctheta_im12 = theta * cp_im12 + (1.0 - theta) * cn_im12;
                double ctheta_ip12 = theta * cp_ip12 + (1.0 - theta) * cn_ip12;

                double u_im12 = 0.5 * (u[i - 1] + u[i]);
                double u_ip12 = 0.5 * (u[i] + u[i + 1]);

                A.coeffRef(i, i - 1) = dx * dtinv * mass[0]
                    - 0.5 * u_im12 * theta;
                A.coeffRef(i, i) = dtinv_pseu + dx * dtinv * mass[1]
                    - 0.5 * u_im12 * theta + 0.5 * u_ip12 * theta;
                A.coeffRef(i, i + 1) = dx * dtinv * mass[2]
                    + 0.5 * u_ip12 * theta;
                rhs[i] = -(
                      dx * dtinv * mass[0] * (cp[i - 1] - cn[i - 1])
                    + dx * dtinv * mass[1] * (cp[i] - cn[i])
                    + dx * dtinv * mass[2] * (cp[i + 1] - cn[i + 1])
                    + u_ip12 * ctheta_ip12 - u_im12 * ctheta_im12
                    );
            }
            {
                //
                // wwest boundary
                //
                int i = 0;

                double cp_i = cp[i];       // = c^{n+1,p}_{i}
                double cp_ip1 = cp[i + 1];       // = c^{n+1,p}_{i+1}
                double cp_ip2 = cp[i + 2];       // = c^{n+1,p}_{i+2}

                // Dirichlet
                // double bc_location_ess = 1.0;
                // w_bc = bc_val(bc_location_ess);
                w_bc[0] = 1. / 12.;
                w_bc[1] = 10. / 12.;
                w_bc[2] = 1. / 12.;
                w_bc[0] = 0.5 * (1. + alpha_bc);
                w_bc[1] = 0.5 * (1. - 2. * alpha_bc);
                w_bc[2] = 0.5 * (alpha_bc);
                A.coeffRef(i, i    ) = w_bc[0];
                A.coeffRef(i, i + 1) = w_bc[1];
                A.coeffRef(i, i + 2) = w_bc[2];
                rhs[i] = +bc0 - (w_bc[0] * cp_i + w_bc[1] * cp_ip1 + w_bc[2] * cp_ip2);  // if u>0 this is upwind
            }
            {
                //
                // eeast boundary
                //
                int i = nx - 1;
                // bc_location_nat = 0.5;  // boundary location between node i+1 and i
                // w_bc = bc_val(bc_location_ess);
                w_bc[0] = 11. / 24.;
                w_bc[1] = 14. / 24.;
                w_bc[2] = -1. / 24.;
                w_bc[0] = 1. + alpha_bc;
                w_bc[1] = 1. - 2. * alpha_bc;
                w_bc[2] = alpha_bc;

                double cn_i = cn[i];       // = h^{n}_{i}
                double cn_im1 = cn[i - 1];       // = h^{n}_{i-1}
                double cn_im2 = cn[i - 2];       // = h^{n}_{i-2}
                double cp_i = cp[i];       // = h^{n+1,p}_{i}
                double cp_im1 = cp[i - 1];       // = h^{n+1,p}_{i-1}
                double cp_im2 = cp[i - 2];       // = h^{n+1,p}_{i-2}
                double ctheta_i = theta * cp_i + (1.0 - theta) * cn_i;
                double ctheta_im1 = theta * cp_im1 + (1.0 - theta) * cn_im1;

                // Outflow boundary (natural boundary)
                double dcdt = dtinv * (
                    w_bc[0] * (cp_i - cn_i) +
                    w_bc[1] * (cp_im1 - cn_im1) +
                    w_bc[2] * (cp_im2 - cn_im2)
                    );

                double u_im12 = 0.5 * (u[i] + u[i - 1]);
                double udcdx =  u_im12 * (ctheta_i - ctheta_im1) * dxinv;

                A.coeffRef(i, i) = dtinv * w_bc[0] + u_im12 * theta * dxinv;
                A.coeffRef(i, i - 1) = dtinv * w_bc[1] - u_im12 * theta * dxinv;
                A.coeffRef(i, i - 2) = dtinv * w_bc[2];
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
            for (int i = 0; i < nx; ++i)
            {
                cp[i] += solution[i];
                delta_c[i] = solution[i]; 
                if (dc_max < std::abs(delta_c[i]))
                {
                    dc_max = std::abs(delta_c[i]);
                    dc_maxi = i;
                }
            }

            if (nst <= nst_print_matrix)
            {
                log_file << "=== Matrix ============================================" << std::endl;
                log_file << std::setprecision(4) << std::scientific << Eigen::MatrixXd(A) << std::endl;
                log_file << "=== RHS ===============================================" << std::endl;
                log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
            }
            if (nst <= nst_print_matrix)
            {
                log_file << "=== cn, delta_c, cp ===================================" << std::endl;
                for (int i = 0; i < nx; ++i)
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
        }
        STOP_TIMER(Newton iteration);
        for (int i = 0; i < nx; ++i)
        {
            cn[i] = cp[i];
        }
        if (nst <= nst_print_matrix)
        {
            log_file << "=== cp, cn, cp-cn =====================================" << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                log_file << std::setprecision(15) << std::scientific << cp[i] << ",   " << cn[i] << ",   " << cp[i] - cn[i] << std::endl;
            }
        }
        std::cout << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << ";   " << used_iter << std::endl;
        if (nst <= nst_print_matrix)
        {
            log_file << "=======================================================" << std::endl;
            log_file << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << ";   " << used_iter << std::endl;
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
            map_file->put_time_variable(map_psi_name, nst_map, psi);
            map_file->put_time_variable(map_eq8_name, nst_map, eq8);
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
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
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
        int nx = x.size();
        x_gp[0] = x[nx - 1];
        x_gp[1] = x[nx - 2];
        x_gp[2] = x[nx - 3];
        f_gp[0] = f[0];
        f_gp[1] = f[1];
        f_gp[2] = f[2];
        x_loc = (1.0 - x_loc) * x[nx - 1] + x_loc * x[nx - 2];
    }
    double b_1 = (f_gp[1] - f_gp[0]) / (x_gp[1] - x_gp[0]);
    double b_2 = ((f_gp[2] - f_gp[1]) / (x_gp[2] - x_gp[1]) - (f_gp[1] - f_gp[0]) / (x_gp[1] - x_gp[0])) / (x_gp[2] - x_gp[0]);
    double f_xloc = f_gp[0] + b_1 * (x_loc - x_gp[0]) + b_2 * (x_loc - x_gp[0]) * (x_loc - x_gp[1]);
    double dfdx_xloc = b_1 + b_2 * (x_loc - x_gp[0]) + b_2 * (x_loc - x_gp[1]);
    double d2fdx2_xloc = 2 * b_2;

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
        int nx = x.size();
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
    double p = h00 * f_cc[0] + h10 * (x_cc[1] - x_cc[0]) * dfdx_cc[0] + h01 * f_cc[1] + h11 * (x_cc[1] - x_cc[0]) * dfdx_cc[1];

    // Calculate weights such that p = w[0] * f[0] + w[1] * f[1] + w[0] * f[0] in non - equidistant grids
    // f_cc = 0.5 * (f_gp[0] + f_gp[1]), f_cc[1] = 0.5 * (f_gp[1] + f_gp[0])
    std::vector<double> w(3, 0.0);
    w[0] = h00 * 0.5 - h10 * (x_cc[1] - x_cc[0]) / (x_gp[1] - x_gp[0]);
    w[1] = h00 * 0.5 + h10 * (x_cc[1] - x_cc[0]) / (x_gp[1] - x_gp[0]) + h01 * 0.5 - h11 * (x_cc[1] - x_cc[0]) / (x_gp[0] - x_gp[1]);
    w[2] = h01 * 0.5 + h11 * (x_cc[1] - x_cc[0]) / (x_gp[0] - x_gp[1]);

    return w;
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

