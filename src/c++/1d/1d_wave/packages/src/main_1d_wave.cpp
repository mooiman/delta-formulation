//
// programmer: Jan Mooiman
// Date      : 2025-01-03
// Email     : jan.mooiman@outlook.com
//

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>

// for bicgstab solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <toml.h>


#include "cfts.h"
#include "ugrid1d.h"
#include "perf_timer.h"
#include "regularization.h"
#include "definition_map_file.h"
#include "boundary_condition.h"

enum class BED_LEVEL
{
    NONE = 0,
    FLAT,
    SHOAL,
    SLOPED,
    WAVY,
    WAVY_SLOPED,
    WEIR,
    NR_BED_LEVELS
};

int p_index(int i, int j, int nx);
int read_bed_level(std::string filename, std::vector<double> & value);
int initialize_bed_level(BED_LEVEL, std::vector<double>&, std::vector<double>&, std::string&, double);
int initialize_scalar(double, std::vector<double>&, std::vector<double>&);
double scv(double, double);
double Fabs(double, double);
void GetArguments(long argc, char** argv, std::string *file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string> &);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);


// Solve the linear wave equation
// Continuity equation: d(h)/dt + d(q)/dx = 0
// Momentum equation: d(q)/dt + g h d(zeta)/dx = 0

int main(int argc, char* argv[])
{
    std::map<BED_LEVEL, std::string> bed_level_name;
    bed_level_name[BED_LEVEL::FLAT] = "flat";
    bed_level_name[BED_LEVEL::SHOAL] = "shoal";
    bed_level_name[BED_LEVEL::SLOPED] = "sloped";
    bed_level_name[BED_LEVEL::WAVY] = "wavy";
    bed_level_name[BED_LEVEL::WAVY_SLOPED] = "wavy_sloped";
    bed_level_name[BED_LEVEL::WEIR] = "weir";

    int BC_WEST = 0;
    int BC_EAST = 1;

    int status = -1;
    std::string toml_file_name;

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

    START_TIMERN(Main);
    START_TIMER(Writing log-file);
    BED_LEVEL bed_level_type = BED_LEVEL::NONE;

    std::string out_file;
    std::stringstream ss;

    std::string logging = tbl["Logging"].value_or("None");

    bool stationary = false;
    double sign = 1.0;
    double dt = tbl["Numerics"]["dt"].value_or(double(0.0));  // default stationary simulation
    if (dt == 0.0)
    {
        stationary = true;
    }
    double dx = tbl["Numerics"]["dx"].value_or(double(10.));
    double depth = tbl["Domain"]["depth"].value_or(double(10.));

    std::string geometry_type = tbl["Domain"]["geometry_type"].value_or("None");
    if (geometry_type == "WavyBedGeometry") { bed_level_type = BED_LEVEL::WAVY; }
    else if (geometry_type == "SlopedBedGeometry") { bed_level_type = BED_LEVEL::SLOPED; }
    else if (geometry_type == "SlopedWavyBedGeometry") { bed_level_type = BED_LEVEL::WAVY_SLOPED; }
    else if (geometry_type == "UniformGeometry") {
        bed_level_type = BED_LEVEL::FLAT;
        std::stringstream depth_strm;
        depth_strm << std::fixed << std::setprecision(0) << depth;
        bed_level_name[BED_LEVEL::FLAT] = "flat" + depth_strm.str();
    }
    else if (geometry_type == "WeirBorsboom2000Geometry") { bed_level_type = BED_LEVEL::WEIR; }
    else if (geometry_type == "Shoal") { bed_level_type = BED_LEVEL::SHOAL; }
    else { bed_level_type = BED_LEVEL::NONE; }

    ss << "_dx" << dx << "_dt" << dt;

    std::string inp_bed = "bed_level_" + bed_level_name[bed_level_type];

    out_file = output_dir.string() + inp_bed + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");

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
    double Lx = tbl["Domain"]["Lx"].value_or(double(12000.));

    std::string model_title;

    // Time
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    double tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.));

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    std::vector<std::string> ini_vars;  // get the element as an array
    status = get_toml_array(tbl_chp, "ini_vars", ini_vars);
    double gauss_mu = tbl_chp["gauss_mu"].value_or(double(0.0));  // location of th gaussian hump
    double gauss_sigma = tbl_chp["gauss_sigma"].value_or(double(0.25 * Lx));  // sigma of the Guassian hump
    double a0 = tbl_chp["a0"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary

    // Boundary
    tbl_chp = *tbl["Boundary"].as_table();
    double eps_bc_corr = tbl_chp["eps_bc_corr"].value_or(double(0.0001));  // default 1e-4
    std::vector<std::string> bc_type;
    status = get_toml_array(tbl_chp, "bc_type", bc_type);
    double treg = tbl_chp["treg"].value_or(double(150.));  // Regularization time to reach constant boundary value
    std::vector<std::string> bc_vars;  // get the element as an array
    status = get_toml_array(tbl_chp, "bc_vars", bc_vars);
    std::vector<double> bc_vals;
    status = get_toml_array(tbl_chp, "bc_vals", bc_vals);
    std::vector<bool> bc_absorbing;
    status = get_toml_array(tbl_chp, "bc_absorbing", bc_absorbing);

    // Physics
    tbl_chp = *tbl["Physics"].as_table();
    double g = tbl_chp["g"].value_or(double(10.));  // Gravitational acceleration
    bool momentum_convection = tbl_chp["momentum_convection"].value_or(bool(false));
    bool momentum_bed_shear_stress = tbl_chp["momentum_bed_shear_stress"].value_or(bool(false));
    bool momentum_viscosity = tbl_chp["momentum_viscosity"].value_or(bool(false));
    double chezy_coefficient = tbl_chp["chezy_coefficient"].value_or(double(50.0));
    double visc_const = tbl_chp["viscosity"].value_or(double(0.0));

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    // dt already read
    // dx already read
    double c_psi = tbl_chp["c_psi"].value_or(double(4.));
    int iter_max = tbl_chp["iter_max"].value_or(int(50));
    double theta = tbl_chp["theta"].value_or(double(0.501));  // Implicitness factor (0.5 <= theta <= 1.0)
    double alpha = tbl_chp["alpha"].value_or(double(1./8.));  // Linear (spatial) interpolation coefficient
    double eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));  // stop criterium for Newton iteration
    double eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));  // stop criterium for BiCGStab iteratio
    double eps_fabs = tbl_chp["eps_fabs"].value_or(double(1.0e-2));  // epsilon needed to approximate the abs-function by a continues function
    bool regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    bool regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    bool regularization_time = tbl_chp["regularization_time"].value_or(bool(false));
    bool use_eq8 = tbl_chp["use_eq8"].value_or(bool(false));
    if (!use_eq8 && (regularization_init || regularization_iter || regularization_time))
    {
        output_dir = current_dir.string() + "/output_eq8_not_used/";
    }

    // Output
    tbl_chp = *tbl["Output"].as_table();
    double dt_his = tbl_chp["dt_his"].value_or(double(1.0));  // write interval to his-file
    double dt_map = tbl_chp["dt_map"].value_or(double(0.0));  // write interval to his-file

    REGULARIZATION* regularization = new REGULARIZATION(iter_max);

    double dxinv = 1./dx;                                 // invers grid size [m]
    int nx = int(Lx * dxinv) + 1 + 2 ;                    // number of nodes; including 2 virtual points

    int total_time_steps = int((tstop - tstart) / dt) + 1;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    int wrihis;                                           // write interval to his-file
    int wrimap;                                           // write interval to map-file
    if (stationary)
    {
        dt = 0.0;                                         // Time step size [s]
        dtinv = 0.0;                                      // stationary solution
        eps_bc_corr = 1.0;
        iter_max = 2.5 * iter_max;
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
    log_file << "=== Used input variables ==============================" << std::endl;
    log_file << "[Domain]" << std::endl;
    log_file << "Lx = " << Lx << std::endl;
    log_file << "depth = " << depth << std::endl;
    log_file << "geometry_type = " << geometry_type << std::endl;

    log_file << std::endl << "[Time]" << std::endl;
    log_file << "tstart = " << tstart << std::endl;
    log_file << "tstop = " << tstop << std::endl;

    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "gauss_mu = " << gauss_mu << std::endl;
    log_file << "gauss_sigma = " << gauss_sigma << std::endl;
    log_file << "a0 = " << a0 << std::endl;

    log_file << std::endl << "[Boundary]" << std::endl;
    log_file << "bc_type = ";
    for (int i = 0; i < bc_type.size() - 1; ++i) { log_file << bc_type[i] << ", ";}
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
    log_file << "g = " << g << std::endl;
    log_file << "momentum_convection = " << momentum_convection << std::endl;
    log_file << "momentum_bed_shear_stress = " << momentum_bed_shear_stress << std::endl;
    log_file << "momentum_viscosity = " << momentum_viscosity << std::endl;
    log_file << "chezy_coefficient = " << chezy_coefficient << std::endl;
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
    log_file << "eps_abs_function  = " << eps_fabs << std::endl;
    log_file << "use_eq8  = " << use_eq8 << std::endl;
    log_file << "regularization_init = " << regularization_init << std::endl;
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    log_file << "regularization_time = " << regularization_time << std::endl;
    
    log_file << std::endl << "[Output]" << std::endl;
    log_file << "dt_his = " << dt_his << std::endl;
    log_file << "dt_map = " << dt_map << std::endl;

    std::vector<double> x(nx, 0.);                      // x-coordinate
    std::vector<double> y(nx, 0.);                      // y-coordinate
    std::vector<double> zb(nx, -10.0);                  // regularized bed level
    std::vector<double> zb_ini(nx, -10.0);              // initial bed level

    std::vector<double> rhs_viscosity(2 * nx, 0.);      // rhs of the momentum equation for viscosity
    std::vector<double> tmp1(nx, 0.);                 // help array
    std::vector<double> tmp2(nx, 0.);                 // help array

    std::vector<double> qn(nx, 0.);                     // flow at (n)
    std::vector<double> hn(nx, 0.);                     // total depth at (n)
    std::vector<double> qp(nx, 0.);                     // flow at (n+1,p), previous iteration
    std::vector<double> hp(nx, 0.);                     // total depth at (n+1,p), previous iteration
    std::vector<double> delta_h(nx, 0.);
    std::vector<double> delta_q(nx, 0.);
    std::vector<double> s_giv(nx, 0.);  // given initial water level
    std::vector<double> u_giv(nx, 0.);  // given initial velocity
    std::vector<double> s(nx, 0.);  // smoothed water level
    std::vector<double> u(nx, 0.);  // water velocity, needed for postprocessing
    std::vector<double> riemann_pos(nx, 0.);  // Riemann invarinat going to the right, needed for postprocessing
    std::vector<double> riemann_neg(nx, 0.);  // Riemann invarinat going to the left, needed for postprocessing
    std::vector<double> visc_given(nx, visc_const);  // Initialize viscosity array with given value
    std::vector<double> visc_reg(nx, visc_const);  // Initialize given viscosity array with regularized value
    std::vector<double> visc(nx, visc_const);  // Viscosity array used for computation, adjusted for cell peclet number
    std::vector<double> cf(nx, 0.);  // Bed sheart stress coefficient
    std::vector<double> pe(nx, 0.);  // peclet number [-]
    std::vector<double> psi(nx, 0.);  //
    std::vector<double> eq8(nx, 0.);  //
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix

    Eigen::VectorXd solution(2 * nx);               // solution vector [h, q]^{n}
    Eigen::VectorXd rhs(2 * nx);                // RHS vector [h, q]^n

    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;
    
    double alpha_bc = 2. * alpha - 1. / 2.;
    std::vector<double> w_nat(3, 0.0);
    w_nat[0] = 0.5 * (1.0 + alpha_bc);
    w_nat[1] = 0.5 * (1.0 - 2.0 * alpha_bc);
    w_nat[2] = 0.5 * alpha_bc;
    std::vector<double> w_ess(3, 0.0);
    w_ess[0] = 1. / 12.;
    w_ess[1] = 10. / 12.;
    w_ess[2] = 1. / 12.;

    //initialize water level
    std::cout << "    Initialisation" << std::endl;

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i - 1) * dx - Lx / 2;
    }
    status = 0;
    double s_offset = 0.0;
    double zb_slope = 0.0;
    for (int i = 0; i < nx; i++)
    {
        s_giv[i] = s_offset;
        if (ini_vars[0] == "zeta")
        {
            s_giv[i] = 2.0 * a0 * std::exp( -(x[i] - gauss_mu) * (x[i] - gauss_mu) / (2.* gauss_sigma * gauss_sigma) ) + s_offset;  // Initial water level, A.Slob 1983 uses 200. instead of 500.
        }
        s[i] = s_giv[i];
        u_giv[i] = 0.0;
        if (ini_vars[0] == "viscosity")
        {
            // heat equation will be used to check the implementation of the viscosity term
            double nu = visc_const;
            double pi = M_PI;
            double t0 = 10.0;
            u_giv[i] = 1.0 / std::sqrt(4. * pi * nu * t0) * std::exp(-x[i] * x[i] / (4.0 * nu * t0));
        }
        u[i] = u_giv[i];
        // Bed shear coefficient
        cf[i] = g / (chezy_coefficient * chezy_coefficient);
    }
    status = initialize_bed_level(bed_level_type, x, zb_ini, model_title, depth);
    if (regularization_init)
    {
        (void)regularization->given_function(zb, psi, eq8, zb_ini, dx, c_psi, use_eq8);
    }
    else
    {
        for (int i = 0; i < zb_ini.size(); ++i)
        {
            zb[i] = zb_ini[i];
        }
    }
    double min_zb = *std::min_element(zb.begin(), zb.end());
    for (int i = 0; i < nx; i++)
    {
        hn[i] = s[i] - zb[i];  // Initial water depth
        hp[i] = hn[i];  // initial water depth
        qn[i] = hn[i] * u[i];  // initial flow flux
        qp[i] = qn[i];  // initial flow flux
        riemann_pos[i] = std::abs(qn[i] + std::sqrt(g * hn[i]) * hn[i]);
        riemann_neg[i] = std::abs(qn[i] - std::sqrt(g * hn[i]) * hn[i]);
    }

    log_file << "-------------------------------------------------------" << std::endl;
    log_file << "Nodes  : " << nx << std::endl;
    log_file << "Volumes: " << (nx-1) << std::endl;
    log_file << "CFL    : " << std::sqrt(g * std::abs(min_zb)) * dt / dx << std::endl;
    STOP_TIMER(Writing log-file);  // but two write statements are not timed
    if (status != 0) {
        log_file.close();
        exit(1);
    }
    log_file << "=======================================================" << std::endl;


    //status = initialize_scalar(alpha, s_giv, s);

    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        (void)regularization->given_function(visc_reg, psi, eq8, visc_given, dx, c_psi, use_eq8);
        STOP_TIMER(Regularization_init);
    }

    // Merge delta depth and delta flux in one solution-vector (initial equal to zero)
    for (int i = 0; i < nx; i++)
    {
        solution[2 * i] = 0.0;  // Delta h; continuity
        solution[2 * i + 1] = 0.0;  // Delta q; momentum
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
    map_names.push_back("eq8_1d");
    map_names.push_back("pe_id");

    map_file = create_map_file(nc_mapfilename, model_title, x, map_names);

    // Put data on map file
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    map_file->put_time_variable(map_names[0], nst_map, hn);
    map_file->put_time_variable(map_names[1], nst_map, qn);
    map_file->put_time_variable(map_names[2], nst_map, s);
    map_file->put_time_variable(map_names[3], nst_map, u);
    map_file->put_time_variable(map_names[4], nst_map, zb);
    map_file->put_time_variable(map_names[5], nst_map, delta_h);
    map_file->put_time_variable(map_names[6], nst_map, delta_q);
    map_file->put_time_variable(map_names[7], nst_map, visc_reg);
    map_file->put_time_variable(map_names[8], nst_map, visc);
    map_file->put_time_variable(map_names[9], nst_map, psi);
    map_file->put_time_variable(map_names[10], nst_map, eq8);
    map_file->put_time_variable(map_names[11], nst_map, pe);

    ////////////////////////////////////////////////////////////////////////////
    // Create time history file
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

    std::string his_h_name("Water_depth");
    std::string his_q_name("Water_flux");
    std::string his_s_name("Water_level");
    std::string his_u_name("Water_velocity");
    std::string his_zb_name("Bed_level");
    std::string his_riemann_pos_name("Riemann_right_going");
    std::string his_riemann_neg_name("Riemann_left_going");
    his_file->add_variable(his_h_name, "", "Water depth", "m");
    his_file->add_variable(his_q_name, "", "Water flux", "m2 s-1");
    his_file->add_variable(his_s_name, "", "Water level", "m");
    his_file->add_variable(his_u_name, "", "Water velocity", "m s-1");
    his_file->add_variable(his_zb_name, "", "Bed level", "m");
    //his_file->add_variable(his_riemann_pos_name, "", "Rieman(+)", "m2 s-1");
    //his_file->add_variable(his_riemann_neg_name, "", "Rieman(-)", "m2 s-1");

    std::string his_newton_iter_name("newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "", "Newton iterations", "-");
    std::string his_lin_solv_iter_name("lin_solver_iterations");
    his_file->add_variable_without_location(his_lin_solv_iter_name, "", "Lin. solver iterations", "-");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, time);

    std::vector<double> his_values = { hn[i_left], hn[i_mid_left], hn[i_mid], hn[i_mid_right], hn[i_right] };
    his_file->put_variable(his_h_name, nst_his, his_values);
    his_values = { qn[i_left], qn[i_mid_left], qn[i_mid], qn[i_mid_right], qn[i_right] };
    his_file->put_variable(his_q_name, nst_his, his_values);
    his_values = { s[i_left], s[i_mid_left], s[i_mid], s[i_mid_right], s[i_right] };
    his_file->put_variable(his_s_name, nst_his, his_values);
    his_values = { u[i_left], u[i_mid_left], u[i_mid], u[i_mid_right], u[i_right] };
    his_file->put_variable(his_u_name, nst_his, his_values);
    his_values = { zb[i_left], zb[i_mid_left], zb[i_mid], zb[i_mid_right], zb[i_right] };
    his_file->put_variable(his_zb_name, nst_his, his_values);

    //his_values = { riemann_pos[i_left], riemann_pos[i_mid_left], riemann_pos[i_mid], riemann_pos[i_mid_right], riemann_pos[i_right] };
    //his_file->put_variable(his_riemann_pos_name, nst_his, his_values);
    //his_values = { riemann_neg[i_left], riemann_neg[i_mid_left], riemann_neg[i_mid], riemann_neg[i_mid_right], riemann_neg[i_right] };
    //his_file->put_variable(his_riemann_neg_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    his_values = { 0 };
    his_file->put_variable(his_lin_solv_iter_name, nst_his, his_values);

    ////////////////////////////////////////////////////////////////////////////

    Eigen::SparseMatrix<double> A(2 * nx, 2 * nx);
    for (int i = 0; i < 2 * nx; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
    }

    std::cout << "Start time-loop" << std::endl;
    std::cout << std::fixed << std::setprecision(2) << "tstart= " << tstart + time << ";   tstop=" << tstart + tstop << ";   dt= " << dt << std::endl;
 
    STOP_TIMER(Initialization);
   // Start time loop
    double dh_max = 0.0;
    int dh_maxi = 0;
    double dq_max = 0.0;
    int dq_maxi = 0;
    START_TIMER(Time loop);
    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        time = dt * double(nst);
        // Set the right-hand side (rhs) of the equations (at t = nst*dt)
        //std::cout << std::setprecision(5) << "Time: " << dt*double(nst) << std::endl;

        double hn_im2, hn_im1, hn_i, hn_ip1, hn_ip2;  // previous timestep (n)
        double hn_im12, hn_ip12;
        double hp_im2, hp_im1, hp_i, hp_ip1, hp_ip2;  // previous iteration (p)
        double hp_im12, hp_ip12;

        double qn_im2, qn_im1, qn_i, qn_ip1, qn_ip2;  // previous timestep (n)
        double qn_im12, qn_ip12;
        double qp_im2, qp_im1, qp_i, qp_ip1, qp_ip2;  // previous iteration (p)
        double qp_im12, qp_ip12;

        double htheta_i;
        double htheta_im1, htheta_im2, htheta_ip1, htheta_ip2;
        double htheta_im12, htheta_ip12;
        double qtheta_i;
        double qtheta_im1, qtheta_im2, qtheta_ip1, qtheta_ip2;
        double qtheta_im12, qtheta_ip12;


        double bc0;
        double bc1;
        double wh_bnd = 0.0;  // west depth-boundary
        double wq_bnd = 0.0;  // west discharge discharge
        double wz_bnd = 0.0;  // west zeta-boundary
        double wu_bnd = 0.0;  // west u-boundary
        double eh_bnd = 0.0;  // east depth-boundary
        double eq_bnd = 0.0;  // east discharge boundary
        double ez_bnd = 0.0;  // east zeta-boundary
        double eu_bnd = 0.0;  // east u-boundary
        int select = 1;  // given value at both sides
        (void) boundary_condition(bc0, bc1, bc_vals[BC_WEST], bc_vals[BC_EAST], time, treg, select);
        if (bc_vars[BC_WEST] == "h")
        {
            wh_bnd = bc0;
        }
        if (bc_vars[BC_WEST] == "q")
        {
            wq_bnd = bc0;
        }
        if (bc_vars[BC_WEST] == "zeta")
        {
            wz_bnd = bc0;
        }
        if (bc_vars[BC_WEST] == "u")
        {
            wu_bnd = bc0;
        }

        if (bc_vars[BC_EAST] == "h")
        {
            eh_bnd = bc1;
        }
        if (bc_vars[BC_EAST] == "q")
        {
            eq_bnd = bc1;
        }
        if (bc_vars[BC_EAST] == "zeta")
        {
            ez_bnd = bc1;
        }
        if (bc_vars[BC_EAST] == "u")
        {
            eu_bnd = bc1;
        }

        int used_newton_iter = 0;
        int used_lin_solv_iter = 0;
        START_TIMER(Newton iteration);
        for (int iter = 0; iter < iter_max; ++iter)
        {
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
                hn_i = hn[i];       // = h^{n}_{i}
                hn_ip1 = hn[i + 1]; // = h^{n}_{i+1}

                qn_im1 = qn[i - 1]; // = q^{n}_{i-1}
                qn_i = qn[i];       // = q^{n}_{i}
                qn_ip1 = qn[i + 1]; // = q^{n}_{i+1}

                hn_im12 = 0.5 * (hn_i + hn_im1);
                hn_ip12 = 0.5 * (hn_ip1 + hn_i);

                qn_im12 = 0.5 * (qn_i + qn_im1);
                qn_ip12 = 0.5 * (qn_ip1 + qn_i);

                hp_im1 = hp[i - 1]; // = h^{n+1,p}_{i-1}
                hp_i = hp[i];       // = h^{n+1,p}_{i}
                hp_ip1 = hp[i + 1]; // = h^{n+1,p}_{i+1}

                qp_im1 = qp[i - 1]; // = q^{n+1,p}_{i-1}
                qp_i = qp[i];  // = q^{n+1,p}_{i}
                qp_ip1 = qp[i + 1]; // = q^{n+1,p}_{i+1}

                hp_im12 = 0.5 * (hp_i + hp_im1);  // = h^{n+1,p}_{i-1/2}
                hp_ip12 = 0.5 * (hp_ip1 + hp_i);  // = h^{n+1,p}_{i+1/2}

                qp_im12 = 0.5 * (qp_i + qp_im1);  // = q^{n+1,p}_{i-1/2}
                qp_ip12 = 0.5 * (qp_ip1 + qp_i);  // = q^{n+1,p}_{i+1/2}

                htheta_i = theta * hp_i + (1.0 - theta) * hn_i;
                htheta_im1 = theta * hp_im1 + (1.0 - theta) * hn_im1;
                htheta_ip1 = theta * hp_ip1 + (1.0 - theta) * hn_ip1;
                htheta_im12 = theta * hp_im12 + (1.0 - theta) * hn_im12;
                htheta_ip12 = theta * hp_ip12 + (1.0 - theta) * hn_ip12;

                qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                qtheta_im1 = theta * qp_im1 + (1.0 - theta) * qn_im1;
                qtheta_ip1 = theta * qp_ip1 + (1.0 - theta) * qn_ip1;
                qtheta_im12 = theta * qp_im12 + (1.0 - theta) * qn_im12;
                qtheta_ip12 = theta * qp_ip12 + (1.0 - theta) * qn_ip12;

                int ph = 2 * p_index(i, 0, nx);  // continuity equation
                int ph_e = 2 * p_index(i + 1, 0, nx);  // continuity equation
                int ph_w = 2 * p_index(i - 1, 0, nx);  // continuity equation
                int pq = ph + 1;  // q=hu-momentum equation
                int pq_e = ph_e + 1;  // q=hu-momentum equation
                int pq_w = ph_w + 1;  // q=hu-momentum equation
                //
                // continuity equation (dh/dt ... = 0)
                //
                rhs[ph] = 0.0;
                // mass-matrix continuity
                A.coeffRef(ph, ph_w) = dtinv * dx * mass[0];
                A.coeffRef(ph, ph) = dtinv * dx * mass[1];
                A.coeffRef(ph, ph_e) = dtinv * dx * mass[2];
                // flow flux
                A.coeffRef(ph, ph_w + 1) = -0.5 * theta;
                A.coeffRef(ph, ph + 1) = 0.5 * theta - 0.5 * theta;
                A.coeffRef(ph, ph_e + 1) = 0.5 * theta;
                //
                rhs[ph] -=
                    dtinv * dx * mass[0] * (hp_im1 - hn_im1) +
                    dtinv * dx * mass[1] * (hp_i - hn_i) +
                    dtinv * dx * mass[2] * (hp_ip1 - hn_ip1);
                // flux
                rhs[ph] += -(qtheta_ip12 - qtheta_im12);
                //
                // momentum (dq/dt ... = 0)
                //
                rhs[pq] = 0.0;
                // mass-matrix momentum
                A.coeffRef(pq, ph_w + 1) = dtinv * dx * mass[0];
                A.coeffRef(pq, ph + 1) = dtinv * dx * mass[1];
                A.coeffRef(pq, ph_e + 1) = dtinv * dx * mass[2];
                // 
                rhs[pq] -=
                    dtinv * dx * mass[0] * (qp_im1 - qn_im1) +
                    dtinv * dx * mass[1] * (qp_i - qn_i) +
                    dtinv * dx * mass[2] * (qp_ip1 - qn_ip1);
                //
                // pressure
                //
                double ztheta_im1 = htheta_im1 + zb[i - 1];
                double ztheta_i = htheta_i + zb[i];
                double ztheta_ip1 = htheta_ip1 + zb[i + 1];
                // Contribution due to Delta h^{n+1, p+1}
                A.coeffRef(pq, ph_w) = 0.125 * theta * g * (ztheta_i - ztheta_im1);
                A.coeffRef(pq, ph) = 0.375 * theta * (ztheta_i - ztheta_im1) + 0.375 * theta * g * (ztheta_ip1 - ztheta_i);
                A.coeffRef(pq, ph_e) = 0.125 * theta * g * (ztheta_ip1 - ztheta_i);
                // contribution due to Delta zeta^{n+1, p+1} 
                A.coeffRef(pq, ph_w) += -0.125 * theta * g * (htheta_im1 + 3. * htheta_i);
                A.coeffRef(pq, ph) += 0.125 * theta * g * (htheta_im1 + 3. * htheta_i) - 0.125 * theta * g * (htheta_ip1 + 3. * htheta_i);
                A.coeffRef(pq, ph_e) += 0.125 * theta * g * (htheta_ip1 + 3. * htheta_i);

                double rhs_pressure = -0.5 * g * 0.25 * (
                    (1.0 * htheta_im1 + 3.0 * htheta_i) * (ztheta_i - ztheta_im1) +
                    (3.0 * htheta_i + 1.0 * htheta_ip1) * (ztheta_ip1 - ztheta_i)
                    );
                rhs[pq] += rhs_pressure;
                //
                if (momentum_convection)
                {
                    //
                    // convection
                    //
                    A.coeffRef(pq, ph_w) += +0.5 * theta * qtheta_im12 * qtheta_im12 / (htheta_im12 * htheta_im12);
                    A.coeffRef(pq, ph) += -0.5 * theta * qtheta_ip12 * qtheta_ip12 / (htheta_ip12 * htheta_ip12) + 0.5 * qtheta_im12 * qtheta_im12 / (htheta_im12 * htheta_im12);
                    A.coeffRef(pq, ph_e) += -0.5 * theta * qtheta_ip12 * qtheta_ip12 / (htheta_ip12 * htheta_ip12);
                    //
                    A.coeffRef(pq, ph_w + 1) += -theta * qtheta_im12 / htheta_im12;
                    A.coeffRef(pq, ph + 1) += -theta * qtheta_im12 / htheta_im12 + theta * qtheta_ip12 / htheta_ip12;
                    A.coeffRef(pq, ph_e + 1) += +theta * qtheta_ip12 / htheta_ip12;
                    //
                    double rhs_convection = -(qtheta_ip12 * qtheta_ip12 / htheta_ip12 - qtheta_im12 * qtheta_im12 / htheta_im12);
                    rhs[pq] += rhs_convection;
                }
                if (momentum_bed_shear_stress)
                {
                    double cf_im1 = cf[i - 1];
                    double cf_i = cf[i];
                    double cf_ip1 = cf[i + 1];
                    double htheta_im14 = scv(htheta_i, htheta_im1);
                    double qtheta_im14 = scv(qtheta_i, qtheta_im1);
                    double abs_qtheta_im14 = scv(Fabs(qtheta_i, eps_fabs), Fabs(qtheta_im1, eps_fabs));
                    double htheta_ip14 = scv(htheta_i, htheta_ip1);
                    double qtheta_ip14 = scv(qtheta_i, qtheta_ip1);
                    double abs_qtheta_ip14 = scv(Fabs(qtheta_i, eps_fabs), Fabs(qtheta_ip1, eps_fabs));
                    double cf_im14 = scv(cf_i, cf_im1);
                    double cf_ip14 = scv(cf_i, cf_ip1);
                    //
                    // bed_shear_stress
                    //
                    A.coeffRef(pq, ph_w) += -0.5 * 0.25 * theta * cf_im14 * abs_qtheta_im14 * 2. * qtheta_im14 / (htheta_im14 * htheta_im14 * htheta_im14);
                    A.coeffRef(pq, ph) += -0.5 * 0.25 * theta * cf_im14 * abs_qtheta_im14 * 2. * qtheta_im14 / (htheta_im14 * htheta_im14 * htheta_im14)
                        - 0.5 * 0.25 * theta * cf_ip14 * abs_qtheta_ip14 * 2. * qtheta_ip14 / (htheta_ip14 * htheta_ip14 * htheta_ip14);
                    A.coeffRef(pq, ph_e) += -0.5 * 0.25 * theta * cf_ip14 * abs_qtheta_ip14 * 2. * qtheta_ip14 / (htheta_ip14 * htheta_ip14 * htheta_ip14);
                    //
                    double J1_im14 = cf_im14 * qtheta_im14 / (htheta_im14 * htheta_im14) * dxinv * (Fabs(qtheta_i, eps_fabs) - Fabs(qtheta_im1, eps_fabs));
                    double J1_ip14 = cf_ip14 * qtheta_ip14 / (htheta_ip14 * htheta_ip14) * dxinv * (Fabs(qtheta_ip1, eps_fabs) - Fabs(qtheta_i, eps_fabs));
                    double J2_im14 = cf_im14 * abs_qtheta_im14 / (htheta_im14 * htheta_im14);
                    double J2_ip14 = cf_ip14 * abs_qtheta_ip14 / (htheta_ip14 * htheta_ip14);
                    A.coeffRef(pq, ph_w + 1) += dx * 0.5 * 0.25 * (theta * J1_im14 + theta * J2_im14);
                    A.coeffRef(pq, ph + 1) += dx * 0.5 * (0.75 * (theta * J1_im14 + theta * J2_im14) + 0.75 * (theta * J1_ip14 + theta * J2_ip14));
                    A.coeffRef(pq, ph_e + 1) += dx * 0.5 * 0.25 * (theta * J1_ip14 + theta * J2_ip14);
                    //
                    double rhs_bed_stress = -(
                        0.5 * dx * cf_im14 * abs_qtheta_im14 * qtheta_im14 / (htheta_im14 * htheta_im14) +
                        0.5 * dx * cf_ip14 * abs_qtheta_ip14 * qtheta_ip14 / (htheta_ip14 * htheta_ip14)
                        );
                    rhs[pq] += rhs_bed_stress;
                }
                if (momentum_viscosity)
                {
                    //
                    // viscosity
                    // 
                    double visc_im12 = -0.5 * (visc[i - 1] + visc[i]);
                    double visc_ip12 = -0.5 * (visc[i] + visc[i + 1]);

                    double A_im12 = visc_im12 * theta;
                    double A_ip12 = visc_ip12 * theta;
                    double B_im12 = -visc_im12 * theta / (htheta_im12) * dxinv * (htheta_i - htheta_im1);
                    double B_ip12 = -visc_ip12 * theta / (htheta_ip12) * dxinv * (htheta_ip1 - htheta_i);
                    double C_im12 = visc_im12 * theta * qtheta_im12 / (htheta_im12 * htheta_im12) * dxinv * (htheta_i - htheta_im1);
                    double C_ip12 = visc_ip12 * theta * qtheta_ip12 / (htheta_ip12 * htheta_ip12) * dxinv * (htheta_ip1 - htheta_i);
                    double D_im12 = -visc_im12 * theta * qtheta_im12 / htheta_im12;
                    double D_ip12 = -visc_ip12 * theta * qtheta_ip12 / htheta_ip12;
                    // 
                    // nu theta d(dzeta)/dx
                    A.coeffRef(pq, ph_w + 1) +=  dxinv * A_im12;
                    A.coeffRef(pq, ph + 1)   += -dxinv * A_im12;
                    A.coeffRef(pq, ph + 1)   += -dxinv * A_ip12;
                    A.coeffRef(pq, ph_e + 1) +=  dxinv * A_ip12;
                    //
                    // - nu theta q/h dq
                    A.coeffRef(pq, ph_w + 1) += -0.5 * B_im12;
                    A.coeffRef(pq, ph + 1)   += -0.5 * B_im12;
                    A.coeffRef(pq, ph + 1)   +=  0.5 * B_ip12;
                    A.coeffRef(pq, ph_e + 1) +=  0.5 * B_ip12;
                    //
                    // nu theta q/h^2 d(h)/dx dh
                    A.coeffRef(pq, ph_w) += -0.5 * C_im12;
                    A.coeffRef(pq, ph)   += -0.5 * C_im12;
                    A.coeffRef(pq, ph)   +=  0.5 * C_ip12;
                    A.coeffRef(pq, ph_e) +=  0.5 * C_ip12;
                    //
                    // - nu theta q/h d(dh)/dx
                    A.coeffRef(pq, ph_w) +=  dxinv * D_im12;
                    A.coeffRef(pq, ph)   += -dxinv * D_im12;
                    A.coeffRef(pq, ph)   += -dxinv * D_ip12;
                    A.coeffRef(pq, ph_e) +=  dxinv * D_ip12;
                    //
                    rhs_viscosity[pq] = -(
                           visc_ip12 * dxinv * (qtheta_ip1 - qtheta_i) - visc_ip12 * qtheta_ip12 / htheta_ip12 * dxinv * (htheta_ip1 - htheta_i)
                        - (visc_im12 * dxinv * (qtheta_i - qtheta_im1) - visc_im12 * qtheta_im12 / htheta_im12 * dxinv * (htheta_i - htheta_im1)
                          )
                        );
                    rhs[pq] += rhs_viscosity[pq];
                }
            }
            {
                //
                // wwest boundary
                //
                int i = 0;
                int ph = 2 * p_index(i, 0, nx);  // continuity equation
                int ph_e = 2 * p_index(i + 1, 0, nx);  // continuity equation (east)
                int ph_ee = 2 * p_index(i + 2, 0, nx);  // continuity equation (east-east)
                int pq = ph + 1;  // u-momentum equation 

                w_ess[0] = w_nat[0];
                w_ess[1] = w_nat[1];
                w_ess[2] = w_nat[2];

                hn_i = hn[i];       // = h^{n}_{i}
                hn_ip1 = hn[i + 1];       // = h^{n}_{i+1}
                hn_ip2 = hn[i + 2];       // = h^{n}_{i+2}
                hp_i = hp[i];       // = h^{n+1,p}_{i}
                hp_ip1 = hp[i + 1];       // = h^{n+1,p}_{i+1}
                hp_ip2 = hp[i + 2];       // = h^{n+1,p}_{i+2}
                htheta_i = theta * hp_i + (1.0 - theta) * hn_i;
                htheta_ip1 = theta * hp_ip1 + (1.0 - theta) * hn_ip1;
                htheta_ip2 = theta * hp_ip2 + (1.0 - theta) * hn_ip2;

                qn_i = qn[i];       // = q^{n}_{i}
                qn_ip1 = qn[i + 1];       // = q^{n}_{i+1}
                qn_ip2 = qn[i + 2];       // = q^{n}_{i+2}
                qp_i = qp[i];       // = q^{n+1,p}_{i}
                qp_ip1 = qp[i + 1];       // = q^{n+1,p}_{i+1}
                qp_ip2 = qp[i + 2];       // = q^{n+1,p}_{i+1}
                qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                qtheta_ip1 = theta * qp_ip1 + (1.0 - theta) * qn_ip1;
                qtheta_ip2 = theta * qp_ip2 + (1.0 - theta) * qn_ip2;
                //
                // Essential boundary condition
                //
                hn_ip12 = w_ess[0] * hn_i + w_ess[1] * hn_ip1 + w_ess[2] * hn_ip2;
                hp_ip12 = w_ess[0] * hp_i + w_ess[1] * hp_ip1 + w_ess[2] * hp_ip2;
                htheta_ip12 = w_ess[0] * htheta_i + w_ess[1] * htheta_ip1 + w_ess[2] * htheta_ip2;
                qn_ip12 = w_ess[0] * qn_i + w_ess[1] * qn_ip1 + w_ess[2] * qn_ip2;
                qp_ip12 = w_ess[0] * qp_i + w_ess[1] * qp_ip1 + w_ess[2] * qp_ip2;
                qtheta_ip12 = w_ess[0] * qtheta_i + w_ess[1] * qtheta_ip1 + w_ess[2] * qtheta_ip2;
                double zb_ip12 = w_ess[0] * zb[i] + w_ess[1] * zb[i + 1] + w_ess[2] * zb[i + 2];

                A.coeffRef(ph, ph) = 0.0;
                A.coeffRef(ph, ph_e) = 0.0;
                A.coeffRef(ph, ph_ee) = 0.0;
                //
                A.coeffRef(ph, ph + 1) = 0.0;
                A.coeffRef(ph, ph_e + 1) = 0.0;
                A.coeffRef(ph, ph_ee + 1) = 0.0;
                //
                rhs[ph] = 0.0;
                //
                double h_given = wz_bnd - zb_ip12;
                double h_infty = h_given;  // s_offset - zb_ip12;
                double c_wave = std::sqrt(g * h_infty);
                if (bc_absorbing[BC_WEST])
                {
                    if (bc_type[BC_WEST] == "mooiman")
                    {
                        // q
                        A.coeffRef(ph, ph + 1) = w_ess[0];
                        A.coeffRef(ph, ph_e + 1) = w_ess[1];
                        A.coeffRef(ph, ph_ee + 1) = w_ess[2];
                        //  c_wave * h
                        A.coeffRef(ph, ph) = w_ess[0] * c_wave;
                        A.coeffRef(ph, ph_e) = w_ess[1] * c_wave;
                        A.coeffRef(ph, ph_ee) = w_ess[2] * c_wave;
                        //
                        rhs[ph] = -(qp_ip12  + c_wave * (hp_ip12 - h_infty))
                            + 2. * c_wave * wz_bnd
                            + 2. * h_infty * wu_bnd
                            + 2. * wq_bnd
                            + 2. * wh_bnd
                            ;
                    }
                    else
                    {
                        //
                        // momentum + c_wave * continuity (ingoing signal)
                        // 
                        double con_fac = c_wave;
                        if (momentum_convection) { con_fac = c_wave + qp_ip12 / hp_ip12; }
                        // 
                        // momentum part
                        // 
                        A.coeffRef(ph, ph + 1)    = dtinv * w_ess[0];
                        A.coeffRef(ph, ph_e + 1)  = dtinv * w_ess[1];
                        A.coeffRef(ph, ph_ee + 1) = dtinv * w_ess[2];
                        //
                        // continuity part (added and multiplied by c_wave)
                        //
                        A.coeffRef(ph, ph)    = dtinv * con_fac * w_ess[0];
                        A.coeffRef(ph, ph_e)  = dtinv * con_fac * w_ess[1];
                        A.coeffRef(ph, ph_ee) = dtinv * con_fac * w_ess[2];
                        //
                        double dhdt = dtinv * (hp_ip12 - hn_ip12);
                        double dqdt = dtinv * (qp_ip12 - qn_ip12);
                        rhs[ph] = - (dqdt + con_fac * dhdt);

                        double corr_term = 0.0;
                        if (bc_vars[BC_WEST] == "zeta")
                        {
                            A.coeffRef(ph, ph)    += dtinv * w_ess[0] + eps_bc_corr * w_ess[0];
                            A.coeffRef(ph, ph_e)  += dtinv * w_ess[1] + eps_bc_corr * w_ess[1];
                            A.coeffRef(ph, ph_ee) += dtinv * w_ess[2] + eps_bc_corr * w_ess[2];
                            corr_term = -dhdt - eps_bc_corr * (hp_ip12 - (wz_bnd - zb_ip12));
                            rhs[ph] += corr_term;
                        }
                        if (bc_vars[BC_WEST] == "q")
                        {
                            A.coeffRef(ph, ph + 1)    += dtinv * w_ess[0] + eps_bc_corr * w_ess[0];
                            A.coeffRef(ph, ph_e + 1)  += dtinv * w_ess[1] + eps_bc_corr * w_ess[1];
                            A.coeffRef(ph, ph_ee + 1) += dtinv * w_ess[2] + eps_bc_corr * w_ess[2];
                            corr_term = -dqdt - eps_bc_corr * (qp_ip12 - wq_bnd);
                            rhs[ph] += corr_term;
                        }
                    }
                }
                else
                {
                    // Not absorbing
                    double corr_term = 0.0;
                    w_ess[0] = 1. / 12.;
                    w_ess[1] = 10. / 12.;
                    w_ess[2] = 1. / 12.;

                    if (bc_vars[BC_WEST] == "zeta")
                    {
                        A.coeffRef(ph, ph) = eps_bc_corr * w_ess[0];
                        A.coeffRef(ph, ph_e) = eps_bc_corr * w_ess[1];
                        A.coeffRef(ph, ph_ee) = eps_bc_corr * w_ess[2];
                        corr_term = -eps_bc_corr * (hp_ip12 - (wz_bnd - zb_ip12));
                        rhs[ph] += corr_term;
                    }
                    if (bc_vars[BC_WEST] == "q")
                    {
                        A.coeffRef(ph, ph + 1) = eps_bc_corr * w_ess[0];
                        A.coeffRef(ph, ph_e + 1) = eps_bc_corr * w_ess[1];
                        A.coeffRef(ph, ph_ee + 1) = eps_bc_corr * w_ess[2];
                        corr_term = -eps_bc_corr * (qp_ip12 - wq_bnd);
                        rhs[ph] += corr_term;
                    }
                }
                // 
                // Natural boundary condition
                // 
                hn_ip12 = w_nat[0] * hn_i + w_nat[1] * hn_ip1 + w_nat[2] * hn_ip2;
                hp_ip12 = w_nat[0] * hp_i + w_nat[1] * hp_ip1 + w_nat[2] * hp_ip2;
                htheta_ip12 = w_nat[0] * htheta_i + w_nat[1] * htheta_ip1 + w_nat[2] * htheta_ip2;

                qn_ip12 = w_nat[0] * qn_i + w_nat[1] * qn_ip1 + w_nat[2] * qn_ip2;
                qp_ip12 = w_nat[0] * qp_i + w_nat[1] * qp_ip1 + w_nat[2] * qp_ip2;
                qtheta_ip12 = w_nat[0] * qtheta_i + w_nat[1] * qtheta_ip1 + w_nat[2] * qtheta_ip2;

                A.coeffRef(pq, ph) = 0.0;
                A.coeffRef(pq, ph_e) = 0.0;
                A.coeffRef(pq, ph_ee) = 0.0;
                //
                A.coeffRef(pq, ph + 1) = 0.0;
                A.coeffRef(pq, ph_e + 1) = 0.0;
                A.coeffRef(pq, ph_ee + 1) = 0.0;
                //
                rhs[pq] = 0.0;
                //
                // momentum - c_wave * continuity
                // 
                double dhdt = dtinv * (hp_i - hn_i) * w_nat[0]
                    + dtinv * (hp_ip1 - hn_ip1) * w_nat[1]
                    + dtinv * (hp_ip2 - hn_ip2) * w_nat[2];
                double dqdx = dxinv * (qtheta_ip1 - qtheta_i);
                double dqdt = dtinv * (qp_i - qn_i) * w_nat[0]
                    + dtinv * (qp_ip1 - qn_ip1) * w_nat[1]
                    + dtinv * (qp_ip2 - qn_ip2) * w_nat[2];
                double dzetadx = dxinv * g * htheta_ip12 * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]);
                // 
                // momentum part dq/dt + gh d(zeta)/dt
                // 
                A.coeffRef(pq, ph) += w_nat[0] * dxinv * theta * g * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]) - dxinv * theta * g * htheta_ip12;
                A.coeffRef(pq, ph_e) += w_nat[1] * dxinv * theta * g * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]) + dxinv * theta * g * htheta_ip12;
                A.coeffRef(pq, ph_ee) += w_nat[2] * dxinv * theta * g * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]);
                A.coeffRef(pq, ph + 1) += dtinv * w_nat[0];
                A.coeffRef(pq, ph_e + 1) += dtinv * w_nat[1];
                A.coeffRef(pq, ph_ee + 1) += dtinv * w_nat[2];
                rhs[pq] += -(dqdt + dzetadx);
                if (momentum_convection) // 
                {
                    double aa = - dxinv * 2. * qtheta_ip12 / (htheta_ip12 * htheta_ip12) * (qtheta_ip1 - qtheta_i)
                        + dxinv * 2. * (qtheta_ip12 * qtheta_ip12) / (htheta_ip12 * htheta_ip12 * htheta_ip12) * (htheta_ip1 - htheta_i);
                    double bb = dxinv * 2. / htheta_ip12 * (qtheta_ip1 - qtheta_i) - dxinv * 2. * qtheta_ip12 / (htheta_ip12 * htheta_ip12) * (htheta_ip1 - htheta_i);
                    double cc = -(qtheta_ip12 * qtheta_ip12) / (htheta_ip12 * htheta_ip12);
                    double dd = 2. * qtheta_ip12 / htheta_ip12;

                    A.coeffRef(pq, ph) += theta * aa * w_nat[0] + dxinv * theta * cc;
                    A.coeffRef(pq, ph_e) += theta * aa * w_nat[1] - dxinv * theta * cc;
                    A.coeffRef(pq, ph_ee) += theta * aa * w_nat[2];
                    A.coeffRef(pq, ph + 1) += theta * bb * w_nat[0] - dxinv * theta * dd;
                    A.coeffRef(pq, ph_e + 1) += theta * bb * w_nat[1] + dxinv * theta * dd;
                    A.coeffRef(pq, ph_ee + 1) += theta * bb * w_nat[2];
                    rhs[pq] += -(
                          dxinv * dd * (qtheta_ip1 - qtheta_i) +
                        + dxinv * cc * (htheta_ip1 - htheta_i)
                        );
                }
                if (momentum_bed_shear_stress) // 
                {
                    double cf_i = cf[i];
                    double cf_ip1 = cf[i + 1];

                    double cf_ip12 = scv(cf_i, cf_ip1);
                    double htheta_ip12 = scv(htheta_i, htheta_ip1);
                    double qtheta_ip12 = scv(qtheta_i, qtheta_ip1);
                    double abs_qtheta_ip12 = scv(Fabs(qtheta_i, eps_fabs), Fabs(qtheta_ip1, eps_fabs));

                    A.coeffRef(pq, ph)   += -0.5 * theta * cf_ip12 * 2. * qtheta_ip12 * abs_qtheta_ip12 / (htheta_ip12 * htheta_ip12 * htheta_ip12);
                    A.coeffRef(pq, ph_e) += -0.5 * theta * cf_ip12 * 2. * qtheta_ip12 * abs_qtheta_ip12 / (htheta_ip12 * htheta_ip12 * htheta_ip12);

                    double J1_ip12 = cf_ip12 * qtheta_ip12 / (htheta_ip12 * htheta_ip12) * dxinv * (Fabs(qtheta_ip1, eps_fabs) - Fabs(qtheta_i, eps_fabs));
                    double J2_ip12 = cf_ip12 * abs_qtheta_ip12 / (htheta_ip12 * htheta_ip12);
                    A.coeffRef(pq, ph + 1) += 0.5 * theta * (J1_ip12 + J2_ip12);
                    A.coeffRef(pq, ph_e + 1) += 0.5 * theta * (J1_ip12 + J2_ip12);
                    rhs[pq] += -(
                        cf_ip12 * qtheta_ip12 * abs_qtheta_ip12 / (htheta_ip12 * htheta_ip12)
                        );
                }
                if (momentum_viscosity) // 
                {
                }
                //
                // continuity part (added and multiplied by -c_wave)
                //
                double con_fac = c_wave;
                if (momentum_convection) { con_fac = c_wave - qp_ip12 / hp_ip12; }
                A.coeffRef(pq, ph) += -con_fac * dtinv * w_nat[0];
                A.coeffRef(pq, ph_e) += -con_fac * dtinv * w_nat[1];
                A.coeffRef(pq, ph_ee) += -con_fac * dtinv * w_nat[2];
                A.coeffRef(pq, ph + 1) += -con_fac * -dxinv * theta;
                A.coeffRef(pq, ph_e + 1) += -con_fac * dxinv * theta;
                A.coeffRef(pq, ph_ee + 1) += 0.0;
                rhs[pq] += con_fac * (dhdt + dqdx);
            }
            {
                //
                // eeast boundary
                //
                // Right boundary (rhs only) for both continuity and momentum
                //
                int i = nx - 1;
                int ph = 2 * p_index(i, 0, nx);  // continuity equation
                int ph_w = 2 * p_index(i - 1, 0, nx);  // continuity equation
                int ph_ww = 2 * p_index(i - 2, 0, nx);  // continuity equation
                int pq = ph + 1;  // u-momentum equation

                w_ess[0] = w_nat[0];
                w_ess[1] = w_nat[1];
                w_ess[2] = w_nat[2];

                hn_i = hn[i];       // = h^{n}_{i}
                hn_im1 = hn[i - 1];       // = h^{n}_{i-1}
                hn_im2 = hn[i - 2];       // = h^{n}_{i-2}
                hp_i = hp[i];       // = h^{n+1,p}_{i}
                hp_im1 = hp[i - 1];       // = h^{n+1,p}_{i-1}
                hp_im2 = hp[i - 2];       // = h^{n+1,p}_{i-2}
                htheta_i = theta * hp_i + (1.0 - theta) * hn_i;
                htheta_im1 = theta * hp_im1 + (1.0 - theta) * hn_im1;
                htheta_im2 = theta * hp_im2 + (1.0 - theta) * hn_im2;

                qn_i = qn[i];       // = q^{n}_{i}
                qn_im1 = qn[i - 1];       // = q^{n}_{i-1}
                qn_im2 = qn[i - 2];       // = q^{n}_{i-1}
                qp_i = qp[i];       // = q^{n+1,p}_{i}
                qp_im1 = qp[i - 1];       // = q^{n+1,p}_{i-1}
                qp_im2 = qp[i - 2];       // = q^{n+1,p}_{i-2}
                qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                qtheta_im1 = theta * qp_im1 + (1.0 - theta) * qn_im1;
                qtheta_im2 = theta * qp_im2 + (1.0 - theta) * qn_im2;
                //
                // Essential boundary condition
                //
                hn_im12 = w_ess[0] * hn_i + w_ess[1] * hn_im1 + w_ess[2] * hn_im2;
                hp_im12 = w_ess[0] * hp_i + w_ess[1] * hp_im1 + w_ess[2] * hp_im2;
                htheta_im12 = w_ess[0] * htheta_i + w_ess[1] * htheta_im1 + w_ess[2] * htheta_im2;
                qn_im12 = w_ess[0] * qn_i + w_ess[1] * qn_im1 + w_ess[2] * qn_im2;
                qp_im12 = w_ess[0] * qp_i + w_ess[1] * qp_im1 + w_ess[2] * qp_im2;
                qtheta_im12 = w_ess[0] * qtheta_i + w_ess[1] * qtheta_im1 + w_ess[2] * qtheta_im2;

                double zb_im12 = w_ess[0] * zb[i] + w_ess[1] * zb[i - 1] + w_ess[2] * zb[i - 2];

                A.coeffRef(ph, ph) = 0.0;
                A.coeffRef(ph, ph_w) = 0.0;
                A.coeffRef(ph, ph_ww) = 0.0;
                //
                A.coeffRef(ph, ph + 1) = 0.0;
                A.coeffRef(ph, ph_w + 1) = 0.0;
                A.coeffRef(ph, ph_ww + 1) = 0.0;
                //
                rhs[ph] = 0.0;

                double h_given = ez_bnd - zb_im12;
                double h_infty = h_given;  // s_offset - zb_im12;
                double c_wave = std::sqrt(g * h_infty);

                if (bc_absorbing[BC_EAST])
                {
                    if (bc_type[BC_EAST] == "mooiman")
                    {
                        // q
                        A.coeffRef(ph, ph + 1) = w_ess[0];
                        A.coeffRef(ph, ph_w + 1) = w_ess[1];
                        A.coeffRef(ph, ph_ww + 1) = w_ess[2];
                        // c_wave * h
                        A.coeffRef(ph, ph) = -w_ess[0] * c_wave;
                        A.coeffRef(ph, ph_w) = -w_ess[1] * c_wave;
                        A.coeffRef(ph, ph_ww) = -w_ess[2] * c_wave;
                        //
                        rhs[ph] = -(qp_im12 - c_wave * (hp_im12 - h_infty))
                            - 2. * c_wave * ez_bnd
                            - 2. * h_infty * eu_bnd
                            - 2. * eq_bnd
                            - 2. * eh_bnd
                            ;
                    }
                    else
                    {
                        //
                        // momentum - c_wave * continuity (ingoing signal)
                        // 
                        double con_fac = c_wave;
                        if (momentum_convection) { con_fac = c_wave - qp_ip12 / hp_ip12; }
                        // 
                        // momentum part
                        // 
                        A.coeffRef(ph, ph + 1) = dtinv * w_ess[0];
                        A.coeffRef(ph, ph_w + 1) = dtinv * w_ess[1];
                        A.coeffRef(ph, ph_ww + 1) = dtinv * w_ess[2];
                        //
                        // continuity part (added and multiplied by -c_wave)
                        //
                        A.coeffRef(ph, ph) = dtinv * -con_fac * w_ess[0];
                        A.coeffRef(ph, ph_w) = dtinv * -con_fac * w_ess[1];
                        A.coeffRef(ph, ph_ww) = dtinv * -con_fac * w_ess[2];
                        //
                        double dhdt = dtinv * (hp_im12 - hn_im12);
                        double dqdt = dtinv * (qp_im12 - qn_im12);
                        rhs[ph] = -(dqdt - con_fac * dhdt);

                        double corr_term = 0.0;
                        if (bc_vars[BC_EAST] == "zeta")
                        {
                            if (stationary) { sign = -1.0; }
                            A.coeffRef(ph, ph)    += dtinv * w_ess[0] + eps_bc_corr * w_ess[0];
                            A.coeffRef(ph, ph_w)  += dtinv * w_ess[1] + eps_bc_corr * w_ess[1];
                            A.coeffRef(ph, ph_ww) += dtinv * w_ess[2] + eps_bc_corr * w_ess[2];
                            corr_term = -dhdt + sign * eps_bc_corr * (hp_im12 - (ez_bnd - zb_im12));
                            rhs[ph] += corr_term;
                            sign = 1.0;
                        }
                        if (bc_vars[BC_EAST] == "q")
                        {
                            if (stationary) { sign = -1.0; }
                            A.coeffRef(ph, ph + 1)    += dtinv * w_ess[0] - eps_bc_corr * w_ess[0];
                            A.coeffRef(ph, ph_w + 1)  += dtinv * w_ess[1] - eps_bc_corr * w_ess[1];
                            A.coeffRef(ph, ph_ww + 1) += dtinv * w_ess[2] - eps_bc_corr * w_ess[2];
                            corr_term = -dqdt - sign * eps_bc_corr * (qp_im12 - eq_bnd);
                            rhs[ph] += corr_term;
                            sign = 1.0;
                        }
                    }
                }
                else
                {
                    // Not absorbing
                    double corr_term = 0.0;
                    w_ess[0] = 1. / 12.;
                    w_ess[1] = 10. / 12.;
                    w_ess[2] = 1. / 12.;
                    if (bc_vars[BC_EAST] == "zeta")
                    {
                        A.coeffRef(ph, ph) = eps_bc_corr * w_ess[0];
                        A.coeffRef(ph, ph_w) = eps_bc_corr * w_ess[1];
                        A.coeffRef(ph, ph_ww) = eps_bc_corr * w_ess[2];
                        corr_term = -eps_bc_corr * (hp_im12 - (ez_bnd - zb_im12));
                        rhs[ph] = corr_term;
                    }
                    if (bc_vars[BC_EAST] == "q")
                    {
                        A.coeffRef(ph, ph + 1) = eps_bc_corr * w_ess[0];
                        A.coeffRef(ph, ph_w + 1) = eps_bc_corr * w_ess[1];
                        A.coeffRef(ph, ph_ww + 1) = eps_bc_corr * w_ess[2];
                        corr_term = -eps_bc_corr * (qp_im12 - eq_bnd);
                        rhs[ph] = corr_term;
                    }
                }
                //
                // Natural boundary condition
                //  
                hn_im12 = w_nat[0] * hn_i + w_nat[1] * hn_im1 + w_nat[2] * hn_im2;
                hp_im12 = w_nat[0] * hp_i + w_nat[1] * hp_im1 + w_nat[2] * hp_im2;
                htheta_im12 = w_nat[0] * htheta_i + w_nat[1] * htheta_im1 + w_nat[2] * htheta_im2;

                qn_im12 = w_nat[0] * qn_i + w_nat[1] * qn_im1 + w_nat[2] * qn_im2;
                qp_im12 = w_nat[0] * qp_i + w_nat[1] * qp_im1 + w_nat[2] * qp_im2;
                qtheta_im12 = w_nat[0] * qtheta_i + w_nat[1] * qtheta_im1 + w_nat[2] * qtheta_im2;
                //            
                A.coeffRef(pq, ph) = 0.0;
                A.coeffRef(pq, ph_w) = 0.0;
                A.coeffRef(pq, ph_ww) = 0.0;
                //
                A.coeffRef(pq, ph + 1) = 0.0;
                A.coeffRef(pq, ph_w + 1) = 0.0;
                A.coeffRef(pq, ph_ww + 1) = 0.0;
                //
                rhs[pq] = 0.0;
                //
                // momentum + c_wave * continuity
                // 
                double dhdt = dtinv * (hp_i - hn_i) * w_nat[0]
                    + dtinv * (hp_im1 - hn_im1) * w_nat[1]
                    + dtinv * (hp_im2 - hn_im2) * w_nat[2];
                double dqdx = dxinv * (qtheta_i - qtheta_im1);
                double dqdt = dtinv * (qp_i - qn_i) * w_nat[0]
                    + dtinv * (qp_im1 - qn_im1) * w_nat[1]
                    + dtinv * (qp_im2 - qn_im2) * w_nat[2];
                double dzetadx = dxinv * g * htheta_im12 * (htheta_i + zb[i] - htheta_im1 - zb[i - 1]);
                // 
                // momentum part dq/dt + gh d(zeta)/dx
                // 
                A.coeffRef(pq, ph) += w_nat[0] * dxinv * theta * g * (htheta_i + zb[i] - htheta_im1 - zb[i - 1]) + dxinv * theta * g * htheta_im12;
                A.coeffRef(pq, ph_w) += w_nat[1] * dxinv * theta * g * (htheta_i + zb[i] - htheta_im1 - zb[i - 1]) - dxinv * theta * g * htheta_im12;
                A.coeffRef(pq, ph_ww) += w_nat[2] * dxinv * theta * g * (htheta_i + zb[i] - htheta_im1 - zb[i - 1]);
                A.coeffRef(pq, ph + 1) += dtinv * w_nat[0];
                A.coeffRef(pq, ph_w + 1) += dtinv * w_nat[1];
                A.coeffRef(pq, ph_ww + 1) += dtinv * w_nat[2];
                rhs[pq] += -(dqdt + dzetadx);
                if (momentum_convection) // 
                {
                    double aa = - dxinv * 2. * qtheta_im12 / (htheta_im12 * htheta_im12) * (qtheta_i - qtheta_im1)
                        + dxinv * 2. * (qtheta_im12 * qtheta_im12) / (htheta_im12 * htheta_im12 * htheta_im12) * (htheta_i - htheta_im1);
                    double bb = dxinv * 2. / htheta_im12 * (qtheta_i - qtheta_im1) - dxinv * 2. * qtheta_im12 / (htheta_im12 * htheta_im12) * (htheta_i - htheta_im1);
                    double cc = -(qtheta_im12 * qtheta_im12) / (htheta_im12 * htheta_im12);
                    double dd = 2. * qtheta_im12 / htheta_im12;

                    A.coeffRef(pq, ph) += theta * aa * w_nat[0] - dxinv * theta * cc;
                    A.coeffRef(pq, ph_w) += theta * aa * w_nat[1] + dxinv * theta * cc;
                    A.coeffRef(pq, ph_ww) += theta * aa * w_nat[2];
                    A.coeffRef(pq, ph + 1) += theta * bb * w_nat[0] + dxinv * theta * dd;
                    A.coeffRef(pq, ph_w + 1) += theta * bb * w_nat[1] - dxinv * theta * dd;
                    A.coeffRef(pq, ph_ww + 1) += theta * bb * w_nat[2];
                    rhs[pq] += -(
                        dxinv * dd * (qtheta_i - qtheta_im1)
                        + dxinv * cc * (htheta_i - htheta_im1)
                        );
                }
                if (momentum_bed_shear_stress) // 
                {
                    double cf_im1 = cf[i - 1];
                    double cf_i = cf[i];
                    double htheta_im12 = scv(htheta_i, htheta_im1);
                    double qtheta_im12 = scv(qtheta_i, qtheta_im1);
                    double abs_qtheta_im12 = scv(Fabs(qtheta_i, eps_fabs), Fabs(qtheta_im1, eps_fabs));
                    double cf_im12 = scv(cf_i, cf_im1);
                    //
                    // bed_shear_stress
                    //
                    A.coeffRef(pq, ph  ) += -0.5 * theta * cf_im12 * 2. * qtheta_im12 * abs_qtheta_im12 / (htheta_im12 * htheta_im12 * htheta_im12);
                    A.coeffRef(pq, ph_w) += -0.5 * theta * cf_im12 * 2. * qtheta_im12 * abs_qtheta_im12 / (htheta_im12 * htheta_im12 * htheta_im12);
                    //
                    double J1_im12 = cf_im12 * qtheta_im12 / (htheta_im12 * htheta_im12) * dxinv * (Fabs(qtheta_i, eps_fabs) - Fabs(qtheta_im1, eps_fabs));
                    double J2_im12 = cf_im12 * abs_qtheta_im12 / (htheta_im12 * htheta_im12);
                    A.coeffRef(pq, ph + 1) += 0.5 * theta * (J1_im12 + J2_im12);
                    A.coeffRef(pq, ph_w + 1) += 0.5 * theta * (J1_im12 + J2_im12);
                    //
                    double rhs_bed_stress = -(
                        cf_im12 * qtheta_im12 * abs_qtheta_im12 / (htheta_im12 * htheta_im12)
                        );
                    rhs[pq] += rhs_bed_stress;
                }
                if (momentum_viscosity) // 
                {
                }
                //
                // continuity part (added and multiplied by +c_wave)
                //
                double con_fac = c_wave;
                if (momentum_convection) { con_fac = c_wave + qp_im12 / hp_im12; }
                A.coeffRef(pq, ph) += con_fac * dtinv * w_nat[0];
                A.coeffRef(pq, ph_w) += con_fac * dtinv * w_nat[1];
                A.coeffRef(pq, ph_ww) += con_fac * dtinv * w_nat[2];
                A.coeffRef(pq, ph + 1) += con_fac * dxinv * theta;
                A.coeffRef(pq, ph_w + 1) += con_fac * -dxinv * theta;
                A.coeffRef(pq, ph_w + 1) += 0.0;
                rhs[pq] += -con_fac * (dhdt + dqdx);
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
            Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
            solver.compute(A);
            //solution = solver.solve(rhs);
            solver.setTolerance(eps_bicgstab);
            solution = solver.solveWithGuess(rhs, solution);
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
                         << "    estimated error:" << solver.error() << std::endl;
            }

            // The new solution is the previous iterant plus the delta
            dh_max = 0.0;
            dh_maxi = 0;
            dq_max = 0.0;
            dq_maxi = 0;
            for (int i = 0; i < nx; ++i)
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

            if (regularization_iter)
            {
                START_TIMER(Regularization_iter_loop);
                if (momentum_viscosity)
                {
                    for (int i = 0; i < nx; ++i)
                    {
                        u[i] = qp[i] / hp[i];
                    }
                    //(void)regularization->first_derivative(psi, visc_reg, u, dx);
                    (void)regularization->given_function(tmp1, psi, eq8, u, dx, c_psi, use_eq8);
                    for (int i = 0; i < nx; ++i)
                    {
                        qp[i] = tmp1[i] * hp[i];
                        visc[i] = visc_reg[i] * std::abs(psi[i]);
                        pe[i] = qp[i] / hp[i] * dx / visc[i];
                    }
                }
                STOP_TIMER(Regularization_iter_loop);
            }
            if (logging == "matrix")
            {
                log_file << "=== Matrix ============================================" << std::endl;
                log_file << Eigen::MatrixXd(A) << std::endl;
                //log_file << "=== Eigen values ======================================" << std::endl;
                //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).eigenvalues() << std::endl;
                log_file << "=== RHS, solution  ====================================" << std::endl;
                for (int i = 0; i < 2 * nx; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << rhs[i] << "' " << solution[i] << std::endl;
                }
                log_file << "=== hp, qp ============================================" << std::endl;
                for (int i = 0; i < nx; ++i)
                {
                    log_file << std::setprecision(8) << std::scientific << hp[i] << ", " << qp[i] << std::endl;
                }
                log_file << "=======================================================" << std::endl;
            }
            used_lin_solv_iter = std::max(used_lin_solv_iter, (int) solver.iterations());
            if (dh_max < eps_newton && dq_max < eps_newton)
            {
                used_newton_iter = iter;
                break;
            }
        }
        STOP_TIMER(Newton iteration);
        if (stationary)
        {
            std::cout << "stationary solution " << std::endl;
            log_file << "stationary solution " << std::endl;
            log_file << std::setprecision(8) << std::scientific
                << "    Iter: " << used_newton_iter + 1
                << "    Delta h^{n + 1,p + 1}: " << dh_max << " at: " << dh_maxi
                << "    Delta q^{n + 1,p + 1}: " << dq_max << " at: " << dq_maxi
                << std::endl;
        }
        else
        {
            if (std::fmod(time, 900.) == 0)
            {
                std::cout << std::fixed << std::setprecision(2) << tstart + time << ";   " << tstart + tstop << std::endl;
            }
            if (logging == "iterations" || logging == "matrix")
            {
                log_file << "time [sec]: " << std::setprecision(2) << std::scientific << time
                    << std::setprecision(8) << std::scientific
                    << "    Newton iterations  : " << used_newton_iter + 1
                    << "    Delta h^{n + 1,p + 1}: " << dh_max << " at: " << dh_maxi
                    << "    Delta q^{n + 1,p + 1}: " << dq_max << " at: " << dq_maxi
                    << std::endl;
            }

        }
        if (used_newton_iter + 1 == iter_max)
        {
            if (dh_max > eps_newton || dq_max > eps_newton)
            {
                log_file << "    ----    maximum number of iterations reached, probably not converged" << std::endl;
            }
        }
        for (int i = 0; i < nx; ++i)
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
            if (momentum_viscosity)
            {
                for (int i = 0; i < nx; ++i)
                {
                    u[i] = qp[i] / hp[i];
                }
                //(void)regularization->first_derivative(psi, visc_reg, u, dx);
                (void)regularization->given_function(tmp1, psi, eq8, u, dx, c_psi, use_eq8);
                for (int i = 0; i < nx; ++i)
                {
                    qp[i] = tmp1[i] * hp[i];
                    visc[i] = visc_reg[i] * std::abs(psi[i]);
                    pe[i] = qp[i] / hp[i] * dx / visc[i];
                }
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
            map_file->put_time_variable(map_names[7], nst_map, visc_reg);
            map_file->put_time_variable(map_names[8], nst_map, visc);
            map_file->put_time_variable(map_names[9], nst_map, psi);
            map_file->put_time_variable(map_names[10], nst_map, eq8);
            map_file->put_time_variable(map_names[11], nst_map, pe);
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
            std::vector<double> his_values = { hn[i_left], hn[i_mid_left], hn[i_mid], hn[i_mid_right],  hn[i_right] };
            his_file->put_variable(his_h_name, nst_his, his_values);
            his_values = { qn[i_left], qn[i_mid_left], qn[i_mid], qn[i_mid_right],  qn[i_right] };
            his_file->put_variable(his_q_name, nst_his, his_values);
            his_values = { s[i_left], s[i_mid_left], s[i_mid], s[i_mid_right],  s[i_right] };
            his_file->put_variable(his_s_name, nst_his, his_values);
            his_values = { u[i_left], u[i_mid_left], u[i_mid], u[i_mid_right],  u[i_right]};
            his_file->put_variable(his_u_name, nst_his, his_values);
            his_values = { zb[i_left], zb[i_mid_left], zb[i_mid], zb[i_mid_right], zb[i_right] };
            his_file->put_variable(his_zb_name, nst_his, his_values);

            //his_values = { riemann_pos[i_left], riemann_pos[i_left], riemann_pos[i_mid], riemann_pos[i_right], riemann_pos[i_right] };
            //his_file->put_variable(his_riemann_pos_name, nst_his, his_values);
            //his_values = { riemann_neg[i_left], riemann_neg[i_left], riemann_neg[i_mid], riemann_neg[i_right], riemann_neg[i_right] };
            //his_file->put_variable(his_riemann_neg_name, nst_his, his_values);

            his_values = { double(used_newton_iter + 1) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            //his_values = { double(used_lin_solv_iter) / double(used_newton_iter + 1) };
            his_values = { double(used_lin_solv_iter)};
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

    return 0;
}
int p_index(int i, int j, int nx)
{
    return j * nx + i;
}

int initialize_scalar(double alpha, std::vector<double>& value_in, std::vector<double>& value_out)
{
    int nx = value_in.size();
    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector
    Eigen::VectorXd rhs(nx);                // RHS vector

    int i = 0;
    A.coeffRef(i, i) = -alpha;
    A.coeffRef(i, i + 1) = alpha;
    rhs[i] = value_in[i];
    for (int i = 1; i < nx-1; ++i)
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

    for (int i = 0; i < nx; ++i)
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

    int nx = value.size();
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
    if (nrow != nx)
    {
        std::cout << "Dimension of bed_filename does not match;  " << nrow << " != " << nx << std::endl;
        return 1;
    }
    for (int k = 0; k < nx; ++k)
    {
        input_file >> dummy >> value[k];
    }
    input_file.close();
    return 0;
}

int initialize_bed_level(BED_LEVEL bed_type, std::vector<double>& x, std::vector<double>& zb, std::string & model_title, double depth)
{
    int status = 0;
    model_title = "Not specified";
    if (bed_type == BED_LEVEL::NONE)
    {
        std::cout << "Bed level type BED_LEVEL::NONE not supported" << std::endl;
        status = 1;
    }
    else if (bed_type == BED_LEVEL::FLAT)
    {
        std::stringstream depth_strm;
        depth_strm << std::fixed << std::setprecision(2) << depth;
        model_title = "Flat bed level at " + depth_strm.str() + " [m]";
        for (int i = 0; i < x.size(); ++i)
        {
            zb[i] = -depth;
        }
    }
    else if (bed_type == BED_LEVEL::SHOAL)
    {
        model_title = "Shoaling bedlevel from -10 [m] to -2.5 [m]";
        double slope_begin = 0.0;
        double slope_end = 1000.;
        for (int i = 0; i < x.size(); ++i)
        {
            if (x[i] < slope_begin)
                zb[i] = -10.0;
            else if(x[i] < slope_end)
                zb[i] = -10. + (x[i] - slope_begin) / (slope_end - slope_begin) * 7.5;
            else
                zb[i] = -2.5;
        }
    }
    else if (bed_type == BED_LEVEL::SLOPED)
    {
        model_title = "Sloped bed level 10 [cm] per [km] (par 5.1, Platzek2019)";
        // From -3.9 [m] to -4.0 [m]
        double ib = -1e-04;  // 0.1 mm per meter
        double slope_begin = x[1];
        double slope_end = x[x.size() - 2];
        for (int i = 0; i < x.size(); ++i)
        {
            zb[i] = -3.9 + (x[i] - slope_begin) * ib;
        }
    }
    else if (bed_type == BED_LEVEL::WEIR)
    {
        // Borsboom_development1Derrorminmovingadaptgridmethod_AdaptMethodLinesCRC2001.pdf
        model_title = "Weir: from -12 [m] to -5 [m] and from -5 [m] to -10 [m] (Borsboom_presCASATUE2023)";
        double x0 = x[1];
        double zb_def = -12.0;

        double slope_up_begin = 200.;
        double slope_up_end = 250.;
        double slope_down_begin = 350.;
        double slope_down_end = 450.;
        double zb_begin = zb_def;
        double zb_weir = zb_def + 7.0;
        double zb_end = zb_def + 2.0;
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] - x0;
            if (x[i] < slope_up_begin)
                zb[i] = zb_begin;
            else if (x[i] < slope_up_end)
                zb[i] = zb_begin + (x[i] - slope_up_begin) / (slope_up_end - slope_up_begin) * (zb_weir - zb_begin);
            else if (x[i] < slope_down_begin)
                zb[i] = zb_weir;
            else if (x[i] < slope_down_end)
                zb[i] = zb_weir + (x[i] - slope_down_begin) / (slope_down_end - slope_down_begin) * (zb_end - zb_weir);
            else if (x[i] >= slope_down_end)
                zb[i] = zb_end;
        }
    }
    else if (bed_type == BED_LEVEL::WAVY)
    {
        model_title = "Wavy bed level (par 5.2, Platzek2019)";
        int nx = x.size();
        double Lx = x[nx - 1] - x[0];
        double dx = Lx / (nx - 1);
        double zb_ref = -4.0;  // Mean depth for the wave bed level
        double Ab = 0.3;    // Amplitude of the bed forms
        double Nb = 25.0;  // Number of waves in bed level in the domain
        double fb = 2.0 * M_PI * Nb / ((nx-3) * dx);
        for (int i = 0; i < x.size(); ++i)
        {
            zb[i] = zb_ref - Ab * cos(fb * (x[i]- x[1]));
        }
        zb[0] = zb[1];
        zb[nx - 1] = zb[nx - 2];
        int a = 1;
    }
    else if (bed_type == BED_LEVEL::WAVY_SLOPED)
    {
        model_title = "Sloped wavy bed level (par 5.3, Platzek2019)";
        int nx = x.size();
        double Lx = x[nx - 1] - x[0];
        double dx = Lx / (nx - 1);
        double zb_ref = -3.9;  // Mean depth for the wave bed level
        double Ab = 0.3;    // Amplitude of the bed forms
        double Nb = 25.0;  // Number of waves in bed level in the domain
        double fb = 2.0 * M_PI * Nb / ((nx - 3) * dx);
        double ib = -1e-04;  // 0.1 mm per meter
        for (int i = 0; i < x.size(); ++i)
        {
            zb[i] = zb_ref - Ab * cos(fb * (x[i] - x[1])) + ib * (x[i] - x[1]);
        }
        zb[0] = zb[1];
        zb[nx - 1] = zb[nx - 2];

    }
    else
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Bed level type not supported" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
        status = 1; 
    }
    return status;
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
    std::string * file_name)   
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
int get_toml_array(toml::table tbl, std::string keyw, std::vector<std::string> & values)
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
int get_toml_array(toml::table tbl, std::string keyw, std::vector<double> & values)
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

