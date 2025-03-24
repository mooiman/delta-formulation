#define _USE_MATH_DEFINES
#include <ios>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <filesystem>

//
// Air pollution:
//     Hundsdorfer and verwer 2003, pg 7
//
//  du_1/dt = k1 u_3 - k2_u1
//  du_2/dt = k1 u_3 - k3 u_2 u_4 + sigma_2
//  du_3/dt = k3 u_2 u_4 - k1 u_3
//  du_4/dt = k2 u_1 - k3 u_2 u_4
//

#include "Eigen/IterativeLinearSolvers"
#include <toml.h>

#include "cfts.h"
#include "perf_timer.h"

double rhsfev(int ispecie, const std::vector<double>& u, const std::vector<double>& k);
void runge_kutta_4(double t0, double& t1, const std::vector<double>& u0, std::vector<double>& u1, const std::vector<double>& k, double dt);
void GetArguments(long argc, char** argv, std::string* file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

double k0(double t0);
void explicit_euler(double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k, 
    double dt, double sigma2);
void runge_kutta_4 (double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k,
    double dt, double sigma2);
double rhsfev(int, std::vector<double>&, std::vector<double>&, double);

enum class TIME_INTEGRATOR
{
    EULER_EXPLICIT,
    RUNGE_KUTTA_4,
    FULLY_IMPLICIT,
};

int main(int argc, char* argv[]) {
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
        std::cout << tbl << "\n";
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

    std::cout << "----------------------------" << std::endl;
    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Current directory   : " << std::filesystem::absolute(current_dir) << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;
    std::cout << "----------------------------" << std::endl;

    std::string out_file;
    std::stringstream ss;
    ss << "air_pollution";
    out_file = output_dir.string() + ss.str();

    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Two way chemical reaction, <no time integrator specified>");

    Eigen::SparseMatrix<double> A(4, 4);
    Eigen::VectorXd rhs(4);                // RHS vector [u1, u2, u3, u4]^n
    Eigen::VectorXd solution(4);           // solution vector [u1, u2, u3, u4]^{n}

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

    // Time
    log_file << "[Time]" << std::endl;
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    log_file << "tstart = " << tstart << std::endl;
    double tstop = tbl_chp["tstop"].value_or(double(300.));  // default 5 minutes
    log_file << "tstop = " << tstop << std::endl;

    // Numerics
    enum TIME_INTEGRATOR time_integration_type;
    log_file << std::endl << "[Numerics]" << std::endl;
    tbl_chp = *tbl["Numerics"].as_table();
    double dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary
    log_file << "dt = " << dt << std::endl;
    std::string integration = tbl_chp["integration"].value_or("None");
    if (integration == "delta_formulation") { time_integration_type = TIME_INTEGRATOR::FULLY_IMPLICIT; }
    else if (integration == "runge_kutta_4") { time_integration_type = TIME_INTEGRATOR::RUNGE_KUTTA_4; }
    log_file << "integration = " << integration << std::endl;
    int iter_max = tbl_chp["iter_max"].value_or(int(50));
    log_file << "iter_max  = " << iter_max << std::endl;
    double theta = tbl_chp["theta"].value_or(double(0.501));  // Implicitness factor (0.5 <= theta <= 1.0)
    log_file << "theta  = " << theta << std::endl;
    double eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));  // stop criterium for Newton iteration
    log_file << "eps_newton  = " << eps_newton << std::endl;
    double eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));  // stop criterium for BiCGStab iteration
    log_file << "eps_bicgstab  = " << eps_bicgstab << std::endl;

    // Physics
    log_file << std::endl << "[Physics]" << std::endl;
    tbl_chp = *tbl["Physics"].as_table();
    double k1 = tbl_chp["k1"].value_or(double(1.0));  // as in Hundsdorfer and Verwer
    log_file << "k1 = " << k1 << std::endl;
    double k2 = tbl_chp["k2"].value_or(double(10.));  // as in Hundsdorfer and Verwer 
    log_file << "k2 = " << k2 << std::endl;
    double k3 = tbl_chp["k3"].value_or(double(1.0));  // as in Hundsdorfer and Verwer
    log_file << "k3 = " << k3 << std::endl;
    double sigma2 = tbl_chp["sigma"].value_or(double(1.0e-07));  // as in Hundsdorfer and Verwer
    log_file << "sigma2 = " << sigma2 << std::endl;

    // Initial
    log_file << std::endl << "[Initial]" << std::endl;
    tbl_chp = *tbl["Initial"].as_table();
    double u0_ini = tbl_chp["O"].value_or(double(1.0));  // 0.d0 ! O Zuurstof-molekuul (Oxygen-molecule)
    log_file << "O = " << u0_ini << std::endl;
    double u1_ini = tbl_chp["NO"].value_or(double(10.));  // 1.3d08 ! NO Stikstof-oxide
    log_file << "NO = " << u1_ini << std::endl;
    double u2_ini = tbl_chp["NO2"].value_or(double(1.0));  // 5.0d11 ! NO2 Stikstof-di-oxide
    log_file << "NO2 = " << u2_ini << std::endl;
    double u3_ini = tbl_chp["O3"].value_or(double(1.0));  // 8.0d11 ! O3 = Ozon
    log_file << "O3 = " << u3_ini << std::endl;

    log_file << "=======================================================" << std::endl;

    std::vector<double> un(4, 0.);
    std::vector<double> up(4, 0.);
    std::vector<double> k(3, 0.);

    if (time_integration_type == TIME_INTEGRATOR::EULER_EXPLICIT)
    {
        model_title = "Air pollution, Euler explicit";
    }
    if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
    {
        model_title = "Air pollution, Runge-Kutta 4 (RK4)";
    }
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        model_title = "Air pollution, Delta formulation";
    }

    START_TIMERN(Main);

    bool stationary = false;
    if (dt == 0.0)
    {
        stationary = true;
    }

    double dtinv = 0.0;
    int total_time_steps = int((tstop - tstart) / dt) + 1;
    if (stationary) 
    { 
        dt = 0.0;                                         // Time step size [s]
        dtinv = 0.0;                                      // stationary solution
        theta = 1.0;                                      // Stationary solution
        tstop = 1.;
        total_time_steps = 2;                             // initial step (step 1), stationary result (step 2)
    }
    else
    {
        dtinv = 1. / dt;
    }
    int wrihis = 1; // int(std::max(dt, 1.0 / dt));      // write interval to his-file (every delta t)

    un[0] = u0_ini;     
    un[1] = u1_ini; 
    un[2] = u2_ini;
    un[3] = u3_ini; 
    k[0] = 0.0;     // Time dependent
    k[1] = k2;
    k[2] = k3; 

    double t0 = tstart * dt;
    double t1 =  t0;

    ////////////////////////////////////////////////////////////////////////////
    // Define time history file
    std::cout << "    Define his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations (corrected for virtual points)
    double x = 0.0;
    double y = 0.0;
    std::vector<double> x_obs = { x };
    std::vector<double> y_obs = { y };

    std::vector<std::string> obs_stations;
    obs_stations.push_back("Measurement station");
    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_O_name("O");
    std::string his_NO_name("NO");
    std::string his_NO2_name("NO2");
    std::string his_O3_name("O3");
    std::string his_newton_name("newton");
    std::string his_bicgstab_name("bicgstab");
    his_file->add_variable(his_O_name, "", "atomic oxygen", "-");
    his_file->add_variable(his_NO_name, "", "nitrogen oxide", "-");
    his_file->add_variable(his_NO2_name, "", "nitrogen dioxide", "-");
    his_file->add_variable(his_O3_name, "", "ozon", "-");
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        his_file->add_variable_without_location(his_newton_name, "", "Newton iteration", "-");
        his_file->add_variable_without_location(his_bicgstab_name, "", "BiCGSTab iteration", "-");
    }

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, double(0) * dt);

    int used_newton_iter = 0;
    int used_bicgstab_iter = 0;

    std::vector<double> his_values(1);
    his_values[0] = un[0];
    his_file->put_variable(his_O_name, nst_his, his_values);
    his_values[0] = un[1];
    his_file->put_variable(his_NO_name, nst_his, his_values);
    his_values[0] = un[2];
    his_file->put_variable(his_NO2_name, nst_his, his_values);
    his_values[0] = un[3];
    his_file->put_variable(his_O3_name, nst_his, his_values);
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        his_values[0] = used_newton_iter;
        his_file->put_variable(his_newton_name, nst_his, his_values);
        his_values[0] = used_bicgstab_iter;
        his_file->put_variable(his_bicgstab_name, nst_his, his_values);
    }

    for (int i = 0; i < 4; ++i)
    {
        up[i] = un[i];
    }

    START_TIMER(Time loop);
    std::cout << "Start time-loop" << std::endl;
    std::cout << std::setprecision(5) << double(0) * dt << "; " << double(total_time_steps - 1) * dt << std::endl;

    for (int nst = 1; nst < total_time_steps; ++nst) {
        if (fmod(nst, total_time_steps/100) == 0)
        {
            std::cout << double(nst) * dt << "; " << double(total_time_steps - 1) * dt << std::endl;
        }
        t1 = double(nst) * dt;
        if (time_integration_type == TIME_INTEGRATOR::EULER_EXPLICIT)
        {
            k[0] = k0(t0);
            explicit_euler(t0, t1, un, up, k, dt, sigma2);
        }
        if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
        {
            k[0] = k0(t0);
            runge_kutta_4(t0, t1, un, up, k, dt, sigma2);
        }
        if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
        {
            START_TIMER(Newton iteration);
            ////////////////////////////////////////////////////////////////////////////
            std::vector<double> tmp(4);

            k[0] = k0(t0);
            for (int i = 0; i < 4; ++i)
            {
                A.coeffRef(i, i) = 1.0;
                rhs[i] = 0.0;
                solution[i] = 0.0;
            }

            used_newton_iter = 0;
            used_bicgstab_iter = 0;

            START_TIMER(Newton iteration);
            for (int iter = 0; iter < iter_max; ++iter)
            {
                for (int i = 0; i < 4; ++i)
                {
                    A.coeffRef(i, i) = dtinv * 1.0;
                    rhs[i] = -(dtinv * (up[i] - un[i]));
                }
                // O
                A.coeffRef(0, 0) +=  theta * k[1];
                A.coeffRef(0, 1)  =  0.0;
                A.coeffRef(0, 2)  = -theta * k[0];
                A.coeffRef(0, 3)  =  0.0;
                //tmp[0] = rhsfev(0, up, k, sigma2);
                rhs[0] += rhsfev(0, up, k, sigma2);
                // NO
                A.coeffRef(1, 0)  = 0.0;
                A.coeffRef(1, 1) +=  theta * k[2] * up[3];
                A.coeffRef(1, 2)  = -theta * k[0];
                A.coeffRef(1, 3)  =  theta * k[2] * up[1];
                //tmp[1] = rhsfev(1, up, k, sigma2);
                rhs[1] += rhsfev(1, up, k, sigma2);
                // N02
                A.coeffRef(2, 0)  = 0.0;
                A.coeffRef(2, 1)  = -theta * k[2] * up[3];
                A.coeffRef(2, 2) +=  theta * k[0];
                A.coeffRef(2, 3)  = -theta * k[2] * up[1];
                //tmp[2] = rhsfev(2, up, k, sigma2);
                rhs[2] += rhsfev(2, up, k, sigma2);
                // O3
                A.coeffRef(3, 0)  = -theta * k[1];
                A.coeffRef(3, 1)  = theta * k[2] * up[3];
                A.coeffRef(3, 2)  = 0.0;
                A.coeffRef(3, 3) += theta * k[2] * up[1];
                //tmp[3] = rhsfev(3, up, k, sigma2);
                rhs[3] += rhsfev(3, up, k, sigma2);

                Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
                solver.compute(A);
                //solution = solver.solve(rhs);
                solver.setTolerance(eps_bicgstab);
                solution = solver.solveWithGuess(rhs, solution);  // solution contains u^{n+1,p+1} - u^{n+1,p}
                used_bicgstab_iter += int(solver.iterations());

                double dc_max = 0.0;
                for (int i = 0; i < 4; ++i)
                {
                    if (dc_max < std::abs(solution[i]))
                    {
                        dc_max = std::abs(solution[i]);
                    }
                }
                for (int i = 0; i < 4; ++i)
                {
                    up[i] += solution[i];
                }
                if (dc_max < eps_newton || iter == iter_max - 1)
                {
                    used_newton_iter += iter+1;
                    break;
                }
            }
            ////////////////////////////////////////////////////////////////////////////
            STOP_TIMER(Newton iteration);
        }
        if (std::fmod(nst, wrihis) == 0)
        {
            nst_his++;
            START_TIMER(Writing his-file);
            if (stationary)
            {
                his_file->put_time(nst_his, double(nst));
            }
            else
            {
                his_file->put_time(nst_his, double(nst) * dt);
            }
            his_values[0] = up[0];
            his_file->put_variable(his_O_name, nst_his, his_values);
            his_values[0] = up[1];
            his_file->put_variable(his_NO_name, nst_his, his_values);
            his_values[0] = up[2];
            his_file->put_variable(his_NO2_name, nst_his, his_values);
            his_values[0] = up[3];
            his_file->put_variable(his_O3_name, nst_his, his_values);
            if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
            {
                his_values[0] = used_newton_iter;
                his_file->put_variable(his_newton_name, nst_his, his_values);
                his_values[0] = used_bicgstab_iter / used_newton_iter;
                his_file->put_variable(his_bicgstab_name, nst_his, his_values);
            }
            STOP_TIMER(Writing his-file);
        }
        t0 = t1;
        for (int i = 0; i < 4; ++i)
        {
            un[i] = up[i];
        }
    }
    STOP_TIMER(Time loop);
    STOP_TIMER(main);
    PRINT_TIMER(timing_filename.data());

    (void)his_file->close();

    if (time_integration_type == TIME_INTEGRATOR::EULER_EXPLICIT)
    {
        std::cout << "-------------------------------" << std::endl;
        std::cout << "Time integrator: EULER_EXPLICIT" << std::endl;
        std::cout << "-------------------------------" << std::endl;
    }
    if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
    {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Time integrator: RUNGE_KUTTA_4" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        std::cout << "-------------------------------" << std::endl;
        std::cout << "Time integrator: FULLY_IMPLICIT" << std::endl;
        std::cout << "Total number of Newton iterations: " << used_newton_iter/total_time_steps << std::endl;
        std::cout << "Total number of linear iterations: " << used_bicgstab_iter/total_time_steps << std::endl;
        std::cout << "-------------------------------" << std::endl;
    }
    std::cout << "Press Enter to finish";
    std::cin.ignore();

    return 0;
}
//------------------------------------------------------------------------------
double k0(double t0)
{
    // compute time dependent k-coefficient
    double th = fmod(t0 / 3600.0, 24.0);
    double kres = 1.0e-40; // at night
    if (th >= 4.0 && th <= 20.0) 
    {
        // daylight period
        kres = 1.0e-5 * exp(7.0 * pow(sin(M_PI / 16 * (th - 4)), 0.2));
    }
    return kres;
}
//------------------------------------------------------------------------------
void explicit_euler(double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k, double dt, double sigma2) {
    double at;
    std::vector<double> au(4), k1(4), k2(4), k3(4), k4(4);

    at = t0;
    for (int i = 0; i < 4; ++i) {
        au[i] = u0[i];
    }
    for (int i = 0; i < 4; ++i) {
        k1[i] = rhsfev(i, au, k, sigma2); // Assuming rhsfev is defined elsewhere
    }

    for (int i = 0; i < 4; ++i) {
        u1[i] = u0[i] + dt * k1[i];
    }
    t1 = t0 + dt;
}
//------------------------------------------------------------------------------
void runge_kutta_4(double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k, double dt, double sigma2) {
    double at, nul1, nul2;
    std::vector<double> au(4), k1(4), k2(4), k3(4), k4(4);

    at = t0;
    for (int i = 0; i < 4; ++i) {
        au[i] = u0[i];
    }
    for (int i = 0; i < 4; ++i) {
        k1[i] = rhsfev(i, au, k, sigma2);
    }
    nul1 = k1[0] + k1[2] + k1[3];
    nul2 = k1[1] + k1[2] - sigma2;

    at = t0 + 0.5 * dt;
    for (int i = 0; i < 4; ++i) {
        au[i] = u0[i] + 0.5 * dt * k1[i];
    }
    for (int i = 0; i < 4; ++i) {
        k2[i] = rhsfev(i, au, k, sigma2);
    }
    nul1 = k2[0] + k2[2] + k2[3];
    nul2 = k2[1] + k2[2] - sigma2;

    at = t0 + 0.5 * dt;
    for (int i = 0; i < 4; ++i) {
        au[i] = u0[i] + 0.5 * dt * k2[i];
    }
    for (int i = 0; i < 4; ++i) {
        k3[i] = rhsfev(i, au, k, sigma2);
    }
    nul1 = k3[0] + k3[2] + k3[3];
    nul2 = k3[1] + k3[2] - sigma2;

    at = t0 + dt;
    for (int i = 0; i < 4; ++i) {
        au[i] = u0[i] + dt * k3[i];
    }
    for (int i = 0; i < 4; ++i) {
        k4[i] = rhsfev(i, au, k, sigma2);
    }
    nul1 = k4[0] + k4[2] + k4[3];
    nul2 = k4[1] + k4[2] - sigma2;

    for (int i = 0; i < 4; ++i) {
        u1[i] = u0[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    t1 = t0 + dt;
}
//------------------------------------------------------------------------------
double rhsfev(int ispecie, std::vector<double>& u, std::vector<double>& k, double sigma2) {
    double result = 0.0;

    if (ispecie == 0) {
        result = k[0] * u[2] - k[1] * u[0];
    }
    else if (ispecie == 1) {
        result = k[0] * u[2] - k[2] * u[1] * u[3] + sigma2;
    }
    else if (ispecie == 2) {
        result = k[2] * u[1] * u[3] - k[0] * u[2];
    }
    else if (ispecie == 3) {
        result = k[1] * u[0] - k[2] * u[1] * u[3];
    }

    return result;
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
