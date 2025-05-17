#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <algorithm>
#include <filesystem>

//
//
// Brusselator problem:
//     https://mate.unipv.it/~boffi/teaching/download/Brusselator.pdf
//
//  dx/dt = 1 − (b + 1)x + ax2
//  dy/dt = bx − ax2
//

#include <Eigen/IterativeLinearSolvers>
#include <toml.h>

#include "cfts.h"
#include "perf_timer.h"

double rhsfev(int ispecie, const std::vector<double>& u, const std::vector<double>& k);
void runge_kutta_4(double t0, double& t1, const std::vector<double>& u0, std::vector<double>& u1, const std::vector<double>& k, double dt);
void GetArguments(long argc, char** argv, std::string* file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

enum class TIME_INTEGRATOR
{
    EULER_EXPLICIT,
    RUNGE_KUTTA_4,
    FULLY_IMPLICIT,
};

int main(int argc, char* argv[]) {
    int status = -1;
    std::string toml_file_name("--- not-defined ---");
    
    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path start_dir;
    std::filesystem::path output_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();
    start_dir = std::filesystem::current_path();

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
            //std::cin.ignore();
            exit(1);
        }
        tbl = toml::parse_file(toml_file_name);
        std::cout << tbl << "\n";
        output_dir = start_dir;
        output_dir += "/output/";
        std::filesystem::create_directory(output_dir);
    }
    else
    {
        std::cout << "No \'toml\' file is read." << std::endl;
        output_dir = ".";
    }

    std::cout << "----------------------------" << std::endl;
    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Start directory     : " << start_dir << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;
    std::cout << "----------------------------" << std::endl;

    std::string out_file;
    std::stringstream ss;
    ss << "brusselator";
    out_file = output_dir.string() + ss.str();

    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Two way chemical reaction, <no time integrator specified>");

    Eigen::SparseMatrix<double> A(2, 2);
    Eigen::VectorXd rhs(2);                // RHS vector [u1, u2, u3, u4]^n
    Eigen::VectorXd solution(2);           // solution vector [u1, u2, u3, u4]^{n}

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

    log_file << "=======================================================" << std::endl;

    std::vector<double> un(2, 0.);
    std::vector<double> dundt(2, 0.);
    std::vector<double> up(2, 0.);
    std::vector<double> dupdt(2, 0.);
    std::vector<double> u_ini(2, 0.);
    std::vector<double> u_ana(2, 0.);
    std::vector<double> k(2, 0.);

    if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
    {
        model_title = "Brusselator, Runge-Kutta 4 (RK4)";
    }
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        model_title = "Brusselator, Delta formulation";
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

    un[0] = 0.0;
    un[1] = 0.0;
    u_ini[0] = un[0];
    u_ini[1] = un[1];
    u_ana[0] = un[0];
    u_ana[1] = un[1];
    k[0] = k1;  // i.e. the 'a' in the equations
    k[1] = k2;  // i.e. the 'b' in the equations

    double t0 = tstart * dt;
    double t1 = t0;

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

    std::string his_u1_name("u1");
    std::string his_u2_name("u2");
    std::string his_newton_name("newton");
    std::string his_bicgstab_name("bicgstab");
    his_file->add_variable(his_u1_name, "", "U1", "-");
    his_file->add_variable(his_u2_name, "", "U2", "-");
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
    his_file->put_variable(his_u1_name, nst_his, his_values);
    his_values[0] = un[1];
    his_file->put_variable(his_u2_name, nst_his, his_values);
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        his_values[0] = used_newton_iter;
        his_file->put_variable(his_newton_name, nst_his, his_values);
        his_values[0] = used_bicgstab_iter;
        his_file->put_variable(his_bicgstab_name, nst_his, his_values);
    }

    for (int i = 0; i < 2; ++i)
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

        u_ana[0] = 0.0;
        u_ana[1] = 0.0;

        if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
        {
            runge_kutta_4(t0, t1, un, up, k, dt);
        }

        if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
        {
            START_TIMER(Newton iteration);
            ////////////////////////////////////////////////////////////////////////////
            double dc_max = 0.0;
            std::vector<double> tmp(2);

            for (int i = 0; i < 2; ++i)
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
                for (int i = 0; i < 2; ++i)
                {
                    A.coeffRef(i, i) = dtinv * 1.0;
                    rhs[i] = -(dtinv * (up[i] - un[i]));
                }
                // u0
                A.coeffRef(0, 0) += -theta * (- k[1] - 1.0 + 2.0 * k[0] * up[0] * up[1]);
                A.coeffRef(0, 1)  = -theta * k[0] * up[0] * up[0];
                rhs[0] += rhsfev(0, up, k);
                // u1
                A.coeffRef(1, 0) = -theta * (k[1] - 2.0 * up[0] * up[1]);
                A.coeffRef(1, 1) += -theta * (-k[0] * up[0] * up[0]);
                rhs[1] += rhsfev(1, up, k);

                Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
                solver.compute(A);
                //solution = solver.solve(rhs);
                solver.setTolerance(eps_bicgstab);
                solution = solver.solveWithGuess(rhs, solution);  // solution contains u^{n+1,p+1} - u^{n+1,p}
                used_bicgstab_iter += int(solver.iterations());

                dc_max = 0.0;
                for (int i = 0; i < 2; ++i)
                {
                    if (dc_max < std::abs(solution[i]))
                    {
                        dc_max = std::abs(solution[i]);
                    }
                }
                for (int i = 0; i < 2; ++i)
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
            his_file->put_variable(his_u1_name, nst_his, his_values);
            his_values[0] = up[1];
            his_file->put_variable(his_u2_name, nst_his, his_values);
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
        for (int i = 0; i < 2; ++i)
        {
            un[i] = up[i];
        }
    }
    STOP_TIMER(Time loop);
    STOP_TIMER(Main);
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
        std::cout << "Total number of Newton iterations: " << used_newton_iter / total_time_steps << std::endl;
        std::cout << "Total number of linear iterations: " << used_bicgstab_iter / total_time_steps << std::endl;
        std::cout << "-------------------------------" << std::endl;
    }
    //std::cout << "Press Enter to finish";
    //std::cin.ignore();

    return 0;
}
//------------------------------------------------------------------------------
void runge_kutta_4(double t0, double& t1, const std::vector<double>& u0, std::vector<double>& u1, const std::vector<double>& k, double dt) {
    double at = t0;
    std::vector<double> au = u0;
    std::vector<double> k1(2), k2(2), k3(2), k4(2);

    for (int i = 0; i < 2; ++i) {
        k1[i] = rhsfev(i, au, k);
    }
    at = t0 + 0.5 * dt;
    for (int i = 0; i < 2; ++i) {
        au[i] = u0[i] + 0.5 * dt * k1[i];
    }
    for (int i = 0; i < 2; ++i) {
        k2[i] = rhsfev(i, au, k);
    }
    at = t0 + 0.5 * dt;
    for (int i = 0; i < 2; ++i) {
        au[i] = u0[i] + 0.5 * dt * k2[i];
    }
    for (int i = 0; i < 2; ++i) {
        k3[i] = rhsfev(i, au, k);
    }
    at = t0 + dt;
    for (int i = 0; i < 2; ++i) {
        au[i] = u0[i] + dt * k3[i];
    }
    for (int i = 0; i < 2; ++i) {
        k4[i] = rhsfev(i, au, k);
    }
    for (int i = 0; i < 2; ++i) {
        u1[i] = u0[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    t1 = t0 + dt;
}
//------------------------------------------------------------------------------
double rhsfev(int ispecie, const std::vector<double>& u, const std::vector<double>& k) {
    double result = 0.0;

    if (ispecie == 0) {
        result = 1. - (k[1] + 1.0) * u[0] + k[0] * u[0] * u[0] * u[1];
    }
    else if (ispecie == 1) {
        result = k[1] * u[0] - k[0] * u[0] * u[0] * u[1];
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
