#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <algorithm>
#include <filesystem>
#include <thread>

//
// Two way chemical reaction:
//     Hundsdorfer and Verwer 2003, pg 26 and pg 35
//
//  du_1/dt = -k1 u_1 + k2 u_2
//  du_2/dt =  k1 u_1 - k2 u_2
//

#include <Eigen/IterativeLinearSolvers>
#include <toml.h>

#include "cfts.h"
#include "compile_date_and_time.h"
#include "data_input_struct.h"
#include "main_version.h"
#include "perf_timer.h"
#include "read_input_toml_file.h"

double rhsfev(int ispecie, const std::vector<double>& u, const std::vector<double>& k);
void runge_kutta_4(double t0, double& t1, const std::vector<double>& u0, std::vector<double>& u1, const std::vector<double>& k, double dt);
void GetArguments(long argc, char** argv, std::filesystem::path & file_name);
int write_used_input(struct _data_input data, std::ofstream & log_file);


int main(int argc, char* argv[]) {
    int status = -1;
    std::filesystem::path toml_file_name("---not-defined---");

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
        std::cout << "Git commit time     : " << getgitdatestring_main() << std::endl;
        std::cout << "Git commit hash     : " << getgitbuildstring_main() << std::endl;
        std::cout << "Git branch          : " << getgitbranchstring_main() << std::endl;
        std::cout << "Git repository URL  : " << getgiturlstring_main() << std::endl;
        std::cout << std::endl;
        std::cout << "usage: two_way_chemical_reaction.exe --toml <input_file>" << std::endl;
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
    std::cout << "Git commit time     : " << getgitdatestring_main() << std::endl;
    std::cout << "Git commit hash     : " << getgitbuildstring_main() << std::endl;
    std::cout << "Git branch          : " << getgitbranchstring_main() << std::endl;
    std::cout << "Git repository URL  : " << getgiturlstring_main() << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Start directory     : " << start_dir << std::endl;
    std::cout << "Input directory     : " << input_dir << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;

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

    ss << "two_way_chemical_reaction";
    out_file = output_dir.string() + ss.str();

    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("Two way chemical reaction, <no time integrator specified>");

    std::ofstream log_file;
    log_file.open(log_filename);

    // Initial
    double u0_init = input_data.initial.u0_ini;
    double u1_init = input_data.initial.u1_ini;

    // Numerics
    double dt = input_data.numerics.dt;
    std::string time_integration = input_data.numerics.time_integration;
    int iter_max = input_data.numerics.iter_max;
    double theta = input_data.numerics.theta;  // Implicitness factor (0.5 <= theta <= 1.0)
    double eps_newton = input_data.numerics.eps_newton;  // stop criterium for Newton iteration
    double eps_bicgstab = input_data.numerics.eps_bicgstab;  // stop criterium for BiCGStab iteration

    // Physics
    double k1 = input_data.physics.k1;
    double k2 = input_data.physics.k2;

    // Output
    double dt_his = input_data.output.dt_his;
    double dt_screen = input_data.output.dt_screen;

    int nst_to_his = std::max(1, int(dt_his/dt));
    int nst_to_screen = std::max(1, int(dt_screen/dt));

    // Time
    double tstart = input_data.time.tstart;
    double tstop = input_data.time.tstop;

    std::vector<double> un(2, 0.);
    std::vector<double> dundt(2, 0.);
    std::vector<double> up(2, 0.);
    std::vector<double> dupdt(2, 0.);
    std::vector<double> u_ini(2, 0.);
    std::vector<double> u_ana(2, 0.);
    std::vector<double> k(2, 0.);

    log_file << "=== Used input variables ==============================" << std::endl;
    status = write_used_input(input_data, log_file);
    log_file << "=======================================================" << std::endl;

    if (time_integration == "runge_kutta_4")
    {
        model_title = "Two way chemical reaction, Runge-Kutta 4 (RK4)";
    }
    else if (time_integration == "delta_formulation")
    {
        model_title = "Two way chemical reaction, Delta formulation";
    }
    else
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Integration type \'" << time_integration << "\' not supported." << std::endl;
        log_file << "Integration type \'" << time_integration << "\' not supported." << std::endl;
        //std::cout << "Press Enter to finish";
        //std::cin.ignore();
        std::chrono::duration<int, std::milli> timespan(3000);
        std::this_thread::sleep_for(timespan);
        exit(1);
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

    un[0] = u0_init;
    un[1] = u1_init;
    u_ini[0] = un[0];
    u_ini[1] = un[1];
    u_ana[0] = un[0];
    u_ana[1] = un[1];
    k[0] = k1;
    k[1] = k2;  

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
    if (time_integration == "delta_formulation")
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
    if (time_integration == "delta_formulation")
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

    Eigen::SparseMatrix<double> A(2, 2);
    Eigen::VectorXd rhs(2);                // RHS vector [u1, u2]^n
    Eigen::VectorXd solution(2);           // solution vector [u1, u2]^{n}

    START_TIMER(Time loop);
    std::cout << "Start time-loop" << std::endl;
    std::cout << std::setprecision(5) << double(0) * dt << "; " << double(total_time_steps - 1) * dt << std::endl;

    for (int nst = 1; nst < total_time_steps; ++nst) {
        if (nst % nst_to_screen == 0)
        {
            std::cout << double(nst) * dt << "; " << double(total_time_steps - 1) * dt << std::endl;
        }
        t1 = double(nst) * dt;

        u_ana[0] = 0.0;
        u_ana[1] = 0.0;

        if (time_integration == "runge_kutta_4")
        {
            runge_kutta_4(t0, t1, un, up, k, dt);
        }

        if (time_integration == "delta_formulation")
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
                A.coeffRef(0, 0) += -theta * -k[0];
                A.coeffRef(0, 1)  = -theta * k[1];
                rhs[0] += rhsfev(0, up, k);
                // u1
                A.coeffRef(1, 0) = -theta * k[0];
                A.coeffRef(1, 1) += -theta * -k[1];
                rhs[1] += rhsfev(1, up, k);

                Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
                solver.compute(A);
                solver.setTolerance(eps_bicgstab);
                solution = solver.solve(rhs);
                //solution = solver.solveWithGuess(rhs, solution);  // solution contains u^{n+1,p+1} - u^{n+1,p}
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
        if (nst % nst_to_his == 0)
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
            if (time_integration == "delta_formulation")
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

    if (time_integration == "euler_explicit")
    {
        std::cout << "-------------------------------" << std::endl;
        std::cout << "Time integrator: EULER_EXPLICIT" << std::endl;
        std::cout << "-------------------------------" << std::endl;
    }
    if (time_integration == "runge_kutta_4")
    {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Time integrator: RUNGE_KUTTA_4" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    if (time_integration == "delta_formulation")
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
        result = -k[0] * u[0] + k[1] * u[1];
    }
    else if (ispecie == 1) {
        result = k[0] * u[0] - k[1] * u[1];
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
