#define _USE_MATH_DEFINES
#include <ios>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <filesystem>

#include "Eigen/IterativeLinearSolvers"
#include "include/netcdf.h"

#include "cfts.h"
#include "perf_timer.h"

enum class TIME_INTEGRATOR
{
    EULER_EXPLICIT,
    RUNGE_KUTTA_4,
    FULLY_IMPLICIT,
};
double k0(double t0);
void explicit_euler(double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k, 
    double dt, double sigma2);
void runge_kutta_4 (double t0, double& t1, std::vector<double>& u0, std::vector<double>& u1, std::vector<double>& k,
    double dt, double sigma2);
double rhsfev(int, std::vector<double>&, std::vector<double>&, double);

int main(int argc, char* argv[]) {
    int status;
    double t0, t1, dt, days, dtinv;
    std::vector<double> un(4, 0.);
    std::vector<double> up(4, 0.);
    std::vector<double> k(3, 0.);
    double tstart, tstop;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path current_dir;
    std::filesystem::path output_dir;

    if (argc == 1)
    {
        exec_file = argv[0];
        exec_dir = exec_file.parent_path();
    }
    else
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "No support of commandline arguments" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
        exit(1);
    }

    current_dir = ".";
    output_dir = current_dir.string() + "/output/";
    std::filesystem::create_directory(output_dir);

    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Current directory   : " << std::filesystem::relative(current_dir) << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;
    std::cout << std::endl;


    std::string filename;

    std::string out_file;
    std::stringstream ss;
    ss << "air_pollution";
    out_file = output_dir.string() + ss.str();
    std::string timing_filename(out_file + "_timing.log");
    std::string his_filename(out_file + "_his.nc");
    std::string model_title("Advection-diffusion-reaction, <no time integrator specified>");

    Eigen::SparseMatrix<double> A(4, 4);
    Eigen::VectorXd rhs(4);                // RHS vector [u1, u2, u3, u4]^n
    Eigen::VectorXd solution(4);           // solution vector [u1, u2, u3, u4]^{n}

    enum TIME_INTEGRATOR time_integration_type = TIME_INTEGRATOR::RUNGE_KUTTA_4;  // set RK4 as default
    time_integration_type = TIME_INTEGRATOR::EULER_EXPLICIT;
    time_integration_type = TIME_INTEGRATOR::RUNGE_KUTTA_4;
    time_integration_type = TIME_INTEGRATOR::FULLY_IMPLICIT;

    if (time_integration_type == TIME_INTEGRATOR::EULER_EXPLICIT)
    {
        model_title = "Advection-diffusion-reaction, Euler explicit";
    }
    if (time_integration_type == TIME_INTEGRATOR::RUNGE_KUTTA_4)
    {
        model_title = "Advection-diffusion-reaction, Runge-Kutta 4 (RK4)";
    }
    if (time_integration_type == TIME_INTEGRATOR::FULLY_IMPLICIT)
    {
        model_title = "Advection-diffusion-reaction, Delta formulation";
    }

    START_TIMERN(main);

    dt = 0.5;
    dtinv = 1. / dt;
    days = 4.0;
    tstart = 0.0; // start om 0 uur van de eerste dag
    tstop = days * 24.0 * 3600.0; // [sec]
    // tstop = 600.0; // [sec]
    int total_time_steps = int((tstop - tstart) / dt) + 1;
    int wrihis = int(std::max(1., 1.0 / dt));      // write interval to his-file (every delta t)

    un[0] = 0.0;     // ! 0.d0 ! O Zuurstof-molekuul (Oxygen-molecule)
    un[1] = 2.0e-1;  // 1.3d08 ! NO Stikstof-oxide
    un[2] = 2.0e-3;  // 5.0d11 ! NO2 Stikstof-di-oxide
    un[3] = 2.0e-1;  // 8.0d11 ! O3 = Ozon
    double sigma2 = 1.0e-7;  // 1.0d06 
    k[0] = 0.0;     // Time dependent
    k[1] = 2.0e-2;  // 1.0d05
    k[2] = 1.0e-3;  // 1.0d-16

    t0 = tstart * dt;

    //data_file << std::setw(21) << "time"
    //        << std::setw(21) << "1:O" 
    //        << std::setw(21) << "2:NO" 
    //        << std::setw(21) << "3:NO2" 
    //        << std::setw(21) << "4:O3" << std::endl;
    //data_file << std::scientific << std::setprecision(13) << t1;
    //for (int i = 0; i < 4; ++i) {
    //    data_file << std::setw(21) << u1[i];
    //}
    //data_file << std::endl;

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
    his_file->add_variable(his_O_name, "", "atomic oxygen", "-");
    his_file->add_variable(his_NO_name, "", "nitrogen oxide", "-");
    his_file->add_variable(his_NO2_name, "", "nitrogen dioxide", "-");
    his_file->add_variable(his_O3_name, "", "ozon", "-");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, double(0) * dt);

    std::vector<double> his_values(1);
    his_values[0] = un[0];
    his_file->put_variable(his_O_name, nst_his, his_values);
    his_values[0] = un[1];
    his_file->put_variable(his_NO_name, nst_his, his_values);
    his_values[0] = un[2];
    his_file->put_variable(his_NO2_name, nst_his, his_values);
    his_values[0] = un[3];
    his_file->put_variable(his_O3_name, nst_his, his_values);

    int used_newton_iter = 0;
    int used_bicgstab_iter = 0;
    for (int i = 0; i < 4; ++i)
    {
        up[i] = un[i];
    }
    START_TIMER(time - loop);
    std::cout << "Start time-loop" << std::endl;
    std::cout << std::setprecision(5) << double(0) * dt << "; " << double(total_time_steps - 1) * dt << std::endl;
    for (int nst = 1; nst < total_time_steps; ++nst) {
        if (fmod(nst, 10000) == 0)
        {
            std::cout << double(nst) * dt << "; " << double(total_time_steps-1) * dt << std::endl;
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
            START_TIMER(Newton iteration initialization);
            ////////////////////////////////////////////////////////////////////////////
            double theta = 0.501;
            double eps_newton = 1e-12;
            double eps_bicgstab = 1e-12;
            double dc_max = 0.0;
            std::vector<double> tmp(4);

            k[0] = k0(t0);

            for (int i = 0; i < 4; ++i)
            {
                A.coeffRef(i, i) = 1.0;
                rhs[i] = 0.0;
                solution[i] = 0.0;
            }
            STOP_TIMER(Newton iteration initialization);

            START_TIMER(Newton iteration);
            int iter_max = 100;
            for (int iter = 0; iter < iter_max; ++iter)
            {
                for (int i = 0; i < 4; ++i)
                {
                    A.coeffRef(i, i) = dtinv * 1.0;
                    //tmp[i] = -(dtinv * (up[i] - un[i]));
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
                dc_max = 0.0;
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
            STOP_TIMER(Newton iteration);
            ////////////////////////////////////////////////////////////////////////////
        }
        if (std::fmod(nst, wrihis) == 0)
        {
            nst_his++;
            START_TIMER(Writing his - file);
            his_file->put_time(nst_his, t1);
            his_values[0] = up[0];
            his_file->put_variable(his_O_name, nst_his, his_values);
            his_values[0] = up[1];
            his_file->put_variable(his_NO_name, nst_his, his_values);
            his_values[0] = up[2];
            his_file->put_variable(his_NO2_name, nst_his, his_values);
            his_values[0] = up[3];
            his_file->put_variable(his_O3_name, nst_his, his_values);
            STOP_TIMER(Writing his - file);
        }
        t0 = t1;
        for (int i = 0; i < 4; ++i)
        {
            un[i] = up[i];
        }
    }
    (void)his_file->close();
    STOP_TIMER(time - loop);

    STOP_TIMER(main);
    PRINT_TIMER(timing_filename.data());
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


