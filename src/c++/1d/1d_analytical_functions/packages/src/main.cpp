//
//
// Programmer: J. Mooiman
// Date      : 2025-03-01
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
#include <thread>

#include <toml.h>

void GetArguments(long argc, char** argv, std::string* file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

#include "ugrid1d.h"
#include "adv_diff_init_concentration.h"
#include "perf_timer.h"
#include "definition_map_file.h"

int main(int argc, char* argv[])
{
    bool stationary = false;
    std::string toml_file_name("--- not-defined ---");

    int status = -1;

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
        status = 0;
    }
    else
    {
        std::cout << "No \'toml\' file is read." << std::endl;
        current_dir = ".";
        output_dir = ".";
    }

    START_TIMERN(main);
    START_TIMER(Simulation initialization);
    SHAPE_CONC shape_conc = SHAPE_CONC::NONE;

    // Domain
    tbl_chp = *tbl["Domain"].as_table();
    double Lx = tbl_chp["Lx"].value_or(double(1.0));

    // Time
    tbl_chp = *tbl["Time"].as_table();
    double tstart = tbl_chp["tstart"].value_or(double(0.0));
    double tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.));

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    std::string shape_type = tbl_chp["shape"].value_or("None");
    if (shape_type == "Constant") { shape_conc = SHAPE_CONC::Constant; }
    else if (shape_type == "Envelope") { shape_conc = SHAPE_CONC::Envelope; }
    else if (shape_type == "Boundary_layer") { shape_conc = SHAPE_CONC::Boundary_layer; }
    else { shape_conc = SHAPE_CONC::NONE; }

    //Physics
    tbl_chp = *tbl["Physics"].as_table();
    std::vector<double> u;
    status = get_toml_array(tbl_chp, "u", u);
    double diffusion = tbl_chp["diffusion"].value_or(double(0.0));

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    double dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary
    if (dt == 0.0) { stationary = true; }
    double dx = tbl_chp["dx"].value_or(double(0.01));

    // Numerics
    tbl_chp = *tbl["Output"].as_table();
    double dt_map = tbl_chp["dt_map"].value_or(double(0.01));


    std::string out_file;
    std::stringstream ss;

    std::string logging = tbl["Logging"].value_or("None");

    std::string model_title("not defined");
    if (shape_conc == SHAPE_CONC::Envelope) {
        model_title = "Analytic wave package";
        ss << "analytic_wave_package";
    }
    if (shape_conc == SHAPE_CONC::Boundary_layer) {
        model_title = "Analytic boundary layer";
        ss << "analytic_boundary_layer";

    }
    out_file = output_dir.string() + ss.str();
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

    //  select 
    //      Constant c, u>0
    //      Sine c(0,t)= sine at west boundary, uniform u>0
    //      Envelope  uniform u>0, sine within cosine envelope

    double dxinv = 1. / dx;
    int nx = int(Lx * dxinv) + 1 + 2; // nr nodes; including 2 virtual points

    int total_time_steps = (int) u.size();

    log_file << "=== Used input variables ==============================" << std::endl;
    log_file << "[Domain]" << std::endl;
    log_file << "Lx = " << Lx << std::endl;

    log_file << std::endl << "[Time]" << std::endl;
    log_file << "tstart = " << tstart << std::endl;
    log_file << "tstop = " << tstop << std::endl;

    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "shape = " << shape_type << std::endl;

    log_file << std::endl << "[Physics]" << std::endl;
    for (int i = 0; i < u.size() - 1; ++i) { log_file << u[i] << ", "; }
    log_file << u[u.size() - 1] << std::endl;
    log_file << "diffusion = " << diffusion << std::endl;

    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "dt  = " << dt << std::endl;
    log_file << "dx  = " << dx << std::endl;

    log_file << "-------------------------------------------------------" << std::endl;
    log_file << "Nodes  : " << nx << std::endl;
    log_file << "Volumes: " << (nx - 1) << std::endl;
    for (int i = 0; i < u.size() - 1; ++i)
    {
        log_file << "CFL    : " << u[i] * dt / dx << ", ";
    }
    log_file << "CFL    : " << u[u.size()-1] * dt / dx << std::endl;

    STOP_TIMER(Writing log-file);  // but two write statements are not timed
    if (status != 0) {
        log_file.close();
        exit(1);
    }
    log_file << "=======================================================" << std::endl;

    std::vector<double> x(nx, 0.);  // x-coordinate
    std::vector<double> y(nx, 0.);  // y-coordinate
    std::vector<double> conc(nx, 0.0);  //

    //initialize x-coordinate
    for (int i = 0; i < nx; i++)
    {
        x[i] = double(i-1) * dx;
    }

    // initial concentration
    double time = tstart + dt * double(0);
    (void)adv_diff_init_concentration(time, u[0], diffusion, x, Lx, shape_conc, conc);

    ////////////////////////////////////////////////////////////////////////////
    // Define map file 
    std::cout << "    Create map-file" << std::endl;
    std::string nc_mapfilename(map_filename);
    UGRID1D* map_file = new UGRID1D();

    std::vector<std::string> map_names;
    map_names.push_back("c_1d");
    map_file = create_map_file(nc_mapfilename, model_title, x, map_names);

    // Put data on map file
    int nst_map = 0;
    map_file->put_time(nst_map, time);
    map_file->put_time_variable(map_names[0], nst_map, conc);

    ////////////////////////////////////////////////////////////////////////////
    // start time loop
    std::cout << "Start loop" << std::endl;

    STOP_TIMER(Simulation initialization);
    for (int nst = 1; nst < total_time_steps; ++nst)
    {
        // Map-files
        time = tstart + dt * double(nst);

        (void)adv_diff_init_concentration(time, u[nst], diffusion, x, Lx, shape_conc, conc);
        // Put data on time map file
        START_TIMER(Writing map-file);
        nst_map++;
        map_file->put_time(nst_map, double(nst));
        map_file->put_time_variable(map_names[0], nst_map, conc);
        STOP_TIMER(Writing map-file);
    } // End of the time loop

    log_file.close();
    (void)map_file->close();
    STOP_TIMER(Writing log-file);

    STOP_TIMER(main);
    PRINT_TIMER(timing_filename.data());

    std::chrono::duration<int, std::milli> timespan(1000);
    std::this_thread::sleep_for(timespan);
    return 0;
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

