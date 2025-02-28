//
// Programmer: J. Mooiman
// Date      : 2023-07-07
//
//

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <filesystem>
#include <toml.h>

// for bicgstab  solver
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "cfts.h"
#include "ugrid2d.h"
#include "perf_timer.h"

void GetArguments(long argc, char** argv, std::string* file_name);
int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

int p_index(int i, int j, int nx);
double c_scv(double, double, double, double);
double dcdx_scvf1(double, double, double, double);
double dcdy_scvf1(double, double, double, double);
double dcdx_scvf2(double, double, double, double);
double dcdy_scvf2(double, double, double, double);

// Solve the linear wave equation
// Continuity equation: d(h)/dt + d(q)/dx = 0
// Momentum equation: d(q)/dt + 1/2  d(h^2)/dx = 0
// Momentum equation: d(r)/dt + 1/2  d(h^2)/dy = 0

int main(int argc, char *argv[])
{
    START_TIMERN(Main);

    bool stationary = false;
    std::string toml_file_name("---not-defined---");
    int status = -1;

    int BC_NORTH = 0;
    int BC_EAST = 1;
    int BC_SOUTH = 2;
    int BC_WEST = 3;

    std::filesystem::path exec_file;
    std::filesystem::path exec_dir;
    std::filesystem::path current_dir;
    std::filesystem::path output_dir;

    exec_file = argv[0];
    exec_dir = exec_file.parent_path();

    current_dir = "./../..";
    output_dir = exec_dir.string() + "../../output/";
    std::filesystem::create_directory(output_dir);

    std::cout << "Executable directory: " << exec_dir << std::endl;
    std::cout << "Current directory   : " << std::filesystem::absolute(current_dir) << std::endl;
    std::cout << "Output directory    : " << output_dir << std::endl;
    std::cout << std::endl;

    toml::table tbl;
    if (argc == 3 )
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

    std::cout << "Start" << std::endl;

    START_TIMER(Initialization);
    
    std::string out_file;
    std::stringstream ss;
    ss << "2d_wave";
    out_file = output_dir.string() + ss.str();
    std::string his_filename(out_file + "_his.nc");
    std::string log_filename(out_file + ".log");
    std::string map_filename(out_file + "_map.nc");
    std::string timing_filename(out_file + "_timing.log");
    std::string model_title("SMART Numerics: Linear Wave, BiCGSTab");

    std::ofstream log_file;
    log_file.open(log_filename);
    std::cout << "=== Input file ========================================" << std::endl;
    std::cout << std::filesystem::absolute(exec_dir) << std::endl;  
    std::cout << std::filesystem::absolute(toml_file_name) << std::endl;  
    std::cout << "=======================================================" << std::endl;
    log_file << "=== Input file ========================================" << std::endl;
    log_file << std::filesystem::absolute(toml_file_name) << std::endl;  
    log_file << "=== Copy of the input file ============================" << std::endl;
    log_file << tbl << std::endl;  // Write input TOML file to log_file
    log_file << "=======================================================" << std::endl;
    log_file << std::endl;
    log_file << "=== Used input variables ==============================" << std::endl;

    // Time
    log_file << std::endl << "[Time]" << std::endl;
    double tstart = tbl["Time"]["tstart"].value_or(double(0.0));
    log_file << "tstart = " << tstart << std::endl;
    double tstop = tbl["Time"]["tstop"].value_or(double(1800.));
    log_file << "tstop = " << tstop << std::endl;
    double dt = tbl["Time"]["dt"].value_or(double(5.0));  // default stationary
    if (dt == 0.0) { stationary = true;  }
    log_file << "dt = " << dt << std::endl;
    
    // Initial
    log_file << std::endl << "[Initial]" << std::endl;
    std::vector<std::string> ini_vars; // { "velocity", "zeta", "viscosity", "zeta_GaussHump" };  // get the element as an array
    auto vars = tbl["Initial"];
    status = get_toml_array(*vars.as_table(), "ini_vars", ini_vars);
    double gauss_mu = vars["gaussian_mu"].value_or(double(0.0));
    log_file << "Gauss_mu = " << gauss_mu << std::endl;
    double gauss_sigma = vars["gaussian_sigma"].value_or(double(1.0));
    log_file << "Gauss_sigma = " << gauss_sigma << std::endl;
    double a0 = tbl["Initial"]["a0"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary

    // Domain
    log_file << std::endl << "[Domain]" << std::endl;
    double Lx = tbl["Domain"]["Lx"].value_or(double(6000.0));
    log_file << "Lx = " << Lx << std::endl;
    double Ly = tbl["Domain"]["Ly"].value_or(double(0.0));
    if (Ly == 0.0) { Ly = Lx; }
    log_file << "Ly = " << Ly << std::endl;

    //Physics
    log_file << std::endl << "[Physics]" << std::endl;
    bool do_continuity = tbl["Physics"]["do_continuity"].value_or(bool(true));  // default, continuity
    log_file << "do_continuity = " << do_continuity << std::endl;
    bool do_q_equation = tbl["Physics"]["do_q_equation"].value_or(bool(true));  // default, q_equation
    log_file << "do_q_equation = " << do_q_equation << std::endl;
    bool do_r_equation = tbl["Physics"]["do_r_equation"].value_or(bool(true));  // default, r_equation
    log_file << "do_r_equation = " << do_r_equation << std::endl;
    if (do_q_equation && do_r_equation)
    {
        model_title = "SMART Numerics 2D: Linear wave equation, BiCGSTab";
    }
    else if (do_q_equation)
    {
        model_title = "SMART Numerics 2D: Only q-equation, BiCGSTab";
    }
    else if (do_r_equation)
    {
        model_title = "SMART Numerics 2D: Only r-equation, BiCGSTab";
    }
    else
    {
        model_title = "SMART Numerics: No q- and no r-equation (=> no waves computed), BiCGSTab";
    }
    bool do_convection = tbl["Physics"]["do_convection"].value_or(bool(false));  // default, no convection
    log_file << "do_convection = " << do_convection << std::endl;
    bool do_viscosity = tbl["Physics"]["do_viscosity"].value_or(bool(false));  // default, no convection
    double viscosity = tbl["Physics"]["viscosity"].value_or(double(0.0));
    if (viscosity == 0.0) { do_viscosity = false;  }
    log_file << "do_viscosity = " << do_viscosity << std::endl;
    if (do_viscosity) { log_file << "viscosity = " << viscosity << std::endl; }

    // boundary conditions
    log_file << std::endl << "[Boundary]" << std::endl;
    std::vector<bool> bc_absorbing; // { booleans };  // get the element as an array
    auto bc_absorb = tbl["Boundary"];
    status = get_toml_array(*bc_absorb.as_table(), "bc_absorbing", bc_absorbing);
    double treg = tbl["Boundary"]["treg"].value_or(double(300.0));
    if (stationary) { treg = 0.0; }
    log_file << "treg = " << treg << std::endl;

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    double dx = tbl["Numerics"]["dx"].value_or(double(60.0)); // Grid size [m]
    log_file << "dx = " << dx << std::endl;
    double dy = tbl["Numerics"]["dy"].value_or(double(00.0)); // Grid size [m]
    if (dy == 0.0) { double dy = dx; }
    log_file << "dy = " << dy << std::endl;
    double theta = tbl["Numerics"]["theta"].value_or(double(0.501));
    if (stationary) { theta = 1.0; }
    log_file << "theta = " << theta << std::endl;
    int iter_max = tbl["Numerics"]["iter_max"].value_or(int(50));
    log_file << "iter_max = " << iter_max << std::endl;
    double eps_newton = tbl["Numerics"]["eps_newton"].value_or(double(1.0e-12));
    log_file << "eps_newton = " << eps_newton << std::endl;
    double eps_bicgstab = tbl["Numerics"]["eps_bicgstab"].value_or(double(1.0e-12));
    log_file << "eps_bicgstab = " << eps_bicgstab << std::endl;
    bool regularization_init = tbl["Numerics"]["regularization_init"].value_or(bool(true));
    log_file << "regularization_init = " << regularization_init << std::endl;
    bool regularization_iter = tbl["Numerics"]["regularization_iter"].value_or(bool(true));
    log_file << "regularization_iter = " << regularization_iter << std::endl;
    bool regularization_time = tbl["Numerics"]["regularization_time"].value_or(bool(true));
    log_file << "regularization_time = " << regularization_time << std::endl;

    // Output
    log_file << std::endl << "[Output]" << std::endl;
    double dt_his = tbl["Output"]["dt_his"].value_or(double(1.0));
    if (stationary)
    {
        log_file << "dt_his = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_his = " << dt_his << std::endl;
    }
    double dt_map = tbl["Output"]["dt_map"].value_or(double(10.0));
    if (stationary)
    {
        log_file << "dt_map = " << 1.0 << std::endl;
    }
    else
    {
        log_file << "dt_map = " << dt_map << std::endl;
    }
    
    double dxinv = 1. / dx;                               // invers grid size [m]
    double dyinv = 1. / dy;                               // invers grid size [m]
    int nx = int(Lx/dx) + 1 + 2;                          // number of nodes in x-direction; including 2 virtual points
    int ny = int(Ly/dy) + 1 + 2;                          // number of nodes in y-direction; including 2 virtual points
    int nxny = nx * ny;                                   // total number of nodes
    double dxdy = dx * dy ;                               // area of control volume
    //if (viscosity == 0.0)
    //{
    //    viscosity = 0.2 * std::sqrt(dx*dx + dy*dy);
    //}

    // Time 
    int total_time_steps = 0;  // Number of time steps [-]
    double dtinv;                                         // Inverse of dt, if dt==0 then stationary solution [1/s]
    int wrihis;                                           // write interval to his-file
    int wrimap;                                           // write interval to map-file

    if (stationary)
    {
        dt = 0.0;
        dtinv = 0.0;                                      // stationary solution
        theta = 1.0;                                      // Stationary solution
        tstop = 1.;
        total_time_steps = 2;                             // initial step (step 1), stationary result (step 2)
        wrihis = 1;                                       // write interval to his-file
        wrimap = 1;                                       // write interval to map-file
    }
    else
    {
        dtinv = 1. / dt;                                  // Inverse of dt [1/s]
        total_time_steps = int((tstop - tstart) * dtinv) + 1;  // Number of time steps [-]
        theta = 0.501;                                    // Implicitness factor (0.5 <= theta <= 1.0)
        wrihis = std::max(int(dt * dtinv), int(dt_his * dtinv));      // write interval to his-file (every delta t)
        wrimap = std::max(int(dt * dtinv), int(dt_map * dtinv));     // write interval to map-file (every 1 sec , or every delta t)
    }

    double alpha = 1. / 8.;                               // Linear (spatial) interpolation coefficient
    double g = 10.0;                                      // Gravitational acceleration

    // Coefficients used in the boundary discretization

    double alpha_bc = 0.0;  // (2. * alpha - 1. / 2.);  // = 2 * 1/8 - 1/2 = -1/4
    alpha_bc = -0.25;

    std::vector<double> x(nxny, 0.);                      // x-coordinate
    std::vector<double> y(nxny, 0.);                      // y-coordinate
    std::vector<double> zb(nxny, 0.);                   // regularized bed level
    std::vector<double> zb_giv(nxny, -10.);               // initial bed level

    std::vector<double> hn(nxny, 0.);                     // water depth at (n)
    std::vector<double> qn(nxny, 0.);                     // q-flux at (n)
    std::vector<double> rn(nxny, 0.);                     // r-flux at (n)
    std::vector<double> hp(nxny, 0.);                     // total depth at (n+1,p), previous iteration
    std::vector<double> qp(nxny, 0.);                     // x-flow at (n+1,p), previous iteration
    std::vector<double> rp(nxny, 0.);                     // y-flow at (n+1,p), previous iteration
    std::vector<double> dh(nxny, 0.);                     // delta for water depth
    std::vector<double> dq(nxny, 0.);                     // delta for q-flux
    std::vector<double> dr(nxny, 0.);                     // delta for r-flux
    std::vector<double> s(nxny, 0.);                     // water level, needed for post-processing
    std::vector<double> u(nxny, 0.);                     // u-velocity, needed for post-processing
    std::vector<double> v(nxny, 0.);                     // v-velocity, needed for post-processing
    std::vector<double> s_giv(nxny, 0.);  // given initial water level
    std::vector<double> u_giv(nxny, 0.);  // given initial velocity
    std::vector<double> v_giv(nxny, 0.);  // given initial velocity
    std::vector<double> delta_h(nxny, 0.);
    std::vector<double> delta_q(nxny, 0.);
    std::vector<double> delta_r(nxny, 0.);
    std::vector<double> mass(3, 0.);  // weighting coefficients of the mass-matrix in x-direction

    Eigen::VectorXd solution(3 * nxny);                         // solution vector [h, q, r]^{n}
    Eigen::VectorXd rhs(3 * nxny);                          // RHS vector [h, q, r]^n

    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    std::vector<double> w_bc(3, 0.0);
    w_bc[0] = 1.0 + alpha_bc;
    w_bc[1] = 1.0 - 2.0 * alpha_bc;
    w_bc[2] = alpha_bc;

    //initialize water level
    std::cout << "Initialisation" << std::endl;
    status = 0;
    //status = read_bed_level(bed_filename, zb);
    double min_zb = *std::min_element(zb.begin(), zb.end());

    log_file << "=======================================================" << std::endl;
    log_file << "Nodes   : " << nx << "x" << ny << "=" << nxny << std::endl;
    log_file << "Elements: " << (nx - 1) * (ny - 1) << std::endl;
    log_file << "Volumes : " << (nx - 2) * (ny - 2) << std::endl;
    log_file << "CFL     : " << std::sqrt(g * std::abs(min_zb)) * dt * std::sqrt(( 1./(dx*dx) + 1./(dy*dy))) << std::endl;
    log_file << "=======================================================" << std::endl;
    std::cout << "    LxLy: " << Lx << "x" << Ly << std::endl;
    std::cout << "    dxdy: " << dx << "x" << dy << std::endl;
    std::cout << "    nxny: " << nx << "x" << ny << "=" << nxny << std::endl;

    //initialize x- and y-coordinate
    int k = 0;
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            k = j * nx + i;
            x[k] = double(i) * dx - Lx / 2 - dx;
            y[k] = double(j) * dy - Ly / 2 - dx;
        }
    }
    //initialize water level
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            k = j * nx + i;
            if (ini_vars[0] == "zeta" && ini_vars[1] == "zeta_GaussHump")
            {
                s_giv[k] = a0 * std::exp(-((x[k] - gauss_mu) * (x[k] - gauss_mu) + (y[k] - gauss_mu) * (y[k] - gauss_mu)) / (2. * gauss_sigma * gauss_sigma));  // initial water level
            }
            //s_giv[k] = 2. * a0 * std::exp(100.0 * log(0.01) * (x[k] * x[k] + y[k] * y[k]) / (Lx * Lx));
            //s_giv[k] = 0.0;  // needed for stationary test case
            //s_giv[k] = 2. * a0 * std::exp(-(x[k] * x[k]) / (500. * 500.));  // initial water level
            //s_giv[k] = 2. * a0 * std::exp(-(y[k] * y[k]) / (500. * 500.));  // initial water level
            //s_giv[k] = -a0 * x[i] * 0.01;  // initial water level fixed gradient
            //s_giv[k] = a0 * y[k] * 0.01;  // initial water level fixed gradient
            u_giv[k] = 0.0;
            v_giv[k] = 0.0;
            zb[k] = zb_giv[k];
            s[k] = s_giv[k];
            u[k] = u_giv[k];
            v[k] = v_giv[k];
        }
    }
    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        //(void)regular->given_function_2d(z_b, psi, eq8, zb_ini, dx, c_psi, use_eq8);
        STOP_TIMER(Regularization_init);
    }
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            k = j * nx + i;
            hn[k] = s[k] - zb[k];  // Initial water depth
            qn[k] = hn[k] * u[k];  // Initial q=hu -velocity
            rn[k] = hn[k] * v[k];  // Initial r=hv -velocity
            hp[k] = hn[k]; 
            qp[k] = qn[k]; 
            rp[k] = rn[k]; 
        }
    }
    if (regularization_init)
    {
        START_TIMER(Regularization_init);
        //(void)regular->given_function(visc_reg, psi, eq8, visc_given, dx, c_psi, use_eq8);
        STOP_TIMER(Regularization_init);
    }

    // Merge water depth (h) and (q,r) flux in one solution-vector
    for (int i = 0; i < 3 * nxny; ++i)
    {
        solution[i] = 0.0;  // Delta h, \Delta q and \Delta r
        rhs[i] = 0.0;  // 
    }

    ////////////////////////////////////////////////////////////////////////////
    // Define map file 
    std::cout << "Define map-file" << std::endl;
    std::string nc_mapfile(map_filename);
    UGRID2D* map_file = new UGRID2D();
    status = map_file->open(nc_mapfile);
    status = map_file->mesh2d();

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
    for (int j = 0; j < ny - 1; ++j)
    {
        for (int i = 0; i < nx - 1; ++i)
        {
            if (node_mask[j * nx + i] <= 1 || node_mask[j * nx + i + 1] <= 1)
            {
                mesh2d_edge_nodes.push_back(j * nx + i);
                mesh2d_edge_nodes.push_back(j * nx + i + 1);
                node_mask[j * nx + i] += 1;
                node_mask[j * nx + i + 1] += 1;
            }

            if (node_mask[j * nx + i + 1] <= 1 || node_mask[(j + 1) * nx + i + 1] <= 1)
            {
                mesh2d_edge_nodes.push_back(j * nx + i + 1);
                mesh2d_edge_nodes.push_back((j + 1) * nx + i + 1);
                node_mask[j * nx + i + 1] += 1;
                node_mask[(j + 1) * nx + i + 1] += 1;
            }

            if (node_mask[(j + 1) * nx + i + 1] <= 1 || node_mask[(j + 1) * nx + i] <= 1)
            {
                mesh2d_edge_nodes.push_back((j + 1) * nx + i + 1);
                mesh2d_edge_nodes.push_back((j + 1) * nx + i);
                node_mask[(j + 1) * nx + i + 1] += 1;
                node_mask[(j + 1) * nx + i] += 1;
            }

            if (node_mask[(j + 1) * nx + i] <= 1 || node_mask[j * nx + i] <= 1)
            {
                mesh2d_edge_nodes.push_back((j + 1) * nx + i);
                mesh2d_edge_nodes.push_back(j * nx + i);
                node_mask[(j + 1) * nx + i] += 1;
                //node_mask[j * nx + i] += 1;
            }
            mesh2d_face_nodes.push_back(j * nx + i);
            mesh2d_face_nodes.push_back(j * nx + i + 1);
            mesh2d_face_nodes.push_back((j + 1) * nx + i + 1);
            mesh2d_face_nodes.push_back((j + 1) * nx + i);
        }
    }
    status = map_file->put_variable_2("mesh2d_edge_nodes", mesh2d_edge_nodes);
    status = map_file->put_variable_4("mesh2d_face_nodes", mesh2d_face_nodes);

    dim_names.clear();
    dim_names.push_back("mesh2d_nNodes");
    status = map_file->add_variable("mesh2d_node_x", dim_names, "projection_x_coordinate", "x", "m", "mesh2D", "node");
    status = map_file->add_variable("mesh2d_node_y", dim_names, "projection_y_coordinate", "y", "m", "mesh2D", "node");
    status = map_file->put_variable("mesh2d_node_x", x);
    status = map_file->put_variable("mesh2d_node_y", y);

    // Compute mass centres of faces
    std::vector<double> xmc;
    std::vector<double> ymc;
    for (int j = 0; j < ny - 1; ++j)
    {
        for (int i = 0; i < nx - 1; ++i)
        {
            int p0 = p_index(i, j, nx);
            int p1 = p_index(i + 1, j, nx);
            int p2 = p_index(i + 1, j + 1, nx);
            int p3 = p_index(i, j + 1, nx);;
            xmc.push_back(0.25 * (x[p0] + x[p1] + x[p2] + x[p3]));
            ymc.push_back(0.25 * (y[p0] + y[p1] + y[p2] + y[p3]));
        }
    }
    dim_names.clear();
    dim_names.push_back("mesh2d_nFaces");
    status = map_file->add_variable("mesh2d_face_x", dim_names, "projection_x_coordinate", "x", "m", "mesh2D", "face");
    status = map_file->add_variable("mesh2d_face_y", dim_names, "projection_y_coordinate", "y", "m", "mesh2D", "face");
    status = map_file->put_variable("mesh2d_face_x", xmc);
    status = map_file->put_variable("mesh2d_face_y", ymc);

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
    status = map_file->add_variable(map_h_name, dim_names, "sea_floor_depth_below_sea_surface", "Water depth", "m", "mesh2D", "node");
    status = map_file->add_variable(map_q_name, dim_names, "", "Water flux", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_q_name, dim_names, "", "Water flux", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dh_name, dim_names, "", "Delta h^{n+1,p+1}", "m", "mesh2D", "node");
    status = map_file->add_variable(map_dq_name, dim_names, "", "Delta q^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_dr_name, dim_names, "", "Delta r^{n+1,p+1}", "m2 s-1", "mesh2D", "node");
    status = map_file->add_variable(map_s_name, dim_names, "sea_surface_height_above_geoid", "WaterLevel", "m", "mesh2D", "node");
    status = map_file->add_variable(map_u_name, dim_names, "sea_water_x_velocity", "Velocity (x)", "m s-1", "mesh2D", "node");
    status = map_file->add_variable(map_v_name, dim_names, "sea_water_y_velocity", "Velocity (y)", "m s-1", "mesh2D", "node");


    // Put data on map file
    int nst_map = 0;
    map_file->put_time(nst_map, double(0) * dt);
    map_file->put_time_variable(map_h_name, nst_map, hn);
    map_file->put_time_variable(map_q_name, nst_map, qn);
    map_file->put_time_variable(map_r_name, nst_map, rn);
    map_file->put_time_variable(map_dh_name, nst_map, dh);
    map_file->put_time_variable(map_dq_name, nst_map, dq);
    map_file->put_time_variable(map_dr_name, nst_map, dr);
    map_file->put_time_variable(map_s_name, nst_map, s);
    map_file->put_time_variable(map_u_name, nst_map, u);
    map_file->put_time_variable(map_v_name, nst_map, v);


    ////////////////////////////////////////////////////////////////////////////
    // Define time history file
    std::cout << "Define his-file" << std::endl;
    std::string nc_hisfile(his_filename);
    CFTS* his_file = new CFTS();
    status = his_file->open(nc_hisfile, model_title);

    // Initialize observation station locations
    int centre = p_index(nx / 2, ny / 2, nx);
    int w_bnd = p_index(1, ny / 2, nx);  // at boundary, not at virtual point
    int e_bnd = p_index(nx - 2, ny / 2, nx);  // at boundary, not at virtual point
    int s_bnd = p_index(nx / 2, 1, nx);
    int n_bnd = p_index(nx / 2, ny - 2, nx);

    int sw_bnd = p_index(1, 1, nx);
    int ne_bnd = p_index(nx -2, ny-2, nx);
    int nw_bnd = p_index(1, ny - 2, nx);
    int se_bnd = p_index(nx - 2, 1, nx);
    int p_a;
    int p_b;
    int p_c;
    int p_d;
    if (Lx < 2500.0)
    {
        p_a = p_index(int((5. / dx) + nx / 2), int((0. / dx) + ny / 2), nx);
        p_b = p_index(int((4. / dx) + nx / 2), int((3. / dx) + ny / 2), nx);
        p_c = p_index(int((3. / dx) + nx / 2), int((4. / dx) + ny / 2), nx);
        p_d = p_index(int((0. / dx) + nx / 2), int((5. / dx) + ny / 2), nx);
    }
    else
    {
        if (int(2500. / dx) * dx != 2500. || int(2500. / dy) * dy != 2500.)
        {
            std::cout << "----------------------------" << std::endl;
            std::cout << "dx=" << dx << " or dy=" << dy << " is not a divider of 2500 [m]" << std::endl;
            std::cout << "Press Enter to finish";
            std::cin.ignore();
            exit(1);
        }
        p_a = p_index(int((2500. / dx) + nx / 2), int((0. / dx) + ny / 2), nx);
        p_b = p_index(int((2000. / dx) + nx / 2), int((1500. / dx) + ny / 2), nx);
        p_c = p_index(int((1500. / dx) + nx / 2), int((2000. / dx) + ny / 2), nx);
        p_d = p_index(int((0. / dx) + nx / 2), int((2500. / dx) + ny / 2), nx);
    }
    std::vector<double> x_obs = { x[p_a], x[p_b], x[p_c], x[p_d], x[centre], x[n_bnd], x[ne_bnd], x[e_bnd], x[se_bnd], x[s_bnd], x[sw_bnd], x[w_bnd], x[nw_bnd] };
    std::vector<double> y_obs = { y[p_a], y[p_b], y[p_c], y[p_d], y[centre], y[n_bnd], y[ne_bnd], y[e_bnd], y[se_bnd], y[s_bnd], y[sw_bnd], y[w_bnd], y[nw_bnd] };

    std::vector<std::string> obs_stations;
    obs_stations.push_back("A");
    obs_stations.push_back("B");
    obs_stations.push_back("C");
    obs_stations.push_back("D");
    obs_stations.push_back("Centre");
    obs_stations.push_back("N station");
    obs_stations.push_back("NE station");
    obs_stations.push_back("E station");
    obs_stations.push_back("SE station");
    obs_stations.push_back("S station");
    obs_stations.push_back("SW station");
    obs_stations.push_back("W station");
    obs_stations.push_back("NW station");
    his_file->add_stations(obs_stations, x_obs, y_obs);
    his_file->add_time_series();

    std::string his_s_name("water_level");
    his_file->add_variable(his_s_name, "sea_surface_height", "Water level", "m");
    std::string his_u_name("u-velocity");
    his_file->add_variable(his_u_name, "sea_water_x_velocity", "Velocity xdir", "m s-1");
    std::string his_v_name("v-velocity");
    his_file->add_variable(his_v_name, "sea_water_y_velocity", "Velocity ydir", "m s-1");

    // Put data on time history file
    int nst_his = 0;
    his_file->put_time(nst_his, double(0) * dt);

    std::vector<double> his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd] };
    his_file->put_variable(his_s_name, nst_his, his_values);

    his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd] };
    his_file->put_variable(his_u_name, nst_his, his_values);

    his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd] };
    his_file->put_variable(his_v_name, nst_his, his_values);

    std::string his_newton_iter_name("his_newton_iterations");
    his_file->add_variable_without_location(his_newton_iter_name, "iterations", "Newton iteration", "-");
    his_values = { std::numeric_limits<double>::quiet_NaN() };
    his_file->put_variable(his_newton_iter_name, nst_his, his_values);

    std::string his_bicgstab_iter_name("his_bicgstab_iterations");
    his_file->add_variable_without_location(his_bicgstab_iter_name, "iterations", "BiCGStab iteration", "-");
    his_values = { std::numeric_limits<double>::quiet_NaN() };
    his_file->put_variable(his_bicgstab_iter_name, nst_his, his_values);

    std::string his_bicgstab_iter_error_name("bicgstab_iteration_error");
    his_file->add_variable_without_location(his_bicgstab_iter_error_name, "iteration_error", "BiCGStab iteration error", "-");
    his_values = { std::numeric_limits<double>::quiet_NaN() };
    his_file->put_variable(his_bicgstab_iter_error_name, nst_his, his_values);

    Eigen::SparseMatrix<double> A(3 * nxny, 3 * nxny);
    for (int i = 0; i < 3 * nxny; ++i) 
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = solution[i];
    }
    
    std::cout << "Start time-loop" << std::endl;
    std::cout << std::fixed << std::setprecision(2) << tstart + dt * double(0) << ";   " << tstart + tstop << std::endl;

    STOP_TIMER(Initialization);
    if (total_time_steps <= 1)
    {
        std::cout << "No time loop performed, due to total_time_steps <= 1" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
    }

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
        double hn_im2, hn_im1, hn_i, hn_ip1, hn_ip2, hn_jm1, hn_jp1;  // previous timestep (n)
        double hn_im12, hn_ip12, hn_jm12, hn_jp12;
        double hp_im2, hp_im1, hp_i, hp_ip1, hp_ip2, hp_jm1, hp_jp1;  // previous iteration (p)
        double hp_im12, hp_ip12, hp_jm12, hp_jp12;

        double qn_im2, qn_im1, qn_i, qn_ip1, qn_ip2, qn_jm1, qn_jp1;  // timestep (n)
        double qn_im12, qn_ip12;
        double qp_im2, qp_im1, qp_i, qp_ip1, qp_ip2, qp_jm1, qp_jp1;  // iteration (p)
        double qp_im12, qp_ip12;

        double rn_im1, rn_i, rn_ip1, rn_jm1, rn_jp1;  // timestep (n)
        double rn_jm12, rn_jp12;
        double rp_im1, rp_i, rp_ip1, rp_jm1, rp_jp1;  // iteration (p)
        double rp_jm12, rp_jp12;

        double htheta_i;
        double htheta_im1, htheta_im2, htheta_ip1, htheta_ip2;
        double htheta_im12, htheta_ip12;
        double htheta_jm1, htheta_jm2, htheta_jp1, htheta_jp2;
        double htheta_jm12, htheta_jp12;
        double qtheta_i;
        double qtheta_im1, qtheta_im2, qtheta_ip1, qtheta_ip2;
        double qtheta_im12, qtheta_ip12;
        double rtheta_i;
        double rtheta_jm1, rtheta_jp1;
        double rtheta_jm12, rtheta_jp12;

        // compute water level and velocities at t=n+1
        int used_iter = 0;
        int used_lin_iter = 0;

        START_TIMER(Newton iteration);
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        for (int iter = 0; iter < iter_max; ++iter)
        {
            if (nst == 1 && iter == 0)
            {
                START_TIMER(Matrix initialization);
            }
            //boundary nodes
            {
                //
                // north boundary
                //
                int j = ny - 1;
                for (int i = 1; i < nx - 1; ++i)
                {
                    int ph = p_index(i, j, nx);  // continuity equation
                    int ph_s = p_index(i, j - 1, nx);  // continuity equation
                    int ph_ss = p_index(i, j - 2, nx);  // continuity equation
                    if (true)
                    {
                        // Dirichlet
                        //
                        // continuity equation
                        //
                        A.coeffRef(3 * ph, 3 * ph_s) = -1.0;
                        A.coeffRef(3 * ph, 3 * ph) = 1.0;
                        rhs[ph] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_s + 1) = 0.5;
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 0.5;
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_s + 2) = 0.5;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = 0.5;
                        rhs[3 * ph + 2] = 0.0;
                    }
                    else
                    {
                        // weakly reflective
                        hn_i = hn[ph];
                        hn_jm1 = hn[ph_s];
                        hn_jm12 = 0.5 * (hn_i + hn_jm1);

                        hp_i = hp[ph];
                        hp_jm1 = hp[ph_s];
                        hp_jm12 = 0.5 * (hp_i + hp_jm1);

                        rn_i = rn[ph];
                        rn_jm1 = rn[ph_s];

                        rp_i = rp[ph];
                        rp_jm1 = rp[ph_s];
                        rp_jm12 = 0.5 * (rp_i + rp_jm1);

                        rtheta_i = theta * rp_i + (1.0 - theta) * rn_i;
                        rtheta_jm1 = theta * rp_jm1 + (1.0 - theta) * rn_jm1;
                        //
                        // continuity equation
                        //
                        double c_wave = std::sqrt(g * hp_jm12);
                        A.coeffRef(3 * ph, 3 * ph) = -0.75 * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph_s) = -0.75 * theta * c_wave;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph + 2) = 0.5 * theta;
                        A.coeffRef(3 * ph, 3 * ph_s + 2) = 0.5 * theta;
                        //
                        rhs[3 * ph] = -rp_jm12 + c_wave * hp_jm12 - c_wave * std::abs(zb[i]);
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_s + 1) = -1.0;
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 1.0;
                        rhs[ph + 1] = 0.0;
                        if (do_convection)
                        {
                            int a = 1;
                        }
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph) = +c_wave * 0.5 * dtinv;
                        A.coeffRef(3 * ph + 2, 3 * ph_s) = +c_wave * 0.5 * dtinv;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = -c_wave * dy * theta;
                        A.coeffRef(3 * ph + 2, 3 * ph_s + 2) = +c_wave * dy * theta;
                        rhs[3 * ph + 2] = +c_wave * (-dtinv * (hp_jm12 - hn_jm12)
                            - dy * rtheta_jm1 + dy * rtheta_i);
                        if (do_convection)
                        {
                            int a = 1;
                        }
                    }
                }
                //
                // east boundary
                //
                int i = nx - 1;
                for (int j = 1; j < ny - 1; ++j)
                {
                    int ph = p_index(i, j, nx);  // continuity equation
                    int ph_w = p_index(i - 1, j, nx);  // continuity equation
                    int ph_ww = p_index(i - 2, j, nx);  // continuity equation
                    if (true)
                    {
                        // Dirichlet
                        //
                        // continuity equation
                        //
                        A.coeffRef(3 * ph, 3 * ph_w) = -1.0;
                        A.coeffRef(3 * ph, 3 * ph) = 1.0;
                        rhs[3 * ph] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 1) = 0.5;
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 0.5;
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_w + 2) = -1.0;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = 1.0;
                        rhs[3 * ph + 2] = 0.0;
                    }
                    else
                    {
                        // Weakly reflective
                        hn_i = hn[ph];
                        hn_im1 = hn[ph_w];
                        hn_im2 = hn[ph_ww];
                        hn_im12 = 0.5 * (hn_i + hn_im1) + 0.5 * alpha_bc * (hn_im2 - 2. * hn_im1 + hn_i);;

                        hp_i = hp[ph];
                        hp_im1 = hp[ph_w];
                        hp_im2 = hp[ph_ww];
                        hp_im12 = 0.5 * (hp_i + hp_im1) + 0.5 * alpha_bc * (hp_im2 - 2. * hp_im1 + hp_i);

                        htheta_i = theta * hp_i + (1.0 - theta) * hn_i;
                        htheta_im1 = theta * hp_im1 + (1.0 - theta) * hn_im1;
                        htheta_im2 = theta * hp_im2 + (1.0 - theta) * hn_im2;
                        htheta_im12 = 0.5 * (htheta_i + htheta_im1) + 0.5 * alpha_bc * (htheta_im2 - 2.0 * htheta_im1 + htheta_i);

                        qn_i = qn[ph];
                        qn_im1 = qn[ph_w];
                        qn_im2 = qn[ph_ww];
                        qn_im12 = 0.5 * (qn_i + qn_im1) + 0.5 * alpha_bc * (qn_im2 - 2. * qn_im1 + qn_i);

                        qp_i = qp[ph];
                        qp_im1 = qp[ph_w];
                        qp_im2 = qp[ph_ww];
                        qp_im12 = 0.5 * (qp_i + qp_im1) + 0.5 * alpha_bc * (qp_im2 - 2. * qp_im1 + qp_i);

                        qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                        qtheta_im1 = theta * qp_im1 + (1.0 - theta) * qn_im1;
                        qtheta_im12 = 0.5 * (qtheta_im1 + qtheta_i);
                        //
                        A.coeffRef(3 * ph, 3 * ph) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_w) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ww) = 0.0;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph + 1) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_w + 1) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ww + 1) = 0.0;
                        //
                        A.coeffRef(3 * ph, 3 * ph + 2) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_w + 2) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ww + 2) = 0.0;
                        //
                        rhs[3 * ph] = 0.0;
                        //
                        //  Essential boundary condition
                        //
                        double h_infty = std::abs(0.5 * (zb[ph_w] + zb[ph]));
                        double c_wave = std::sqrt(g * h_infty);
                        A.coeffRef(3 * ph, 3 * ph_ww) = -0.75 * w_bc[2] * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph_w) = -0.75 * w_bc[1] * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph) = -0.75 * w_bc[0] * theta * c_wave;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph_ww + 1) = 0.5 * w_bc[2] * theta;
                        A.coeffRef(3 * ph, 3 * ph_w + 1) = 0.5 * w_bc[1] * theta;
                        A.coeffRef(3 * ph, 3 * ph + 1) = 0.5 * w_bc[0] * theta;
                        //
                        rhs[3 * ph] = -qp_im12 + c_wave * hp_im12 - c_wave * h_infty;
                        //
                        // Natural boundary condition (q-momentum part)
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_w) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_ww) = 0.0;
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 1) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 1) = 0.0;
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph + 2) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 2) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_ww + 2) = 0.0;
                        //
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // q-momentum + c_wave * continuity
                        //
                        double mom_q_fac = 1.0;
                        double mom_r_fac = 0.0;
                        double con_fac = c_wave;
                        //
                        // q-momentum part
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_ww + 1) = mom_q_fac * 0.5 * dtinv * alpha_bc;
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 1) = mom_q_fac * 0.5 * dtinv * w_bc[1];
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = mom_q_fac * 0.5 * dtinv * w_bc[0];
                        A.coeffRef(3 * ph + 1, 3 * ph_w) = mom_q_fac * dy * -dxinv * 0.5 * theta * g * htheta_im12;
                        A.coeffRef(3 * ph + 1, 3 * ph) = mom_q_fac * dy * dxinv * 0.5 * theta * g * htheta_im12;
                        rhs[3 * ph + 1] = mom_q_fac * (
                            - dtinv * 0.5 * (qp_i - qn_i) * w_bc[0]
                            - dtinv * 0.5 * (qp_im1 - qn_im1) * w_bc[1]
                            - dtinv * 0.5 * (qp_im2 - qn_im2) * alpha_bc
                            - dy * dxinv * g * htheta_im12 * (htheta_i - htheta_im1)
                            );
                        if (do_convection)
                        {
                            int a = 1;
                        }
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_w + 2) = -1.0;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = 1.0;
                        rhs[3 * ph + 2] = 0.0;
                        if (do_convection)
                        {
                            int a = 1;
                        }
                        //
                       // Natural boundary condition (continuity part added and multiplied by +c_wave)
                       //
                        A.coeffRef(3 * ph + 1, ph_ww) += con_fac * 0.5 * dtinv * w_bc[2];
                        A.coeffRef(3 * ph + 1, ph_w) += con_fac * 0.5 * dtinv * w_bc[1];
                        A.coeffRef(3 * ph + 1, ph) += con_fac * 0.5 * dtinv * w_bc[0];
                        // flow flux
                        A.coeffRef(3 * ph + 1, ph_w + 1) += con_fac * dy * -dxinv * theta;
                        A.coeffRef(3 * ph + 1, ph + 1) += con_fac * dy * dxinv * theta;
                        // rhs
                        rhs[3 * ph + 1] += con_fac * (
                            - dtinv * 0.5 * (hp_i - hn_i) * w_bc[0]
                            - dtinv * 0.5 * (hp_im1 - hn_im1) * w_bc[1]
                            - dtinv * 0.5 * (hp_im2 - hn_im2) * w_bc[2]
                            - dy * dxinv * (qtheta_i - qtheta_im1)
                            );
                   }
                }
                //
                // south boundary
                //
                j = 0;
                for (int i = 1; i < nx - 1; ++i)
                {
                    int ph = p_index(i, j, nx);  // continuity equation
                    int ph_n = p_index(i, j + 1, nx);  // continuity equation
                    int ph_nn =  p_index(i, j + 2, nx);  // continuity equation
                    if (true)
                    {
                        // Dirichlet
                        //
                        // continuity equation
                        //
                        A.coeffRef(3 * ph, 3 * ph_n) = 1.0;
                        A.coeffRef(3 * ph, 3 * ph) = -1.0;
                        rhs[3 * ph] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_n + 1) = 1.0;
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = -1.0;
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_n + 2) = 0.5;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = 0.5;
                        rhs[3 * ph + 2] = 0.0;
                    }
                    else
                    {
                        // Weakly reflective
                        hn_i = hn[ph];
                        hn_jp1 = hn[ph_n];
                        hn_jp12 = 0.5 * (hn_i + hn_jp1);

                        hp_i = hp[ph];
                        hp_jp1 = hp[ph_n];
                        hp_jp12 = 0.5 * (hp_i + hp_jp1);

                        rn_i = rn[ph];
                        rn_jp1 = rn[ph_n];

                        rp_i = rp[ph];
                        rp_jp1 = rp[ph_n];
                        rp_jp12 = 0.5 * (rp_i + rp_jp1);

                        rtheta_i = theta * rp_i + (1.0 - theta) * rn_i;
                        rtheta_jp1 = theta * rp_jp1 + (1.0 - theta) * rn_jp1;
                        // 
                        // continuity equation
                        //
                        double c_wave = std::sqrt(g * hp_jp12);
                        A.coeffRef(3 * ph, 3 * ph) = 0.75 * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph_n) = 0.75 * theta * c_wave;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph + 2) = 0.5 * theta;
                        A.coeffRef(3 * ph, 3 * ph_n + 2) = 0.5 * theta;
                        //
                        rhs[3 * ph] = -rp_jp12 - c_wave * hp_jp12 + c_wave * std::abs(zb[i]);
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph_n + 1) = 1.0;
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = -1.0;
                        rhs[3 * ph + 1] = 0.0;
                        if (do_convection)
                        {
                            int a = 1;
                        }
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph) = -c_wave * 0.5 * dtinv;
                        A.coeffRef(3 * ph + 2, 3 * ph_n) = -c_wave * 0.5 * dtinv;
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = c_wave * dy * theta;
                        A.coeffRef(3 * ph + 2, 3 * ph_n + 2) = -c_wave * dy * theta;
                        rhs[3 * ph + 2] = -c_wave * (-dtinv * (hp_jp12 - hn_jp12)
                            - dy * rtheta_jp1 + dy * rtheta_i);
                        if (do_convection)
                        {
                            int a = 1;
                        }
                    }
                }
                //
                // wwest boundary (2D)
                //
                i = 0;
                for (int j = 1; j < ny - 1; ++j)
                {
                    int ph = p_index(i, j, nx);  // continuity equation
                    int ph_e = p_index(i + 1, j, nx);  // q-momentum equation
                    int ph_ee = p_index(i + 2, j, nx);  // r-momentum equation
                    if (!bc_absorbing[BC_WEST])
                    {
                        // Dirichlet
                        //
                        // continuity equation
                        //
                        A.coeffRef(3 * ph, 3 * ph) = -1.0;
                        A.coeffRef(3 * ph, 3 * ph_e) = 1.0;
                        rhs[3 * ph] = 0.0;
                        //
                        // q-momentum
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 0.5;
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1) = 0.5;
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // r-momentum
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = -1.0;
                        A.coeffRef(3 * ph + 2, 3 * ph_e + 2) = 1.0;
                        rhs[3 * ph + 2] = 0.0;
                    }
                    if (bc_absorbing[BC_WEST])
                    {
                        // weakly reflective 
                        hn_i = hn[ph];
                        hn_ip1 = hn[ph_e];
                        hn_ip2 = hn[ph_ee];
                        hn_ip12 = 0.5 * (hn_i + hn_ip1) + 0.5 * alpha_bc * (hn_i - 2. * hn_ip1 + hn_ip2);
                        hp_i = hp[ph];
                        hp_ip1 = hp[ph_e];
                        hp_ip2 = hp[ph_ee];
                        hp_ip12 = 0.5 * (hp_i + hp_ip1) + 0.5 * alpha_bc * (hp_i - 2. * hp_ip1 + hp_ip2);
                        htheta_i = theta * hp_i + (1.0 - theta) * hn_i;
                        htheta_ip1 = theta * hp_ip1 + (1.0 - theta) * hn_ip1;
                        htheta_ip2 = theta * hp_ip2 + (1.0 - theta) * hn_ip2;
                        htheta_ip12 = 0.5 * (htheta_i + htheta_ip1) + 0.5 * alpha_bc * (htheta_i - 2. * htheta_ip1 + htheta_ip2);

                        qn_i = qn[ph];
                        qn_ip1 = qn[ph_e];
                        qn_ip2 = qn[ph_ee];
                        qn_ip12 = 0.5 * (qn_i + qn_ip1) + 0.5 * alpha_bc * (qn_i - 2. * qn_ip1 + qn_ip2);
                        qp_i = qp[ph];
                        qp_ip1 = qp[ph_e];
                        qp_ip2 = qp[ph_ee];
                        qp_ip12 = 0.5 * (qp_i + qp_ip1) + 0.5 * alpha_bc * (qp_i - 2. * qp_ip1 + qp_ip2);
                        qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                        qtheta_ip1 = theta * qp_ip1 + (1.0 - theta) * qn_ip1;
                        qtheta_ip2 = theta * qp_ip2 + (1.0 - theta) * qn_ip2;
                        qtheta_ip12 = 0.5 * (qtheta_i + qtheta_ip1) + 0.5 * alpha_bc * (qtheta_i - 2. * qtheta_ip1 + qtheta_ip2);
                        //
                        A.coeffRef(3 * ph, 3 * ph) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_e) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ee) = 0.0;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph + 1) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_e + 1) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ee + 1) = 0.0;
                        //
                        A.coeffRef(3 * ph, 3 * ph + 2) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_e + 2) = 0.0;
                        A.coeffRef(3 * ph, 3 * ph_ee + 2) = 0.0;
                        //
                        rhs[3 * ph] = 0.0;
                        //
                        // Essential boundary condition
                        //
                        double zb_ip12 = 0.5 * (zb[i] + zb[i + 1]) + 0.5 * alpha_bc * (zb[i] - 2. * zb[i + 1] + zb[i + 2]);
                        double h_infty = - zb_ip12;
                        double c_wave = std::sqrt(g * h_infty);
                        A.coeffRef(3 * ph, 3 * ph)    = 0.5 * w_bc[0] * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph_e)  = 0.5 * w_bc[1] * theta * c_wave;
                        A.coeffRef(3 * ph, 3 * ph_ee) = 0.5 * w_bc[2] * theta * c_wave;
                        // flow flux
                        A.coeffRef(3 * ph, 3 * ph + 1)    = 0.5 * theta * w_bc[0];
                        A.coeffRef(3 * ph, 3 * ph_e + 1)  = 0.5 * theta * w_bc[1];
                        A.coeffRef(3 * ph, 3 * ph_ee + 1) = 0.5 * theta * w_bc[2];
                        //
                        rhs[3 * ph] = -qtheta_ip12 - c_wave * (htheta_ip12 - h_infty);
                        //
                        // Natural boundary condition
                        //
                        // on the q-momentum row
                        A.coeffRef(3 * ph + 1, 3 * ph) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_e) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_ee) = 0.0;
                        // flow flux
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_ee + 1) = 0.0;
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph + 2) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 2) = 0.0;
                        A.coeffRef(3 * ph + 1, 3 * ph_ee + 2) = 0.0;
                        //
                        rhs[3 * ph + 1] = 0.0;
                        //
                        // q-momentum - c_wave * continuity
                        // 
                        double mom_fac = 1.0;
                        double con_fac = -c_wave;
                        //
                        // continuity part, Natural boundary condition
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph)    += con_fac * 0.5 * dtinv * w_bc[0];
                        A.coeffRef(3 * ph + 1, 3 * ph_e)  += con_fac * 0.5 * dtinv * w_bc[1];
                        A.coeffRef(3 * ph + 1, 3 * ph_ee) += con_fac * 0.5 * dtinv * w_bc[2];
                        // flow flux
                        A.coeffRef(3 * ph + 1, 3 * ph + 1)   += con_fac * dy * -dxinv * theta;
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1) += con_fac * dy *  dxinv * theta;
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1) += 0.0;
                        //
                        rhs[3 * ph + 1] += con_fac * (
                            - dtinv * 0.5 * (hp_i   - hn_i  ) * w_bc[0]
                            - dtinv * 0.5 * (hp_ip1 - hn_ip1) * w_bc[1]
                            - dtinv * 0.5 * (hp_ip2 - hn_ip2) * w_bc[2]
                            - dy * dxinv * (qtheta_ip1 -  qtheta_i)
                            );
                        //
                        // q-momentum part, Natural boundary condition
                        //
                        A.coeffRef(3 * ph + 1, 3 * ph)    += 0.5 * w_bc[0] * dxinv * g * theta * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]) - dy * dxinv * theta * g * htheta_ip12;
                        A.coeffRef(3 * ph + 1, 3 * ph_e)  += 0.5 * w_bc[1] * dxinv * g * theta * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]) + dy * dxinv * theta * g * htheta_ip12;
                        A.coeffRef(3 * ph + 1, 3 * ph_ee) += 0.5 * w_bc[2] * dxinv * g * theta * (htheta_ip1 + zb[i + 1] - htheta_i - zb[i]);

                        A.coeffRef(3 * ph + 1, 3 * ph + 1)    += 0.5 * dtinv * w_bc[0];
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1)  += 0.5 * dtinv * w_bc[1];
                        A.coeffRef(3 * ph + 1, 3 * ph_ee + 1) += 0.5 * dtinv * w_bc[2];

                        rhs[3 * ph + 1] += - (
                            + dtinv * 0.5 * (qp_i   - qn_i  ) * w_bc[0]
                            + dtinv * 0.5 * (qp_ip1 - qn_ip1) * w_bc[1]
                            + dtinv * 0.5 * (qp_ip2 - qn_ip2) * w_bc[2]
                            + dy * g * dxinv * htheta_ip12 * (htheta_ip1 + zb[ph] - htheta_i - zb[ph_e])
                            );
                        if (do_convection)
                        {
                            int a = 1;
                        }
                        //
                        // r-momentum, tangential to boundary (dr/dt ... = 0) Natural boundary condition)
                        //
                        // on the r-momentum row
                        A.coeffRef(3 * ph + 2, 3 * ph) = 0.0;
                        A.coeffRef(3 * ph + 2, 3 * ph_e) = 0.0;
                        A.coeffRef(3 * ph + 2, 3 * ph_ee) = 0.0;
                        // flow flux
                        A.coeffRef(3 * ph + 2, 3 * ph + 1) = 0.0;
                        A.coeffRef(3 * ph + 2, 3 * ph_e + 1) = 0.0;
                        A.coeffRef(3 * ph + 2, 3 * ph_ee + 1) = 0.0;
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) +=  0.5;
                        A.coeffRef(3 * ph + 2, 3 * ph_e + 2) +=  0.5;
                        A.coeffRef(3 * ph + 2, 3 * ph_ee + 2) = 0.0;

                        rhs[3 * ph + 2] = 0.0;
                        if (do_convection)
                        {
                            int a = 1;
                        }
                    }
                }
            }
            //corner nodes
            {
                // NE-corner
                int i = nx - 1;
                int j = ny - 1;
                int ph = 3 * p_index(i, j, nx);  // continuity equation
                int pq = 3 * p_index(i, j, nx) + 1;  // q-momentum equation
                int pr = 3 * p_index(i, j, nx) + 2;  // r-momentum equation
                A.coeffRef(ph, ph - 3) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, ph - 3 * nx) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, pq - 3) = 0.5;
                A.coeffRef(pq, pq) = -1.0;
                A.coeffRef(pq, pq - 3 * nx) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, pr - 3) = 0.5;
                A.coeffRef(pr, pr) = -1.0;
                A.coeffRef(pr, pr - 3 * nx) = 0.5;
                rhs[pr] = 0.0;
                // SE-corner
                i = nx - 1;
                j = 0;
                ph = 3 * p_index(i, j, nx);  // continuity equation
                pq = 3 * p_index(i, j, nx) + 1;  // q-momentum equation
                pr = 3 * p_index(i, j, nx) + 2;  // r-momentum equation
                A.coeffRef(ph, ph - 3) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, ph + 3 * nx) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, pq - 3) = 0.5;
                A.coeffRef(pq, pq) = -1.0;
                A.coeffRef(pq, pq + 3 * nx) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, pr - 3) = 0.5;
                A.coeffRef(pr, pr) = -1.0;
                A.coeffRef(pr, pr + 3 * nx) = 0.5;
                rhs[pr] = 0.0;
                // SW-corner
                i = 0;
                j = 0;
                ph = 3 * p_index(i, j, nx);  // continuity equation
                pq = 3 * p_index(i, j, nx) + 1;  // q-momentum equation
                pr = 3 * p_index(i, j, nx) + 2;  // r-momentum equation
                A.coeffRef(ph, ph + 3) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, ph + 3 * nx) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, pq + 3) = 0.5;
                A.coeffRef(pq, pq) = -1.0;
                A.coeffRef(pq, pq + 3 * nx) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, pr + 3) = 0.5;
                A.coeffRef(pr, pr) = -1.0;
                A.coeffRef(pr, pr + 3 * nx) = 0.5;
                rhs[pr] = 0.0;
                // NW-corner
                i = 0;
                j = ny - 1;
                ph = 3 * p_index(i, j, nx);  // continuity equation
                pq = 3 * p_index(i, j, nx) + 1;  // q-momentum equation
                pr = 3 * p_index(i, j, nx) + 2;  // r-momentum equation
                A.coeffRef(ph, ph - 3 * nx) = 0.5;
                A.coeffRef(ph, ph) = -1.0;
                A.coeffRef(ph, ph + 3) = 0.5;
                rhs[ph] = 0.0;
                A.coeffRef(pq, pq - 3 * nx) = 0.5;
                A.coeffRef(pq, pq) = -1.0;
                A.coeffRef(pq, pq + 3) = 0.5;
                rhs[pq] = 0.0;
                A.coeffRef(pr, pr - 3 * nx) = 0.5;
                A.coeffRef(pr, pr) = -1.0;
                A.coeffRef(pr, pr + 3) = 0.5;
                rhs[pr] = 0.0;
            }
            //
            // interior nodes
            //
            for (int i = 1; i < nx - 1; i++)
            {
                for (int j = 1; j < ny - 1; j++)
                {
                    if (i == nx / 2 + 1 && j == ny / 2 + 1)
                    {
                        int c = 1;
                    }
                    //log_file << "Node: " << i << "  " << j << std::endl;
                    int ph = p_index(i, j, nx);  // continuity equation
                    int ph_n = p_index(i, j + 1, nx);  // continuity equation
                    int ph_ne = p_index(i + 1, j + 1, nx);  // continuity equation
                    int ph_e = p_index(i + 1, j, nx);  // continuity equation
                    int ph_se = p_index(i + 1, j - 1, nx);  // continuity equation
                    int ph_s = p_index(i, j - 1, nx);  // continuity equation
                    int ph_sw = p_index(i - 1, j - 1, nx);  // continuity equation
                    int ph_w = p_index(i - 1, j, nx);  // continuity equation
                    int ph_nw = p_index(i - 1, j + 1, nx);  // continuity equation
                    int pq = ph + 1;  // q-momentum equation
                    int pr = ph + 2;  // r-momentum equation

                    hn_im1 = hn[ph_w];
                    hn_i = hn[ph];
                    hn_ip1 = hn[ph_e];
                    hn_jm1 = hn[ph_s];
                    hn_jp1 = hn[ph_n];

                    qn_im1 = qn[ph_w];
                    qn_i = qn[ph];
                    qn_ip1 = qn[ph_e];
                    qn_jm1 = qn[ph_s];
                    qn_jp1 = qn[ph_n];

                    rn_im1 = rn[ph_w];
                    rn_i = rn[ph];
                    rn_ip1 = rn[ph_e];
                    rn_jm1 = rn[ph_s];
                    rn_jp1 = rn[ph_n];

                    hp_im1 = hp[ph_w];
                    hp_i = hp[ph];
                    hp_ip1 = hp[ph_e];
                    hp_jm1 = hp[ph_s];
                    hp_jp1 = hp[ph_n];

                    qp_im1 = qp[ph_w];
                    qp_i = qp[ph];
                    qp_ip1 = qp[ph_e];
                    qp_jm1 = qp[ph_s];
                    qp_jp1 = qp[ph_n];

                    rp_im1 = rp[ph_w];
                    rp_i = rp[ph];
                    rp_ip1 = rp[ph_e];
                    rp_jm1 = rp[ph_s];
                    rp_jp1 = rp[ph_n];

                    htheta_im1 = theta * hp_im1 + (1.0 - theta) * hn_im1;
                    htheta_i   = theta * hp_i + (1.0 - theta) * hn_i;
                    htheta_ip1 = theta * hp_ip1 + (1.0 - theta) * hn_ip1;
                    htheta_ip12 = 0.5 * (htheta_ip1 + htheta_i);
                    htheta_im12 = 0.5 * (htheta_i + htheta_im1);
                    htheta_jp1 = theta * hp_jp1 + (1.0 - theta) * hn_jp1;
                    htheta_jm1 = theta * hp_jm1 + (1.0 - theta) * hn_jm1;
                    htheta_jp12 = 0.5 * (htheta_jp1 + htheta_i);
                    htheta_jm12 = 0.5 * (htheta_i + htheta_jm1);

                    qtheta_im1 = theta * qp_im1 + (1.0 - theta) * qn_im1;
                    qtheta_i = theta * qp_i + (1.0 - theta) * qn_i;
                    qtheta_ip1 = theta * qp_ip1 + (1.0 - theta) * qn_ip1;
                    qtheta_ip12 = 0.5 * (qtheta_ip1 + qtheta_i);
                    qtheta_im12 = 0.5 * (qtheta_i + qtheta_im1);

                    rtheta_jm1 = theta * rp_jm1 + (1.0 - theta) * rn_jm1;
                    rtheta_i = theta * rp_i + (1.0 - theta) * rn_i;
                    rtheta_jp1 = theta * rp_jp1 + (1.0 - theta) * rn_jp1;
                    rtheta_jp12 = 0.5 * (rtheta_jp1 + rtheta_i);
                    rtheta_jm12 = 0.5 * (rtheta_i + rtheta_jm1);

                    double htheta_0 = theta * hp[ph] + (1.0 - theta) * hn[ph];
                    double htheta_n = theta * hp[ph_n] + (1.0 - theta) * hn[ph_n];
                    double htheta_ne = theta * hp[ph_ne] + (1.0 - theta) * hn[ph_ne];
                    double htheta_e = theta * hp[ph_e] + (1.0 - theta) * hn[ph_e];
                    double htheta_se = theta * hp[ph_se] + (1.0 - theta) * hn[ph_se];
                    double htheta_s = theta * hp[ph_s] + (1.0 - theta) * hn[ph_s];
                    double htheta_sw = theta * hp[ph_sw] + (1.0 - theta) * hn[ph_sw];
                    double htheta_w = theta * hp[ph_w] + (1.0 - theta) * hn[ph_w];
                    double htheta_nw = theta * hp[ph_nw] + (1.0 - theta) * hn[ph_nw];

                    double qtheta_0 = theta * qp[ph] + (1.0 - theta) * qn[ph];
                    double qtheta_n = theta * qp[ph_n] + (1.0 - theta) * qn[ph_n];
                    double qtheta_ne = theta * qp[ph_ne] + (1.0 - theta) * qn[ph_ne];
                    double qtheta_e = theta * qp[ph_e] + (1.0 - theta) * qn[ph_e];
                    double qtheta_se = theta * qp[ph_se] + (1.0 - theta) * qn[ph_se];
                    double qtheta_s = theta * qp[ph_s] + (1.0 - theta) * qn[ph_s];
                    double qtheta_sw = theta * qp[ph_sw] + (1.0 - theta) * qn[ph_sw];
                    double qtheta_w = theta * qp[ph_w] + (1.0 - theta) * qn[ph_w];
                    double qtheta_nw = theta * qp[ph_nw] + (1.0 - theta) * qn[ph_nw];

                    double rtheta_0 = theta * rp[ph] + (1.0 - theta) * rn[ph];
                    double rtheta_n = theta * rp[ph_n] + (1.0 - theta) * rn[ph_n];
                    double rtheta_ne = theta * rp[ph_ne] + (1.0 - theta) * rn[ph_ne];
                    double rtheta_e = theta * rp[ph_e] + (1.0 - theta) * rn[ph_e];
                    double rtheta_se = theta * rp[ph_se] + (1.0 - theta) * rn[ph_se];
                    double rtheta_s = theta * rp[ph_s] + (1.0 - theta) * rn[ph_s];
                    double rtheta_sw = theta * rp[ph_sw] + (1.0 - theta) * rn[ph_sw];
                    double rtheta_w = theta * rp[ph_w] + (1.0 - theta) * rn[ph_w];
                    double rtheta_nw = theta * rp[ph_nw] + (1.0 - theta) * rn[ph_nw];
                    //
                    // continuity equation (dh/dt ... = 0)
                    //
                    if (do_continuity)
                    {
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
                        A.coeffRef(3 * ph, 3 * ph_sw) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph, 3 * ph_s) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph, 3 * ph_se) = dtinv * dxdy * alpha * alpha;
                        //
                        A.coeffRef(3 * ph, 3 * ph_w) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph, 3 * ph) = dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph, 3 * ph_e) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        //
                        A.coeffRef(3 * ph, 3 * ph_nw) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph, 3 * ph_n) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph, 3 * ph_ne) = dtinv * dxdy * alpha * alpha;
                        //
                        if (true)
                        {
                            A.coeffRef(3 * ph, 3 * ph_e + 1) = 0.5 * theta * dy;
                            A.coeffRef(3 * ph, 3 * ph + 1) = 0.5 * theta * dy - 0.5 * theta * dy;
                            A.coeffRef(3 * ph, 3 * ph_w + 1) = -0.5 * theta * dy;
                            //
                            A.coeffRef(3 * ph, 3 * ph_n + 2) = 0.5 * theta * dx;
                            A.coeffRef(3 * ph, 3 * ph + 2) = 0.5 * theta * dx - 0.5 * theta * dx;
                            A.coeffRef(3 * ph, 3 * ph_s + 2) = -0.5 * theta * dx;
                        }
                        else
                        {
                            A.coeffRef(3 * ph, 3 * ph + 1) = 0.0;  //
                            //theta*
                            //( 0.5 * dx * 0.125 * 3. - 0.5 * dx * 0.125 * 3.  // south
                            //+ 0.5 * dy * 0.125 * 3. + 0.5 * dy * 0.125 * 3.  // east
                            //+ 0.5 * dx * 0.125 * 3. + 0.5 * dx * 0.125 * 3.  // north
                            //- 0.5 * dy * 0.125 * 3. - 0.5 * dy * 0.125 * 3.)   // west  Is equal to zero, centre point does not count
                            //;
                            {  // brackets only useful for editor
                                A.coeffRef(3 * ph, 3 * ph_sw + 1) = theta * 0.125 * 0.5 * (-dy);  // flux_7
                                A.coeffRef(3 * ph, 3 * ph_s + 1) = theta * 0.125 * 0.5 * (dy - dy); // flux_2 and flux_7
                                A.coeffRef(3 * ph, 3 * ph_se + 1) = theta * 0.125 * 0.5 * (dy);  // flux_2
                                A.coeffRef(3 * ph, 3 * ph_e + 1) = theta * 0.125 * 0.5 * (dy) * (3. + 3.);  // flux_2 and flux_3
                                A.coeffRef(3 * ph, 3 * ph_ne + 1) = theta * 0.125 * 0.5 * (dy);  // flux_3
                                A.coeffRef(3 * ph, 3 * ph_n + 1) = theta * 0.125 * 0.5 * (dy - dy);  // flux_3 and flux_6
                                A.coeffRef(3 * ph, 3 * ph_nw + 1) = theta * 0.125 * 0.5 * (-dy);  // flux_6
                                A.coeffRef(3 * ph, 3 * ph_w + 1) = theta * 0.125 * 0.5 * (-dy) * (3. + 3.);  // flux_6 and flux_7
                            }
                            //
                            {  // brackets only useful for editor
                                A.coeffRef(3 * ph, 3 * ph_sw + 2) = theta * 0.125 * 0.5 * (-dx);  // flux_0
                                A.coeffRef(3 * ph, 3 * ph_s + 2) = theta * 0.125 * 0.5 * (-dx) * (3. + 3.);  // flux_0 and flux_1
                                A.coeffRef(3 * ph, 3 * ph_se + 2) = theta * 0.125 * 0.5 * (-dx);  // flux 1
                                A.coeffRef(3 * ph, 3 * ph_e + 2) = theta * 0.125 * 0.5 * (-dx + dx);  // flux_1 and flux_4
                                A.coeffRef(3 * ph, 3 * ph_ne + 2) = theta * 0.125 * 0.5 * (dx);  // flux_4
                                A.coeffRef(3 * ph, 3 * ph_n + 2) = theta * 0.125 * 0.5 * (dx) * (3. + 3.);  // flux_4 and flux_5
                                A.coeffRef(3 * ph, 3 * ph_nw + 2) = theta * 0.125 * 0.5 * (dx);  // flux_5
                                A.coeffRef(3 * ph, 3 * ph_w + 2) = theta * 0.125 * 0.5 * (-dx + dx);  // flux_0 and flux_5
                            }
                        }
                        //
                        // RHS continuity equation
                        //
                        // eight fluxes
                        //
                        //double flux_0 = 0.5 * dx * dcdx_scvf1(rtheta_0,  rtheta_s, rtheta_w, rtheta_sw);

                        double flux_0y = 0.25 * dx * dcdy_scvf2(rtheta_0, rtheta_s, rtheta_w, rtheta_sw);
                        double flux_1y = 0.25 * dx * dcdy_scvf2(rtheta_0, rtheta_s, rtheta_e, rtheta_se);
                        double flux_2x = 0.25 * dy * dcdx_scvf1(qtheta_e, qtheta_0, qtheta_se, qtheta_s);
                        double flux_3x = 0.25 * dy * dcdx_scvf1(qtheta_e, qtheta_0, qtheta_ne, qtheta_n);
                        double flux_4y = 0.25 * dx * dcdy_scvf2(rtheta_n, rtheta_0, rtheta_ne, rtheta_e);
                        double flux_5y = 0.25 * dx * dcdy_scvf2(rtheta_n, rtheta_0, rtheta_nw, rtheta_w);
                        double flux_6x = 0.25 * dy * dcdx_scvf1(qtheta_0, qtheta_w, qtheta_n, qtheta_nw);
                        double flux_7x = 0.25 * dy * dcdx_scvf1(qtheta_0, qtheta_w, qtheta_s, qtheta_sw);
                        //
                        rhs[3 * ph] = -dtinv * dxdy * alpha * alpha * (hp[ph_sw] - hn[ph_sw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (hp[ph_s] - hn[ph_s]) +
                            -dtinv * dxdy * alpha * alpha * (hp[ph_se] - hn[ph_se]) +
                            //
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (hp[ph_w] - hn[ph_w]) +
                            -dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha) * (hp[ph] - hn[ph]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (hp[ph_e] - hn[ph_e]) +
                            //
                            -dtinv * dxdy * alpha * alpha * (hp[ph_nw] - hn[ph_nw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (hp[ph_n] - hn[ph_n]) +
                            -dtinv * dxdy * alpha * alpha * (hp[ph_ne] - hn[ph_ne]);
                        //
                        double a = -(flux_0y + flux_1y + flux_2x + flux_3x + flux_4y + flux_5y + flux_6x + flux_7x);
                        double b = -dy * (qtheta_ip12 - qtheta_im12) - dx * (rtheta_jp12 - rtheta_jm12);
                        rhs[3 * ph] += b;
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
//                        log_file << std::scientific << std::setprecision(8) << "(" << i << "," << j << ") "
//                            << "dx:  " << (qtheta_ip12 - qtheta_im12) << "  dy:  " << (rtheta_jp12 - rtheta_jm12) << std::endl;
                    }
                    //
                    // q-momentum(dq/dt ... = 0) 
                    //
                    if (do_q_equation)
                    {
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
                        A.coeffRef(3 * ph + 1, 3 * ph_sw + 1) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph + 1, 3 * ph_s + 1) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 1, 3 * ph_se + 1) = dtinv * dxdy * alpha * alpha;
                        //              
                        A.coeffRef(3 * ph + 1, 3 * ph_w + 1) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 1, 3 * ph + 1) = dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 1, 3 * ph_e + 1) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        //            
                        A.coeffRef(3 * ph + 1, 3 * ph_nw + 1) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph + 1, 3 * ph_n + 1) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 1, 3 * ph_ne + 1) = dtinv * dxdy * alpha * alpha;
                        //
                        if (true)
                        {
                            A.coeffRef(3 * ph + 1, 3 * ph_e) = 0.5 * theta * g * dy * htheta_ip12;
                            A.coeffRef(3 * ph + 1, 3 * ph) = 0.5 * theta * g * dy * htheta_im12 - 0.5 * theta * g * dy * htheta_ip12;
                            A.coeffRef(3 * ph + 1, 3 * ph_w) = -0.5 * theta * g * dy * htheta_im12;
                        }
                        else
                        {  // brackets only useful for editor
                            double fac = theta * 0.125 * 0.0625 * 0.5 * dy;
                            A.coeffRef(3 * ph + 1, 3 * ph) = 0.0;  // fac* (-3.0 + 9.0 + 9.0 - 3.0 - 9.0 + 3.0 + 3.0 - 9.0);  // 
                            A.coeffRef(3 * ph + 1, 3 * ph_sw) = fac * (-3.0 + 1.0);  // 
                            A.coeffRef(3 * ph + 1, 3 * ph_s) = 0.0;  // fac* (-1.0 + 3.0 - 3.0 + 1.0); // 
                            A.coeffRef(3 * ph + 1, 3 * ph_se) = fac * (-1.0 + 3.0);  //
                            A.coeffRef(3 * ph + 1, 3 * ph_e) = fac * (-3.0 + 9.0 + 9.0 - 3.0);  // 
                            A.coeffRef(3 * ph + 1, 3 * ph_ne) = fac * (3.0 - 1.0);  // 
                            A.coeffRef(3 * ph + 1, 3 * ph_n) = 0.0;  // fac* (3.0 - 1.0 + 1.0 - 3.0);  //   
                            A.coeffRef(3 * ph + 1, 3 * ph_nw) = fac * (1.0 - 3.0);  // 
                            A.coeffRef(3 * ph + 1, 3 * ph_w) = fac * (-9.0 + 3.0 + 3.0 - 9.0);  //
                        }
                        //
                        // RHS q-momentum equation
                        //
                        double depth_0 = c_scv(htheta_0, htheta_w, htheta_sw, htheta_s);
                        double depth_1 = c_scv(htheta_0, htheta_s, htheta_se, htheta_e);
                        double depth_2 = c_scv(htheta_0, htheta_e, htheta_ne, htheta_n);
                        double depth_3 = c_scv(htheta_0, htheta_n, htheta_nw, htheta_w);

                        double grad_2x = 0.25 * dxdy / dx * dcdx_scvf1(htheta_e, htheta_0, htheta_se, htheta_s);
                        double grad_3x = 0.25 * dxdy / dx * dcdx_scvf1(htheta_e, htheta_0, htheta_ne, htheta_n);
                        double grad_6x = 0.25 * dxdy / dx * dcdx_scvf1(htheta_0, htheta_w, htheta_n, htheta_nw);
                        double grad_7x = 0.25 * dxdy / dx * dcdx_scvf1(htheta_0, htheta_w, htheta_s, htheta_sw);

                        rhs[3 * ph + 1] = -dtinv * dxdy * alpha * alpha * (qp[ph_sw] - qn[ph_sw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (qp[ph_s] - qn[ph_s]) +
                            -dtinv * dxdy * alpha * alpha * (qp[ph_se] - qn[ph_se]) +
                            //
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (qp[ph_w] - qn[ph_w]) +
                            -dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha) * (qp[ph] - qn[ph]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (qp[ph_e] - qn[ph_e]) +
                            //
                            -dtinv * dxdy * alpha * alpha * (qp[ph_nw] - qn[ph_nw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (qp[ph_n] - qn[ph_n]) +
                            -dtinv * dxdy * alpha * alpha * (qp[ph_ne] - qn[ph_ne]);
                        //
                        double a = -g * 0.25 * (depth_0 + depth_1 + depth_2 + depth_3) * (grad_2x + grad_3x + grad_6x + grad_7x);
                        double b = -0.5 * g * dy * (htheta_ip12 * htheta_ip12 - htheta_im12 * htheta_im12);
                        rhs[3 * ph + 1] += b;
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
                    }
                    // 
                    // r-momentum (dr/dt ... = 0)
                    //
                    if (do_r_equation)
                    {
                        A.coeffRef(3 * ph + 2, 3 * ph_sw + 2) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph + 2, 3 * ph_s + 2) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 2, 3 * ph_se + 2) = dtinv * dxdy * alpha * alpha;
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_w + 2) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 2, 3 * ph + 2) = dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 2, 3 * ph_e + 2) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        //
                        A.coeffRef(3 * ph + 2, 3 * ph_nw + 2) = dtinv * dxdy * alpha * alpha;
                        A.coeffRef(3 * ph + 2, 3 * ph_n + 2) = dtinv * dxdy * alpha * (1.0 - 2. * alpha);
                        A.coeffRef(3 * ph + 2, 3 * ph_ne + 2) = dtinv * dxdy * alpha * alpha;
                        //
                        if (true)
                        {
                            A.coeffRef(3 * ph + 2, 3 * ph_n) = 0.5 * theta * g * dx * htheta_jp12;
                            A.coeffRef(3 * ph + 2, 3 * ph) = 0.5 * theta * g * dx * htheta_jm12 - 0.5 * theta * g * dx * htheta_jp12;
                            A.coeffRef(3 * ph + 2, 3 * ph_s) = -0.5 * theta * g * dx * htheta_jm12;
                        }
                        else
                        {  // brackets only useful for editor
                            double fac = theta * 0.125 * 0.0625 * 0.5 * dx;
                            A.coeffRef(3 * ph + 2, 3 * ph) = 0.0;  // fac* (-3.0 + 9.0 + 9.0 - 3.0 - 9.0 + 3.0 + 3.0 - 9.0);  // 
                            A.coeffRef(3 * ph + 2, 3 * ph_sw) = fac * (-3.0 + 1.0);  // 
                            A.coeffRef(3 * ph + 2, 3 * ph_s) = fac * (-9.0 - 9.0 + 3.0 + 3.0); // 
                            A.coeffRef(3 * ph + 2, 3 * ph_se) = fac * (-3.0 + 1.0);  //
                            A.coeffRef(3 * ph + 2, 3 * ph_e) = fac * (-1.0 + 3.0 - 3.0 + 1.0);  // 
                            A.coeffRef(3 * ph + 2, 3 * ph_ne) = fac * (-1.0 + 3.0);  // 
                            A.coeffRef(3 * ph + 2, 3 * ph_n) = fac * (-3.0 - 3.0 + 9.0 + 9.0);  //   
                            A.coeffRef(3 * ph + 2, 3 * ph_nw) = fac * (-1.0 + 3.0);  // 
                            A.coeffRef(3 * ph + 2, 3 * ph_w) = fac * (-1.0 + 3.0 - 3.0 + 1.0);  //
                        }
                        // 
                        // RHS r-momentum equation
                        //
                        double depth_0 = c_scv(htheta_0, htheta_w, htheta_sw, htheta_s);
                        double depth_1 = c_scv(htheta_0, htheta_s, htheta_se, htheta_e);
                        double depth_2 = c_scv(htheta_0, htheta_e, htheta_ne, htheta_n);
                        double depth_3 = c_scv(htheta_0, htheta_n, htheta_nw, htheta_w);

                        double grad_0y = 0.25 * dxdy / dy * dcdy_scvf2(htheta_0, htheta_s, htheta_w, htheta_sw);
                        double grad_1y = 0.25 * dxdy / dy * dcdy_scvf2(htheta_0, htheta_s, htheta_e, htheta_se);
                        double grad_4y = 0.25 * dxdy / dy * dcdy_scvf2(htheta_n, htheta_0, htheta_ne, htheta_e);
                        double grad_5y = 0.25 * dxdy / dy * dcdy_scvf2(htheta_n, htheta_0, htheta_nw, htheta_w);

                        rhs[3 * ph + 2] = -dtinv * dxdy * alpha * alpha * (rp[ph_sw] - rn[ph_sw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (rp[ph_s] - rn[ph_s]) +
                            -dtinv * dxdy * alpha * alpha * (rp[ph_se] - rn[ph_se]) +
                            //
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (rp[ph_w] - rn[ph_w]) +
                            -dtinv * dxdy * (1.0 - 2. * alpha) * (1.0 - 2. * alpha) * (rp[ph] - rn[ph]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (rp[ph_e] - rn[ph_e]) +
                            //
                            -dtinv * dxdy * alpha * alpha * (rp[ph_nw] - rn[ph_nw]) +
                            -dtinv * dxdy * alpha * (1.0 - 2. * alpha) * (rp[ph_n] - rn[ph_n]) +
                            -dtinv * dxdy * alpha * alpha * (rp[ph_ne] - rn[ph_ne]);
                        //
                        double a = -g * 0.25 * (depth_0 + depth_1 + depth_2 + depth_3) * (grad_0y + grad_1y + grad_4y + grad_5y);
                        double b = -0.5 * g * dy * (htheta_jp12 * htheta_jp12 - htheta_jm12 * htheta_jm12);
                        rhs[3 * ph + 2] += b;
                        if (i == nx / 2 + 1 && j == ny / 2 + 1)
                        {
                            int c = 1;
                        }
                    }
                }
            }
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(Matrix initialization);
            }
            if (nst == 1 && iter == 0)
            {
                START_TIMER(BiCGStab initialization);
            }
            else if (nst != 1)
            {
                START_TIMER(BiCGStab);
            }
            solver.compute(A);
            solver.setTolerance(eps_bicgstab);
            solution = solver.solve(rhs);
            //solution = solver.solveWithGuess(rhs, solution);
            if (nst == 1 && iter == 0)
            {
                STOP_TIMER(BiCGStab initialization);
                //log_file << "=======================================================" << std::endl;
                //log_file << Eigen::MatrixXd(A) << std::endl;
                //log_file << "=======================================================" << std::endl;
                //log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
                //log_file << "=======================================================" << std::endl;
            }
            else if (nst != 1)
            {
                STOP_TIMER(BiCGStab);
            }
            //log_file << "time [sec]:" << dt * double(nst) << std::endl;
            //log_file << "    #iterations    :" << solver.iterations() << std::endl;
            //log_file << "    estimated error:" << solver.error() << std::endl;
            // 
            // The new solution is the previous iterant plus the delta
            dh_max = 0.0;
            dq_max = 0.0;
            dr_max = 0.0;
            dh_maxi = 0;
            dq_maxi = 0;
            dr_maxi = 0;
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    k = j * nx + i;
                    // needed for postprecessing
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
            used_iter = iter;
            used_lin_iter += solver.iterations();
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
                int a = 1;
            }
        }
        STOP_TIMER(Newton iteration);

        if (stationary) 
        {
            std::cout << "stationary solution " << std::endl;
            //log_file << "stationary solution " << std::endl;
            //log_file << std::setprecision(8) << std::scientific
            //    << "    Iter: " << used_iter + 1
            //    << "    Delta h^{n + 1,p + 1}: " << dh_max << " at: " << dh_maxi
            //    << "    Delta q^{n + 1,p + 1}: " << dq_max << " at: " << dq_maxi
            //    << "    Delta r^{n + 1,p + 1}: " << dr_max << " at: " << dr_maxi;
        }
        else
        {
            std::cout << std::fixed << std::setprecision(2) << tstart + dt * double(nst) << ";   " << tstart + tstop << std::endl;
            //log_file << std::setprecision(8) << std::scientific
            //    << tstart + dt * double(nst)
            //    << "    Iter: " << used_iter + 1
            //    << "    Delta h^{n + 1,p + 1}: " << dh_max
            //    << "    Delta q^{n + 1,p + 1}: " << dq_max
            //    << "    Delta r^{n + 1,p + 1}: " << dr_max;
         }
        if (used_iter + 1 == iter_max)
        {
            if (dh_max > eps_newton || dq_max > eps_newton || dr_max > eps_newton)
            {
                log_file << "    ----    maximum number of iterations reached, probably not convergenced" << std::endl;
            }
        }
        for (int k = 0; k < nxny; ++k)
        {
            hn[k] = hp[k]; // h, continuity-eq
            qn[k] = qp[k];  // q, momentum-eq
            rn[k] = rp[k];  // q, momentum-eq
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
            std::vector<double> his_values = { s[p_a], s[p_b], s[p_c], s[p_d], s[centre], s[n_bnd], s[ne_bnd], s[e_bnd], s[se_bnd], s[s_bnd], s[sw_bnd], s[w_bnd], s[nw_bnd] };
            his_file->put_variable(his_s_name, nst_his, his_values);

            his_values = { u[p_a], u[p_b], u[p_c], u[p_d], u[centre], u[n_bnd], u[ne_bnd], u[e_bnd], u[se_bnd], u[s_bnd], u[sw_bnd], u[w_bnd], u[nw_bnd] };
            his_file->put_variable(his_u_name, nst_his, his_values);

            his_values = { v[p_a], v[p_b], v[p_c], v[p_d], v[centre], v[n_bnd], v[ne_bnd], v[e_bnd], v[se_bnd], v[s_bnd], v[sw_bnd], v[w_bnd], v[nw_bnd] };
            his_file->put_variable(his_v_name, nst_his, his_values);

            his_values = { double(used_iter) };
            his_file->put_variable(his_newton_iter_name, nst_his, his_values);
            his_values = { double(solver.iterations()) };
            his_file->put_variable(his_bicgstab_iter_name, nst_his, his_values);
            his_values = { solver.error() };
            his_file->put_variable(his_bicgstab_iter_error_name, nst_his, his_values);
            STOP_TIMER(Writing his-file);
        }
    } // End of the time loop
    STOP_TIMER(Time loop);
    START_TIMER(Finalize);
    // Finalization
    status = his_file->close();
    status = map_file->close();
    log_file.close();
    STOP_TIMER(Finalize);
    STOP_TIMER(Main);
    PRINT_TIMER(timing_filename.data());

    return 0;
}
int p_index(int i, int j, int nx)
{
    return j * nx + i;
}
double c_scv(double h1, double h2, double h3, double h4)
{
    // value at subcontrol volume
    return 0.0625 * (9. * h1 + 3. *h2 +  h3 + 3* h4);
}
double dcdx_scvf1(double h1, double h2, double h3, double h4)
{
    // dcdx at subcontrol volume edge
    return 0.25 * (3. * h1 - 3. * h2 + h3 - h4);
}
double dcdy_scvf1(double h1, double h2, double h3, double h4)
{
    // dcdy at subcontrol volume edge
    return 0.5 * (h4 - h1 + h3 - h2);
}
double dcdy_scvf2(double h1, double h2, double h3, double h4)
{
    // dcdy at subcontrol volume edge
    return 0.25 * (3. * h1 - 3. * h2 + h3 - h4);
}
double dcdx_scvf2(double h1, double h2, double h3, double h4)
{
    // dcdy at subcontrol volume edge
    return 0.5 * (h2 - h1 + h3 - h4);
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

