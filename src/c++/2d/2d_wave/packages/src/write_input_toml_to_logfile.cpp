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

#include <string>
#include <filesystem>

#include "data_input_struct.h"
#include "toml.h"

std::string bool_value_as_string(bool value);
std::string format_as_double(double value);


int write_used_input(struct _data_input data, std::ofstream & log_file){

    log_file << "logging = \"" << data.log.logging << "\"  # \"iterations\", \"matrix\", \"pattern\"" << std::endl;

    // Boundary
    log_file << std::endl << "[Boundary]  # north, east, south, west" << std::endl;
    log_file << "    " << "treg         = " << format_as_double(data.boundary.treg) << "  # Regularization time boundary signal" << std::endl;
    log_file << "    " << "eps_bc_corr  = " << format_as_double(data.boundary.eps_bc_corr) << std::endl;

    log_file << "    " << "bc_type      = [\"";
    for (int i = 0; i < data.boundary.bc_type.size() - 1; ++i) { log_file << data.boundary.bc_type[i] << "\", \""; }
    log_file << data.boundary.bc_type[data.boundary.bc_type.size() - 1] << "\"]  # Type \"free_slip\", \"no_slip\", \"borsboom\", \"mooiman\"  " << std::endl;
    
    log_file << "    " << "bc_vars      = [\"";
    for (int i = 0; i < data.boundary.bc_vars.size() - 1; ++i) { log_file << data.boundary.bc_vars[i] << "\", \""; }
    log_file << data.boundary.bc_vars[data.boundary.bc_vars.size() - 1] << "\"]" << std::endl;

    log_file << "    " << "bc_vals      = [";
    for (int i = 0; i < data.boundary.bc_vals.size() - 1; ++i) { log_file << format_as_double(data.boundary.bc_vals[i]) << ", "; }
    log_file << format_as_double(data.boundary.bc_vals[data.boundary.bc_vals.size() - 1]) << "]" << std::endl;

    // Domain
    log_file << std::endl << "[Domain]" << std::endl;
    log_file << "    " << "bed_level_file = \"" << data.domain.bed_level_filename << "\"" << std::endl;
    log_file << "    " << "mesh_file      = \"" << data.domain.grid_filename << "\"" << std::endl;

    // Initial
    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "    " << "ini_vars      = [\"";
    for (int i = 0; i < data.initial.ini_vars.size() - 1; ++i) { log_file << data.initial.ini_vars[i]  << "\", \""; }
    log_file << data.initial.ini_vars[data.initial.ini_vars.size() - 1] << "\"]" << std::endl; 

    log_file << "    " << "gauss_amp     = " << format_as_double(data.initial.gauss_amp) << std::endl;
    log_file << "    " << "gauss_mu      = " << format_as_double(data.initial.gauss_mu) << std::endl;
    log_file << "    " << "gauss_mu_x    = " << format_as_double(data.initial.gauss_mu_x) << std::endl;
    log_file << "    " << "gauss_mu_y    = " << format_as_double(data.initial.gauss_mu_y) << std::endl;
    log_file << "    " << "gauss_sigma   = " << format_as_double(data.initial.gauss_sigma) << std::endl;
    log_file << "    " << "gauss_sigma_x = " << format_as_double(data.initial.gauss_sigma_x) << std::endl;
    log_file << "    " << "gauss_sigma_y = " << format_as_double(data.initial.gauss_sigma_y) << std::endl;

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "    " << "dt                  = " << format_as_double(data.numerics.dt) << "  # Time step size [s], if dt == 0: then stationary problem" << std::endl;
    log_file << "    " << "theta               = " << format_as_double(data.numerics.theta) << "  # Implicitness factor (0.5 <= theta <= 1.0)" << std::endl;
    log_file << "    " << "c_psi               = " << format_as_double(data.numerics.c_psi) << std::endl;
    log_file << "    " << "iter_max            = " << data.numerics.iter_max << "  # Maximum number of nonlinear iterations" << std::endl;
    log_file << "    " << "eps_newton          = " << format_as_double(data.numerics.eps_newton) << std::endl;
    log_file << "    " << "eps_bicgstab        = " << format_as_double(data.numerics.eps_bicgstab) << std::endl;
    log_file << "    " << "eps_abs_function    = " << format_as_double(data.numerics.eps_abs) << std::endl;
    log_file << "    " << "linear_solver       = \"" << data.numerics.linear_solver << "\"" << "  # \"bicgstab\", \"multigrid\"" << std::endl;
    log_file << "    " << "regularization_init = " << bool_value_as_string(data.numerics.regularization_init) << std::endl;
    log_file << "    " << "regularization_iter = " << bool_value_as_string(data.numerics.regularization_iter) << std::endl;
    log_file << "    " << "regularization_time = " << bool_value_as_string(data.numerics.regularization_time) << std::endl;

    // Physics
    log_file << std::endl << "[Physics]" << std::endl;
    log_file << "    " << "do_linear_waves     = " << bool_value_as_string(data.physics.do_linear_waves) << std::endl;
    log_file << "    " << "do_continuity       = " << bool_value_as_string(data.physics.do_continuity) << std::endl;
    log_file << "    " << "do_q_equation       = " << bool_value_as_string(data.physics.do_q_equation) << std::endl;
    log_file << "    " << "do_r_equation       = " << bool_value_as_string(data.physics.do_r_equation) << std::endl;
    log_file << "    " << "do_convection       = " << bool_value_as_string(data.physics.do_convection) << std::endl;
    log_file << "    " << "do_bed_shear_stress = " << bool_value_as_string(data.physics.do_bed_shear_stress) << std::endl;
    log_file << "    " << "do_viscosity        = " << bool_value_as_string(data.physics.do_viscosity) << std::endl;
    log_file << "    " << "chezy_coefficient   = " << format_as_double(data.physics.chezy_coefficient) << std::endl;
    log_file << "    " << "viscosity           = " << format_as_double(data.physics.visc_const) << std::endl;

    // Output
    log_file << std::endl << "[Output]" << std::endl;
    log_file << "    " << "dt_his = " << format_as_double(data.output.dt_his) << std::endl;
    log_file << "    " << "dt_map = " << format_as_double(data.output.dt_map) << std::endl;
    log_file << "    " << "dt_screen = " << format_as_double(data.output.dt_screen) << std::endl;

    // Time
    log_file << std::endl << "[Time]" << std::endl;
    log_file << "    " << "tstart = " << format_as_double(data.time.tstart) << std::endl;
    log_file << "    " << "tstop  = " << format_as_double(data.time.tstop) << std::endl;

    // Observation points
    log_file << std::endl;
    log_file << "# Hard coded observation points are not listed" << std::endl;
    log_file << std::endl;
    for (size_t i = 0; i < data.obs_points.size(); ++i)
    {
        log_file << "[[ObservationPoint]]" << std::endl;
        log_file << "    x    = " << data.obs_points[i].x    << std::endl;
        log_file << "    y    = " << data.obs_points[i].y    << std::endl;
        log_file << "    name = \"" << data.obs_points[i].name << "\"" << std::endl;
        if (i < data.obs_points.size() - 1) { log_file << std::endl; }
    }

    return 0;
}
std::string bool_value_as_string(bool value)
{
    std::string bool_string;
    if (value == 1) { bool_string = "true"; }
    else if (value ==0 ) { bool_string = "false"; }
    return bool_string;
}
std::string format_as_double(double value) {
    std::stringstream ss;

    // Check if the number is effectively an integer
    if (std::floor(value) == value) {
        // Force 1 digit after the decimal point
        ss << std::fixed << std::setprecision(1) << value;
    } else {
        // Print full precision (optional: limit to 6 digits)
        ss << std::defaultfloat << value;
    }

    return ss.str();
}