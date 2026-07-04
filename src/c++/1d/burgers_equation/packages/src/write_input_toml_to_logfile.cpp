//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
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

    log_file << "Logging = \"" << data.log.logging << "\"  # \"iterations\", \"matrix\", \"pattern\"" << std::endl;

    // Boundary
    log_file << std::endl << "[Boundary]" << "  # west, east" << std::endl;
    log_file << "    " << "treg        = " << format_as_double(data.boundary.treg) << "  # Regularization time boundary signal [s]" << std::endl;
    log_file << "    " << "eps_bc_corr = " << format_as_double(data.boundary.eps_bc_corr) << std::endl;

    log_file << "    " << "bc_signals  = [\"";
    for (int i = 0; i < data.boundary.bc_signals.size() - 1; ++i) { log_file << data.boundary.bc_signals[i] << "\", \""; }
    log_file << data.boundary.bc_signals[data.boundary.bc_signals.size() - 1] << "\"]"  << std::endl;

    log_file << "    " << "bc_type     = [\"";
    for (int i = 0; i < data.boundary.bc_type.size() - 1; ++i) { log_file << data.boundary.bc_type[i] << "\", \""; }
    log_file << data.boundary.bc_type[data.boundary.bc_type.size() - 1] << "\"]"  << std::endl;

    log_file << "    " << "bc_vals     = [";
    for (int i = 0; i < data.boundary.bc_vals.size() - 1; ++i) { log_file << format_as_double(data.boundary.bc_vals[i]) << ", "; }
    log_file << format_as_double(data.boundary.bc_vals[data.boundary.bc_vals.size() - 1]) << "]" << std::endl;

    //Domain
    log_file << std::endl << "[Domain]" << std::endl;
    log_file << "    Lx       = " << format_as_double(data.domain.Lx) << "  # Domain length [m]" << std::endl;
    log_file << "    x_origin = " << format_as_double(data.domain.x_origin) << "  # Origin of the x-coordinate [m]" <<  std::endl;

    // Initial
    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "    u_initial = " << format_as_double(data.initial.u_initial) << "  # initial velocity [m]" << std::endl;
    log_file << "    ini_vars  = [\"";
    for (int i = 0; i < data.initial.ini_vars.size() - 1; ++i) { log_file << data.initial.ini_vars[i] << "\", \""; }
    log_file << data.initial.ini_vars[data.initial.ini_vars.size() - 1] << "\"]" << std::endl;

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "    dt                  = " << format_as_double(data.numerics.dt) << "  # Time step size [s], if dt == 0: then stationary problem" << std::endl;
    log_file << "    dx                  = " << format_as_double(data.numerics.dx) << "  # Grid size [m]" << std::endl;
    log_file << "    c_psi               = " << format_as_double(data.numerics.c_psi) << std::endl;
    log_file << "    theta               = " << format_as_double(data.numerics.theta) << "  # Implicitness factor (0.5 <= theta <= 1.0)" << std::endl;
    log_file << "    iter_max            = " << format_as_double(data.numerics.iter_max) << "  # Maximum number of non-linear iterations" << std::endl;
    log_file << "    eps_newton          = " << format_as_double(data.numerics.eps_newton) << std::endl;
    log_file << "    eps_bicgstab        = " << format_as_double(data.numerics.eps_bicgstab) << std::endl;
    log_file << "    regularization_init = " << bool_value_as_string(data.numerics.regularization_init) << std::endl;
    log_file << "    regularization_iter = " << bool_value_as_string(data.numerics.regularization_iter) << std::endl;
    log_file << "    regularization_time = " << bool_value_as_string(data.numerics.regularization_time) << std::endl;

    // Output
    log_file << std::endl << "[Output]" << std::endl;
    log_file << "    dt_his    = " << format_as_double(data.output.dt_his) << std::endl;
    log_file << "    dt_map    = " << format_as_double(data.output.dt_map) << std::endl;
    log_file << "    dt_screen = " << format_as_double(data.output.dt_screen) << std::endl;

    // Physics
    log_file << std::endl << "[Physics]" << std::endl;
    log_file << "    do_convection = " << bool_value_as_string(data.physics.do_convection) << std::endl;
    log_file << "    do_viscosity  = " << bool_value_as_string(data.physics.do_viscosity) << std::endl;
    log_file << "    viscosity     = " << format_as_double(data.physics.visc_const) << std::endl;
    log_file << "    do_source     = " << bool_value_as_string(data.physics.do_source) << std::endl;
    log_file << "    src_type      = " << data.physics.src_type << std::endl;

    //Time
    log_file << std::endl << "[Time]" << std::endl;
    log_file << "    tunit  = \"" << "s" << "\"" << "  # \"s\", \"m\", \"h\", \"d\"" << std::endl;
    log_file << "    tstart = " << format_as_double(data.time.tstart) << std::endl;
    log_file << "    tstop  = " << format_as_double(data.time.tstop) << std::endl;

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