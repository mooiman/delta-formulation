//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D Burgers equation, fully implicit with delta-formulation and Modified Newton iteration 
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

    // Initial
    log_file << std::endl << "[Initial]" << std::endl;
    log_file << "    u0 = " << format_as_double(data.initial.u0_ini) << std::endl; 
    log_file << "    u1 = " << format_as_double(data.initial.u1_ini) << std::endl;

    // Numerics
    log_file << std::endl << "[Numerics]" << std::endl;
    log_file << "    dt               = " << format_as_double(data.numerics.dt) << "  # Time step size [s], if dt == 0: then stationary problem" << std::endl;
    log_file << "    time_integration = " << "\"" << data.numerics.time_integration <<"\"" << "  # \"delta_formulation\", \"runge_kutta_4\", \"euler_explicit\"" << std::endl;
    log_file << "    theta            = " << format_as_double(data.numerics.theta) << "  # Implicitness factor (0.5 <= theta <= 1.0)" << std::endl;
    log_file << "    iter_max         = " << format_as_double(data.numerics.iter_max) << "  # Maximum number of non-linear iterations" << std::endl;
    log_file << "    eps_newton       = " << format_as_double(data.numerics.eps_newton) << std::endl;
    log_file << "    eps_bicgstab     = " << format_as_double(data.numerics.eps_bicgstab) << std::endl;

    // Output
    log_file << std::endl << "[Output]" << std::endl;
    log_file << "    dt_his    = " << format_as_double(data.output.dt_his) << "  # [s]" << std::endl;
    log_file << "    dt_screen = " << format_as_double(data.output.dt_screen) << "  # [s]" << std::endl;

    // Physics
    log_file << std::endl << "[Physics]" << std::endl;
    log_file << "    k1 = " << format_as_double(data.physics.k1) << "  # Rate constant k1" << std::endl;
    log_file << "    k2 = " << format_as_double(data.physics.k2) << "  # Rate constant k2"<< std::endl;

    //Time
    log_file << std::endl << "[Time]" << std::endl;
    log_file << "    tunit  = \"" << "s" << "\"" << "  # \"s\", \"m\", \"h\", \"d\"" << std::endl;
    log_file << "    tstart = " << format_as_double(data.time.tstart) << std::endl;
    log_file << "    tstop  = " << format_as_double(data.time.tstop) << std::endl;

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