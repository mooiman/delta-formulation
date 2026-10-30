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
#include "read_input_toml_file.h"
#include "toml.h"

int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

_data_input read_toml_file(std::filesystem::path & input_dir, std::filesystem::path & toml_file_name){
    input_dir = input_dir;  // needed if a file is given as input variable
    toml::table tbl;
    toml::table tbl_chp;  // table for a chapter
    struct _data_input data{};

    tbl = toml::parse_file(toml_file_name.c_str());

    data.log.logging = tbl["Logging"].value_or("None");

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    data.initial.u0_ini = tbl_chp["O"].value_or(double(1.0));  // 0.d0 ! O Zuurstof-molekuul (Oxygen-molecule)
    data.initial.u1_ini = tbl_chp["NO"].value_or(double(10.));  // 1.3d08 ! NO Stikstof-oxide
    data.initial.u2_ini = tbl_chp["NO2"].value_or(double(1.0));  // 5.0d11 ! NO2 Stikstof-di-oxide
    data.initial.u3_ini = tbl_chp["O3"].value_or(double(1.0));  // 8.0d11 ! O3 = Ozon

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    data.numerics.dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary simulation
    data.numerics.time_integration = tbl_chp["time_integration"].value_or("None");
    data.numerics.iter_max = tbl_chp["iter_max"].value_or(int(50));
    data.numerics.theta = tbl_chp["theta"].value_or(double(0.501));  // Implicitness factor (0.5 <= theta <= 1.0)
    data.numerics.eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));  // stop criterium for Newton iteration
    data.numerics.eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));  // stop criterium for BiCGStab iteration
    
    // Output
    tbl_chp = *tbl["Output"].as_table();
    data.output.dt_his = tbl_chp["dt_his"].value_or(double(1.0));  // write interval to his-file
    data.output.dt_screen = tbl_chp["dt_screen"].value_or(double(60.0));  // time interval counter on screen 

    // Physics
    tbl_chp = *tbl["Physics"].as_table();
    data.physics.k1 = tbl_chp["k1"].value_or(double(1.0));  // as in Hundsdorfer and Verwer
    data.physics.k2 = tbl_chp["k2"].value_or(double(10.));  // as in Hundsdorfer and Verwer 
    data.physics.k3 = tbl_chp["k3"].value_or(double(1.0));  // as in Hundsdorfer and Verwer
    data.physics.sigma = tbl_chp["sigma"].value_or(double(1.0e-07));  // as in Hundsdorfer and Verwer

    // Time
    tbl_chp = *tbl["Time"].as_table();
    data.time.tunit = tbl_chp["tunit"].value_or("s");
    data.time.tstart = tbl_chp["tstart"].value_or(double(0.0));
    data.time.tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.0));

    return data;
}
