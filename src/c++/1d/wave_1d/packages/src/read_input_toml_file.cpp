//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D wave equation, fully implicit with delta-formulation and Modified Newton iteration 
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
    int status;
    toml::table tbl;
    toml::table tbl_chp;  // table for a chapter
    struct _data_input data{};

    tbl = toml::parse_file(toml_file_name.c_str());

    std::string logging = tbl["Logging"].value_or("None");

    // Domain
    tbl_chp = *tbl["Domain"].as_table();
    data.domain.Lx = tbl_chp["Lx"].value_or(double(12000.));
    data.domain.depth = tbl_chp["depth"].value_or(double(10.));
    data.domain.geometry_type = tbl_chp["geometry_type"].value_or("None");

    // Time
    tbl_chp = *tbl["Time"].as_table();
    data.time.tstart = tbl_chp["tstart"].value_or(double(0.0));
    data.time.tstop = tbl_chp["tstop"].value_or(double(2.0 * 3600.));

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    status = get_toml_array(tbl_chp, "ini_vars", data.initial.ini_vars);
    data.initial.gauss_mu = tbl_chp["gauss_mu"].value_or(double(0.0));  // location of th gaussian hump
    data.initial.gauss_sigma = tbl_chp["gauss_sigma"].value_or(-INFINITY);  // sigma of the Guassian hump
    data.initial.gauss_amp = tbl_chp["gauss_amp"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary

    // Boundary
    tbl_chp = *tbl["Boundary"].as_table();
    data.boundary.eps_bc_corr = tbl_chp["eps_bc_corr"].value_or(double(0.0001));  // default 1e-4
    data.boundary.treg = tbl_chp["treg"].value_or(double(150.));  // Regularization time to reach constant boundary value

    status = get_toml_array(tbl_chp, "bc_type", data.boundary.bc_type);
    status = get_toml_array(tbl_chp, "bc_vars", data.boundary.bc_vars);
    status = get_toml_array(tbl_chp, "bc_vals", data.boundary.bc_vals);

    // Physics
    tbl_chp = *tbl["Physics"].as_table();
    data.physics.g = tbl_chp["g"].value_or(double(10.));  // Gravitational acceleration
    data.physics.do_q_equation = tbl_chp["do_q_equation"].value_or(bool(false));
    data.physics.do_convection = tbl_chp["do_convection"].value_or(bool(false));
    data.physics.do_bed_shear_stress = tbl_chp["do_bed_shear_stress"].value_or(bool(false));
    data.physics.do_viscosity = tbl_chp["do_viscosity"].value_or(bool(false));
    data.physics.chezy_coefficient = tbl_chp["chezy_coefficient"].value_or(double(50.0));
    data.physics.visc_const = tbl_chp["viscosity"].value_or(double(0.0));

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    data.numerics.dt = tbl_chp["dt"].value_or(double(0.0));  // default stationary simulation
    data.numerics.dx = tbl_chp["dx"].value_or(double(10.));
    data.numerics.c_psi = tbl_chp["c_psi"].value_or(double(4.));
    data.numerics.iter_max = tbl_chp["iter_max"].value_or(int(50));
    data.numerics.theta = tbl_chp["theta"].value_or(double(0.501));  // Implicitness factor (0.5 <= theta <= 1.0)
    data.numerics.eps_newton = tbl_chp["eps_newton"].value_or(double(1.0e-12));  // stop criterium for Newton iteration
    data.numerics.eps_bicgstab = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));  // stop criterium for BiCGStab iteratio
    data.numerics.eps_abs = tbl_chp["eps_abs"].value_or(double(1.0e-2));  // epsilon needed to approximate the abs-function by a continues function
    data.numerics.regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    data.numerics.regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    data.numerics.regularization_time = tbl_chp["regularization_time"].value_or(bool(false));

    // Output
    tbl_chp = *tbl["Output"].as_table();
    data.output.dt_his = tbl_chp["dt_his"].value_or(double(1.0));  // write interval to his-file
    data.output.dt_map = tbl_chp["dt_map"].value_or(double(0.0));  // write interval to his-file
    data.output.dt_screen = tbl_chp["dt_screen"].value_or(double(60.0));  // time interval counter on screen 

    return data;
}
std::vector<_ObservationPoint> read_observation_points(toml::table tbl)
{
    auto obs_points = tbl["ObservationPoint"].as_array();
    std::vector<_ObservationPoint> points;
    if (obs_points == nullptr) { return points; }

    for (const auto& item : *obs_points) 
    {
        _ObservationPoint pt;
        const auto& table = *item.as_table();
        pt.x = table["x"].value_or(0.0);
        pt.y = table["y"].value_or(0.0);
        pt.name = table["name"].value_or("unknown");
        points.push_back(pt);
    }

    return points;
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
