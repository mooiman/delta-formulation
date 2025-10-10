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
#include "read_input_toml_file.h"
#include "toml.h"

int get_toml_array(toml::table, std::string, std::vector<std::string>&);
int get_toml_array(toml::table, std::string, std::vector<double>&);
int get_toml_array(toml::table, std::string, std::vector<bool>&);

_data_input read_toml_file(std::filesystem::path & input_dir, std::filesystem::path & toml_file_name){
    int status;
    toml::table tbl;
    toml::table tbl_chp;  // table for a chapter
    struct _data_input data{};

    tbl = toml::parse_file(toml_file_name.c_str());

    data.log.logging = tbl["Logging"].value_or("None");

    // Boundary
    tbl_chp = *tbl["Boundary"].as_table();
    data.boundary.eps_bc_corr = tbl_chp["eps_bc_corr"].value_or(double(0.0001));  // default 1e-4
    data.boundary.treg= tbl_chp["treg"].value_or(double(150.0));

    status = get_toml_array(tbl_chp, "bc_type", data.boundary.bc_type);
    status = get_toml_array(tbl_chp, "bc_vars", data.boundary.bc_vars);
    status = get_toml_array(tbl_chp, "bc_vals", data.boundary.bc_vals);
    status = get_toml_array(tbl_chp, "bc_absorbing", data.boundary.bc_absorbing);

    // Domain
    tbl_chp = *tbl["Domain"].as_table();

    std::string grid_filename = tbl_chp["mesh_file"].value_or("--none--");
    std::filesystem::path full_grid_filename = input_dir;
    full_grid_filename += grid_filename;
    data.domain.grid_filename = grid_filename;
    data.domain.full_grid_filename = full_grid_filename;

    std::string bed_level_filename = tbl_chp["bed_level_file"].value_or("--none--");
    std::filesystem::path full_bed_level_filename = input_dir;
    full_bed_level_filename += bed_level_filename;
    data.domain.bed_level_filename = bed_level_filename;
    data.domain.full_bed_level_filename = full_bed_level_filename;

    // Initial
    tbl_chp = *tbl["Initial"].as_table();
    status = get_toml_array(tbl_chp, "ini_vars", data.initial.ini_vars);
    data.initial.gauss_amp  = tbl_chp["gauss_amp"].value_or(double(0.0));   // amplitude of the gaussian hump at the boundary
    data.initial.gauss_mu   = tbl_chp["gauss_mu"].value_or(double(-INFINITY));
    data.initial.gauss_mu_x = tbl_chp["gauss_mu_x"].value_or(double(0.0));
    data.initial.gauss_mu_y = tbl_chp["gauss_mu_y"].value_or(double(0.0));
    data.initial.gauss_sigma = tbl_chp["gauss_sigma"].value_or(double(-INFINITY));
    data.initial.gauss_sigma_x = tbl_chp["gauss_sigma_x"].value_or(double(1.0));
    data.initial.gauss_sigma_y = tbl_chp["gauss_sigma_y"].value_or(double(1.0));

    // Numerics
    tbl_chp = *tbl["Numerics"].as_table();
    data.numerics.dt       = tbl_chp["dt"].value_or(double(0.0));  // dt == 0 => default stationary
    data.numerics.theta    = tbl_chp["theta"].value_or(double(0.501));
    data.numerics.c_psi    = tbl_chp["c_psi"].value_or(double(4.));
    data.numerics.iter_max = tbl_chp["iter_max"].value_or(int(50));
    data.numerics.eps_newton    = tbl_chp["eps_newton"].value_or(double(1.0e-12));
    data.numerics.eps_bicgstab  = tbl_chp["eps_bicgstab"].value_or(double(1.0e-12));
    data.numerics.eps_abs       = tbl_chp["eps_absolute"].value_or(double(1.0e-2));  // epsilon needed to approximate the abs-function by a continues function
    data.numerics.linear_solver = tbl_chp["linear_solver"].value_or("bicgstab");
    data.numerics.regularization_init = tbl_chp["regularization_init"].value_or(bool(false));
    data.numerics.regularization_iter = tbl_chp["regularization_iter"].value_or(bool(false));
    data.numerics.regularization_time = tbl_chp["regularization_time"].value_or(bool(false));

    //Physics
    tbl_chp = *tbl["Physics"].as_table();
    data.physics.g = tbl_chp["g"].value_or(double(9.81));  // Gravitational acceleration
    data.physics.do_continuity   = tbl_chp["do_continuity"].value_or(bool(true));  // default, continuity
    data.physics.do_q_equation   = tbl_chp["do_q_equation"].value_or(bool(true));  // default, q_equation
    data.physics.do_r_equation   = tbl_chp["do_r_equation"].value_or(bool(true));  // default, r_equation
    data.physics.do_convection   = tbl_chp["do_convection"].value_or(bool(false));  // default, no convection
    data.physics.do_linear_waves = tbl_chp["do_linear_waves"].value_or(bool(true));  // default, do linear wave simulation
    
    data.physics.do_viscosity = tbl_chp["do_viscosity"].value_or(bool(false));  // default, no viscosity
    data.physics.visc_const = tbl_chp["viscosity"].value_or(double(0.0001));  // default 1e-4

    data.physics.do_bed_shear_stress = tbl_chp["do_bed_shear_stress"].value_or(bool(false));  // default, no bed shear stress
    data.physics.chezy_coefficient = tbl_chp["chezy_coefficient"].value_or(double(50.0));

    // Output
    tbl_chp = *tbl["Output"].as_table();
    data.output.dt_his = tbl_chp["dt_his"].value_or(double(1.0));  // write interval to his-file
    data.output.dt_map = tbl_chp["dt_map"].value_or(double(0.0));  // write interval to his-file

    // Time
    tbl_chp = *tbl["Time"].as_table();
    data.time.tstart = tbl_chp["tstart"].value_or(double(0.0));
    data.time.tstop = tbl_chp["tstop"].value_or(double(60.));

    // Observation points
    data.obs_points = read_observation_points(tbl);

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
