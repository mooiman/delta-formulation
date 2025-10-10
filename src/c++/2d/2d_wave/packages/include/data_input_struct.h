//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formuation and Modified Newton iteration 
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

#ifndef __DATA_INPUT_STRUCT_H__
#define __DATA_INPUT_STRUCT_H__

#include <string>
#include <filesystem>

// data structures needed to store input data from the master input file (ie toml-file)

struct _ObservationPoint {
    size_t idx;
    double x;
    double y;
    std::string name;
};
struct _log {
    std::string logging;
};
struct _boundary {
    double treg;
    double eps_bc_corr;
    std::vector<double> bc_vals;
    std::vector<std::string> bc_type;
    std::vector<std::string> bc_vars;
    std::vector<bool> bc_absorbing;
};
struct _domain {
    std::string grid_filename;
    std::filesystem::path full_grid_filename;
    std::string bed_level_filename;
    std::filesystem::path full_bed_level_filename;
};
struct _initial {
    double gauss_amp;
    double gauss_mu;
    double gauss_mu_x;
    double gauss_mu_y;
    double gauss_sigma;
    double gauss_sigma_x;
    double gauss_sigma_y;
    std::vector<std::string> ini_vars;
};
struct _numerics {
    double dt;
    double theta;
    double c_psi;
    int iter_max;
    double eps_newton;
    double eps_bicgstab;
    double eps_abs;  // epsilon needed to approximate the abs-function by a continues function
    bool regularization_init;
    bool regularization_iter;
    bool regularization_time;
    std::string linear_solver;
};
struct _physics {
    double visc_const;
    double g;
    bool do_continuity;
    bool do_q_equation;
    bool do_r_equation;
    bool do_convection;
    bool do_bed_shear_stress;
    bool do_linear_waves;
    double chezy_coefficient;
    bool do_viscosity;
};
struct _time {
    double tstart;
    double tstop;
};
struct _output {
    double dt_his;
    double dt_map;
};

struct _data_input{
    struct _log log;
    struct _boundary boundary;
    struct _domain domain;
    struct _initial initial;
    struct _numerics numerics;
    struct _physics physics; 
    struct _time time;
    struct _output output;
    std::vector<_ObservationPoint> obs_points;
};

#endif __DATA_INPUT_STRUCT_H__
