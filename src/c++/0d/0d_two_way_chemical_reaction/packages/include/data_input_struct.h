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
struct _initial {
    double u0_ini;
    double u1_ini;
};
struct _numerics {
    double dt;
    double theta;
    int iter_max;
    double eps_newton;
    double eps_bicgstab;
    std::string time_integration;
};
struct _physics {
    double k1;
    double k2;
};
struct _time {
    std::string tunit;
    double tstart;
    double tstop;
};
struct _output {
    double dt_his;
    double dt_screen;
};

struct _data_input{
    struct _log log;
    struct _initial initial;
    struct _numerics numerics;
    struct _physics physics; 
    struct _time time;
    struct _output output;
};

#endif __DATA_INPUT_STRUCT_H__
