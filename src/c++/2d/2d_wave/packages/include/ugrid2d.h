//
// Programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
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

#ifndef __UGRID2D_H__
#define __UGRID2D_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <format>

class UGRID2D
{
public:
    UGRID2D();
    ~UGRID2D();
    int open(std::string, std::string);
    int close();
    int mesh2d();

    int def_dimensions(int, int, int, int);
    int add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit);
    int add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit, std::string mesh, std::string location);
    int add_attribute(std::string var_name, std::string att_name, std::string att_value);
    int add_attribute(std::string var_name, std::string att_name, double);
    int add_ntw_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string long_name);
    int add_geom_coordinate(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string std_name, std::string long_name, std::string unit);
    int add_node_count(std::string var_name, std::vector<std::string> dim_names, std::string long_name);
    int add_mesh2d_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string long_name);
    int add_time_series(void);
    int put_time(int i, double time);
    int put_time_variable(std::string var_name, int i_time, std::vector<double> values);
    int put_variable(std::string var_name, std::vector<double> values);
    int put_variable(std::string var_name, std::vector<int> values);
    int put_variable_2(std::string var_name, std::vector<int> values);
    int put_variable_4(std::string var_name, std::vector<int> values);

private:
    int m_ncid;
    int m_times;
    std::string m_time_units;

    int set_global_attribute(std::string name, std::string value);
    int set_attribute(std::string var, std::string att_name, double att_val);
    int set_attribute(std::string var, std::string att_name, int att_val);
    int set_attribute(std::string var, std::string att_name, std::string att_val);

};
#endif __UGRID2D_H__
