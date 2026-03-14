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

#ifndef __CFTS_H__
#define __CFTS_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <format>

class CFTS
{
public:
    CFTS();
    ~CFTS();
    int open(std::string, std::string);
    int close();
    int add_stations(std::vector<std::string>, std::vector<double>, std::vector<double>);
    int add_time_series(void);
    int add_variable(std::string var_name, std::string std_name, std::string long_name, std::string unit);
    int add_variable_without_location(std::string var_name, std::string std_name, std::string long_name, std::string unit);
    int put_time(int i, double time);
    int put_variable(std::string var_name, int i_time, std::vector<double> values);

private:
    int m_ncid;
    int m_times;
    std::string m_time_units;

    int set_global_attribute(std::string name, std::string value);
    int set_attribute(std::string var, std::string att_name, double att_val);
    int set_attribute(std::string var, std::string att_name, int att_val);
    int set_attribute(std::string var, std::string att_name, std::string att_val);
};

#endif __CFTS_H__
