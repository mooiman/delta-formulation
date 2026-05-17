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
//   DESCRIPTION
//
#ifndef __OBSERVATION_STATIONSL_H__
#define __OBSERVATION_STATIONSL_H__

#include <algorithm>
#include <vector>

#include "data_input_struct.h"
#include <include/KDtree.hpp>

int def_observation_stations(std::vector<std::string>& obs_stations, std::vector<_ObservationPoint>& obs_points, KDTree xy_tree, std::vector<double> x, std::vector<double> y,
    size_t nx, size_t ny);

int add_hard_coded_observation_points(std::vector<_ObservationPoint>& obs_points, std::vector<double> x, std::vector<double> y, 
    size_t nx, size_t ny);
struct _ObservationPoint add_obs(size_t ptr_obs, double x_obs, double y_obs, std::string obs_name);
std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name);
std::string string_format_with_zeros(double value, int width);
inline size_t idx(size_t i, size_t j, size_t ny);

#endif  // __OBSERVATION_STATIONSL_H__