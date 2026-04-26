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

#include "observation_stations.h"

int  def_observation_stations(std::vector<std::string>& obs_stations, std::vector<_ObservationPoint>& obs_points, KDTree xy_tree, 
    std::vector<double> x, std::vector<double> y, size_t nx, size_t ny)
{
    int status;
    status = add_hard_coded_observation_points(obs_points, x, y, nx, ny);

    int nsig = 0;
    for (size_t i = 0; i < obs_points.size(); ++i)
    {
        nsig = std::max(nsig, std::max( (int)std::log10(obs_points[i].x), (int)std::log10(obs_points[i].x) ));
    }
    //nsig += 1;
    
    for (size_t i = 0; i < obs_points.size(); ++i)
    {
        obs_stations.push_back(setup_obs_name(obs_points[i].x, obs_points[i].y, nsig, obs_points[i].name));
    }
    
    // Find the index of the observation point
    for (size_t i = 0; i < obs_points.size(); ++i)
    {
        std::vector<double> point = {obs_points[i].x, obs_points[i].y};
        obs_points[i].idx = xy_tree.nearest_index(point);
    }
    return status;
}

int add_hard_coded_observation_points(std::vector<_ObservationPoint>& obs_points, std::vector<double> x, std::vector<double> y, size_t nx, size_t ny)
{
    // Initialize observation station locations
    size_t ptr_obs = idx(nx / 2, ny / 2, ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Centre"));

    // boundaries --------------------------------------------------------------
    ptr_obs  = idx(1     , ny / 2, ny);  // at boundary, not at virtual point
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "W"));

    ptr_obs  = idx(nx - 2, ny / 2, ny);  // at boundary, not at virtual point
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "E"));

    ptr_obs  = idx(nx / 2, 1     , ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "S"));

    ptr_obs  = idx(nx / 2, ny - 2, ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "N"));

    // corners -----------------------------------------------------------------
    ptr_obs = idx(1     , 1     , ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "SW"));

    ptr_obs = idx(nx - 2, ny - 2, ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "NE"));

    ptr_obs = idx(1     , ny - 2, ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "NW"));

    ptr_obs = idx(nx - 2, 1     , ny);
    obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "SE"));

    return 0;
}
struct _ObservationPoint add_obs(size_t ptr_obs, double x_obs, double y_obs, std::string obs_name)
{
    struct _ObservationPoint obs_point; 
    obs_point.name = obs_name;
    obs_point.idx = ptr_obs;
    obs_point.x = x_obs;
    obs_point.y = y_obs;
    return obs_point;
}


std::string setup_obs_name(double x_obs, double y_obs, int nsig, std::string obs_name)
{
    std::string ss_x;
    std::string ss_y;
    ss_x = string_format_with_zeros(x_obs, nsig + 3);
    ss_y = string_format_with_zeros(y_obs, nsig + 3);

    //ss_x << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << x_obs;
    //ss_y << std::setfill('0') << std::setw(nsig + 3) << std::fixed << std::setprecision(2) << y_obs;
    return (obs_name + " (" + ss_x + ", " + ss_y + ") ");
}
std::string string_format_with_zeros(double value, int width) 
{
    std::ostringstream oss;
    if (value < 0) {
        oss << '-';
        oss << std::setw(width - 1) << std::setfill('0') << std::fixed << std::setprecision(1) << -value;
    } else {
        oss << std::setw(width) << std::setfill('0') << std::fixed << std::setprecision(1) << value;
    }
    return oss.str();
}
inline size_t idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}
