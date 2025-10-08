//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formulation and Modified Newton iteration 
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

int  def_observation_stations(std::vector<std::string>& obs_stations, std::vector<_ObservationPoint>& obs_points, KDTree xy_tree, std::vector<double> x, std::vector<double> y,
    double Lx, double Ly, double dx, double dy, int nx, int ny)
{
    int status;
    status = add_hard_coded_observation_points(obs_points, x, y, Lx, Ly, dx, dy, nx, ny);

    int nsig = 0;
    for (int i = 0; i < obs_points.size(); ++i)
    {
        nsig = std::max(nsig, std::max( (int)std::log10(obs_points[i].x), (int)std::log10(obs_points[i].x) ));
    }
    //nsig += 1;
    
    for (int i = 0; i < obs_points.size(); ++i)
    {
        obs_stations.push_back(setup_obs_name(obs_points[i].x, obs_points[i].y, nsig, obs_points[i].name));
    }
    
    // Find the index of the observation point
    for (int i = 0; i < obs_points.size(); ++i)
    {
        std::vector<double> point = {obs_points[i].x, obs_points[i].y};
        obs_points[i].idx = xy_tree.nearest_index(point);
    }
    return status;
}

int add_hard_coded_observation_points(std::vector<_ObservationPoint>& obs_points, std::vector<double> x, std::vector<double> y, double Lx, double Ly, double dx, double dy, int nx, int ny)
{

    // Initialize observation station locations
    int ptr_obs = idx(nx / 2, ny / 2, ny);
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

    // halfway stations --------------------------------------------------------
    //ptr_obs  = idx(    nx / 4 + 1,     ny / 2    , ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway W"));
    //
    //ptr_obs  = idx(3 * nx / 4 - 1,     ny / 2    , ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway E"));
    //
    //ptr_obs  = idx(    nx / 2    ,     ny / 4 + 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway S"));
    //
    //ptr_obs  = idx(    nx / 2    , 3 * ny / 4 - 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway N"));
    //
    //ptr_obs = idx(    nx / 4 + 1,     ny / 4 + 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway SW"));
    //
    //ptr_obs = idx(3 * nx / 4 - 1, 3 * ny / 4 - 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway NE"));
    //
    //ptr_obs = idx(    nx / 4 + 1, 3 * ny / 4 - 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway NW"));
    //
    //ptr_obs = idx(3 * nx / 4 - 1,     ny / 4 + 1, ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "Halfway SE"));

    // stations on equal distance of centre (ie 2500.0 [m]) ---------------------
    //double x_a = std::min(Lx / 2., 2500.);
    //double x_b = std::min(Lx / 2., 2000.);
    //double x_c = std::min(Lx / 2., 1500.);
    //double x_d = std::min(Lx / 2., 0.);
    //double y_a = std::min(Ly / 2., 0.);
    //double y_b = std::min(Ly / 2., 1500.);
    //double y_c = std::min(Ly / 2., 2000.);
    //double y_d = std::min(Ly / 2., 2500.);
    //
    //ptr_obs = idx(int((y_a / dy) + nx / 2), int((x_a / dx) + ny / 2), ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "D"));
    //
    //ptr_obs = idx(int((y_b / dy) + nx / 2), int((x_b / dx) + ny / 2), ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "C"));
    //
    //ptr_obs = idx(int((y_c / dy) + nx / 2), int((x_c / dx) + ny / 2), ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "B"));
    //
    //ptr_obs = idx(int((y_d / dy) + nx / 2), int((x_d / dx) + ny / 2), ny);
    //obs_points.push_back(add_obs(ptr_obs, x[ptr_obs], y[ptr_obs], "A"));

    return 0;
}
struct _ObservationPoint add_obs(int ptr_obs, double x_obs, double y_obs, std::string obs_name)
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
inline int idx(int i, int j, int ny)
{
    return i * ny + j;
}
