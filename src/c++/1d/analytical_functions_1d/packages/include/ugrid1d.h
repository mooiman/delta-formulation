//
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
#ifndef __UGRID1D_H__
#define __UGRID1D_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <format>

class UGRID1D
{
public:
    UGRID1D();
    ~UGRID1D();
    int open(std::string, std::string);
    int close();
    int network();
    int network_geometry();
    int mesh1d();

    int def_dimensions(int, int, int, int, int);
    int add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit);
    int add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit, std::string mesh, std::string location);
    int add_ntw_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string long_name);
    int add_geom_coordinate(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string std_name, std::string long_name, std::string unit);
    int add_node_count(std::string var_name, std::vector<std::string> dim_names, std::string long_name);
    int add_mesh1d_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string long_name);
    int add_mesh1d_point_on_branch(std::string var_name, std::vector<std::string> dim_names, std::string long_name);
    int add_mesh1d_offset_on_branch(std::string var_name, std::vector<std::string> dim_names, std::string long_name, std::string unit);
    int add_time_series(void);
    int put_time(int i, double time);
    int put_time_variable(std::string var_name, int i_time, std::vector<double> values);
    int put_variable(std::string var_name, std::vector<double> values);
    int put_variable(std::string var_name, std::vector<int> values);
    int put_variable_2(std::string var_name, std::vector<int> values);

private:
    int m_ncid;
    int m_times;
    std::string m_time_units;

    int set_global_attribute(std::string name, std::string value);
    int set_attribute(std::string var, std::string att_name, double att_val);
    int set_attribute(std::string var, std::string att_name, int att_val);
    int set_attribute(std::string var, std::string att_name, std::string att_val);

};
#endif __UGRID1D_H__
