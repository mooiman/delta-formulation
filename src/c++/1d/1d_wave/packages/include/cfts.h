//
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
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
    int set_attribute(std::string var, std::string att_name, int att_val);
    int set_attribute(std::string var, std::string att_name, std::string att_val);

};

#endif __CFTS_H__
