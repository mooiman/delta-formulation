#ifndef __INITIAL_CONDITIONS_H__
#define __INITIAL_CONDITIONS_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>

void initial_conditions(std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny, 
    std::vector<double>& hn, std::vector<double>& qn, std::vector<double>& rn,
    std::vector<double>& hp, std::vector<double>& qp, std::vector<double>& rp, std::vector<double>& zb_giv,
    std::vector<std::string>& ini_vars, double gauss_amp, double gauss_mu, double gauss_sigma);

#endif __INITIAL_CONDITIONS_H__
