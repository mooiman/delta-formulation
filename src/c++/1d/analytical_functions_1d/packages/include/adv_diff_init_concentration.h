#ifndef __ADV_DIFF_INIT_CONCENTRATION_H__
#define __ADV_DIFF_INIT_CONCENTRATION_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>

enum class SHAPE_CONC
{
    NONE = 0,
    Constant,
    Envelope,
    Boundary_layer,
    NR_SHAPES
};

void adv_diff_init_concentration(double time, double u_const, double eps, std::vector<double>& x, double Lx, SHAPE_CONC shape, std::vector<double>& d);

#endif __ADV_DIFF_INIT_CONCENTRATION_H__
