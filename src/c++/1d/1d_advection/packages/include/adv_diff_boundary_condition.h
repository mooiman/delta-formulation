#ifndef __ADV_DIFF_BOUNDARY_CONDITION_H__
#define __ADV_DIFF_BOUNDARY_CONDITION_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>
#include <cmath>

enum class BND_TYPE
{
    NONE = 0,
    Constant,
    Sine,
    NR_BND_TYPES
};

void adv_diff_boundary_condition(double&, double&, double&, double&, double&, double&, BND_TYPE bnd_type);

#endif __ADV_DIFF_BOUNDARY_CONDITION_H__