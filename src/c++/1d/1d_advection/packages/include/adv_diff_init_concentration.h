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
    NR_SHAPES
};

void adv_diff_init_concentration(std::vector<double>& mass, std::vector<double>& x, double Lx, SHAPE_CONC shape, std::vector<double>& d);
void control_volumes(std::vector<double>& u_ana, std::vector<double>& cv, double dx, size_t refine);
void compatible_function(std::vector<double>& mass, std::vector<double>& cv, std::vector<double>& u_ana, std::vector<double>& u_out, 
    double dx, size_t refine);


#endif __ADV_DIFF_INIT_CONCENTRATION_H__
