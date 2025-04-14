#ifndef __COMPATIBLE_FUNCTION_H__
#define __COMPATIBLE_FUNCTION_H__

#define _USE_MATH_DEFINES
#include <cstdlib>
#include <vector>

// for bicgstab  solver
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

void compatible_function(std::vector<double>& x, std::vector<double>& u_out, std::vector<double>& u_in, std::vector<double>& mass);

#endif __COMPATIBLE_FUNCTION_H__
