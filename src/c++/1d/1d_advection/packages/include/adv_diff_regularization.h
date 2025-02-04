#ifndef __ADV_DIFF_REGULARIZATION__
#define __ADV_DIFF_REGULARIZATION__
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

void adv_diff_regularization(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& u_giv, double dx, double c_psi);

#endif __ADV_DIFF_REGULARIZATION__