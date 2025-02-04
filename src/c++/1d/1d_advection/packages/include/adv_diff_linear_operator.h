#ifndef __ADV_DIFF_LINEAR_OPERATOR_H__
#define __ADV_DIFF_LINEAR_OPERATOR_H__

#include <cstdlib>
#include <vector>

void adv_diff_linear_operator(std::vector<double>& la, std::vector<double>& lb, std::vector<double>& lc, const std::vector<double> u, const std::vector<double> eps, const double dx, const int);

#endif __ADV_DIFF_LINEAR_OPERATOR_H__
