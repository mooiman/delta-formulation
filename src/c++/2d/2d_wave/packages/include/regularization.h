//
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
#ifndef __2D_REGULARIZATION_H__
#define __2D_REGULARIZATION_H__

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

class REGULARIZATION
{
public:
    REGULARIZATION();
    REGULARIZATION(int);

    void given_function(int nx, int ny, 
        std::vector<double>& u_out, std::vector<double>& psi, 
        std::vector<double>& err, std::vector<double>& u_giv,
        double dx, double dy, double c_psi, bool without_err);
    void first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx);
private:
    std::unique_ptr<std::vector<double>> solve_eq8(int nx, int ny, double dx, double dy, double c_psi, std::vector<double> u0, std::vector<double> u0_xixi, std::vector<double> u0_etaeta);
    std::unique_ptr<std::vector<double>> solve_eq7(int nx, int ny, double dx, double dy, std::vector<double> psi, std::vector<double> u_giv);

    int m_iter_max;
    int p_index(int i, int j, int nx);
};

#endif __2D_REGULARIZATION_H__