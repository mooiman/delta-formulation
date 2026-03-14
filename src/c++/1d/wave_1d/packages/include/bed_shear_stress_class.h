//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
//    Copyright (C) 2025 Jan Mooiman
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------------

#ifndef __BEDSHEARSTRESS_H__
#define __BEDSHEARSTRESS_H__

#include <vector>
#include <Eigen/Sparse>

class BEDSHEARSTRESS
{
public:
    BEDSHEARSTRESS();
    ~BEDSHEARSTRESS();
    BEDSHEARSTRESS(double theta, double dx, double cf_in, double eps_in, int nx);
    int matrix_rhs(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& rhs,
        std::vector<double>& hp, std::vector<double>& qp,
        std::vector<double>& hn, std::vector<double>& qn);
    void rhs(std::vector<double>& rhs_q, std::vector<double> hn, std::vector<double>& qn);

private:
    inline double scv(double c0, double c1);
    inline double J_10(double& h, double& q);  // right hand side q-equation
    inline double J_11(double& h, double& q);
    inline double J_12(double& h, double& q);
    inline double vecq(double& qp);

    double theta;
    double dx;
    double cf;
    double eps;
    int nx;

    double cv_area;
};
#endif  // __BEDSHEARSTRESS_H__
