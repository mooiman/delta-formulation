#ifndef __CONVECTION_H__
#define __CONVECTION_H__

#include <vector>
#include <Eigen/Sparse>

class CONVECTION
{
public:
    CONVECTION();
    ~CONVECTION();
    CONVECTION(double theta, double dx, double dy, int nx, int ny);
    int matrix_2d_q_eq(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs, 
        std::vector<double> hp, std::vector<double> qp, std::vector<double> rp);
    int rhs_2d_q_eq(Eigen::VectorXd rhs);                          // RHS vector [h, q, r]^{n}

private:
    inline int p_index(int i, int j, int nx);
    void c_scvf_xi(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs, int eq, 
        double c0, double c1, double c2, double c3, int p0, int p1, int p2, int p3);
    void c_scvf_eta(Eigen::SparseMatrix<double> A, Eigen::VectorXd rhs, int c_eq,
        double c0, double c1, double c2, double c3, int p0, int p1, int p2, int p3);

    double theta;
    double dx;
    double dy;
    int nx;
    int ny;

    double hp_im12;
    double hn_im12;
    double hp_ip12;
    double hn_ip12;
    double htheta_ip12;
    double htheta_im12;

    double qp_im12;
    double qn_im12;
    double qp_ip12;
    double qn_ip12;
    double qtheta_ip12;
    double qtheta_im12;
};
#endif  // __CONVECTION_H__

