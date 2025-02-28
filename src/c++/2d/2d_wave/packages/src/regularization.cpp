//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2024-11-29
//   email:      jan.mooiman@outlook.com
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Regularization of given function in 2 dimension
//
//   u - d\dx psi d\dx u - d\dy psi d\dy u = u_giv
// 
//   Based in article:
//       From M. Borsboom
//       Applied Numerical Mathematics 26 (1998)
//       Borsboom_developerrorminadaptgridmethod_ApplNumerMath1998.pdf
//

#include "regularization.h"

REGULARIZATION::REGULARIZATION()
{
    m_iter_max = 100;
}
REGULARIZATION::REGULARIZATION(int iter_max) :
    m_iter_max(iter_max)
{
}

void REGULARIZATION::given_function(int nx, int ny, 
    std::vector<double>& u_out, std::vector<double>& psi, 
    std::vector<double>& eq8, std::vector<double>& u_giv_in, 
    double dx, double dy, double c_psi, bool use_eq8)
{
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    int nxny = nx * ny;
    std::vector<double> u_giv(nxny, 0.);
    std::vector<double> u0(nxny, 0.);
    std::vector<double> u1(nxny, 0.);
    std::vector<double> u0_xixi(nxny, 0.);
    std::vector<double> u0_etaeta(nxny, 0.);
    std::vector<double> tmp(nxny, 0.);

    const auto [u_giv_in_min, u_giv_in_max] = std::minmax_element(u_giv_in.begin(), u_giv_in.end());
    double min_range = 0.0001;
    double u_giv_range = *u_giv_in_max - *u_giv_in_min > min_range ? *u_giv_in_max - *u_giv_in_min : min_range;

    for (int i = 0; i < nxny; ++i)
    {
        // scaling of the input array u_giv_in
        u0[i] = u_giv_in[i] / u_giv_range;
        u_giv[i] = u0[i];
    }
    for (int iter = 0; iter < m_iter_max; ++iter)
    {
        //std::cout << "Iteration: " << iter << std::endl;
        //std::cout << "Compute second derivative" << std::endl;
        double u0xixi_max = 0.0;
        double u0etaeta_max = 0.0;
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                u0_xixi[p(i, j, nx)]   = std::abs((u0[p(i + 1, j, nx)] - 2. * u0[p(i, j, nx)] + u0[p(i - 1, j, nx)]));
                u0_etaeta[p(i, j, nx)] = std::abs((u0[p(i, j + 1, nx)] - 2. * u0[p(i, j, nx)] + u0[p(i, j - 1, nx)]));
                u0xixi_max = std::max(u0xixi_max, std::abs(u0_xixi[i]));
                u0etaeta_max = std::max(u0etaeta_max, std::abs(u0_etaeta[i]));
            }
        }
        for (int i = 0; i < nx; ++i)  // horizontal direction
        {
            int j = 0;
            u0_etaeta[p(i, j, nx)] = u0_etaeta[p(i, j+2, nx)];
            u0_etaeta[p(i, j+1, nx)] = u0_etaeta[p(i, j+2, nx)];
            j = nx - 1;
            u0_etaeta[p(i, j, nx)] = u0_etaeta[p(i, j-2, nx)];
            u0_etaeta[p(i, j-1, nx)] = u0_etaeta[p(i, j-2, nx)];
        }
        for (int j = 0; j < ny; ++j)  // vertical direction
        {
            int i = 0;
            u0_xixi[p(i, j, nx)] = u0_xixi[p(i+2, j, nx)];
            u0_xixi[p(i+1, j, nx)] = u0_xixi[p(i+2, j, nx)];
            i = nx - 1;
            u0_xixi[p(i, j, nx)] = u0_xixi[p(i-2, j, nx)];
            u0_xixi[p(i-1, j, nx)] = u0_xixi[p(i-2, j, nx)];
        }


//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(dx, dy, c_psi, u0, u0_xixi, u0_etaeta));
//------------------------------------------------------------------------------
        if (use_eq8)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    psi[p(i, j, nx)] = c_psi * dx * dx * eq8[p(i, j, nx)];
                }
            }
        }
        else
        {
            for (int i = 0; i < nxny; ++i)
            {
                psi[i] = c_psi * dx * dx;
            }
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(dx, dy, psi, u_giv));
//------------------------------------------------------------------------------

        diff_max1 = 0.0;
        for (int i = 0; i < nxny; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(u0[i] - u_giv[i]));
        }
        if (std::abs(diff_max1 - diff_max0) > 1e-12)
        {
            diff_max0 = diff_max1;
        }
        else
        {
            break;
        }
    }
    for (int i = 0; i < nxny; ++i)
    {
        u_out[i] = u0[i] * u_giv_range;
    }
}

void REGULARIZATION::first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx)
{
    int nx = (int)eps.size();
    std::vector<double> pe(nx, 0.0);
    double pe_threshold = 1.8;

    for (int i = 0; i < nx; ++i)
    {
        pe[i] = std::abs(u[i] * dx / eps[i]);
    }
    for (int i = 0; i < nx; ++i)
    {
        psi[i] = 1.0;
        if (pe[i] > pe_threshold)
        {
            psi[i] = pe[i] / pe_threshold;
        }
    }
    int a = 1;
}
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(double dx, std::vector<double> psi, std::vector<double> u_giv)
{
    int nx = psi.size();
    auto u = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nx, 0.0);

    Eigen::SparseMatrix<double> B(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    for (int i = 1; i < nx - 2; ++i)
    {
        double psi_im12 = 0.5 * (psi[i] + psi[i - 1]);
        double psi_ip12 = 0.5 * (psi[i + 1] + psi[i]);
        //psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
        //psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

        //B.coeffRef(i, i - 2) = 0.0;
        B.coeffRef(i, i - 1) = dx * 1. / 8. - psi_im12 / dx;
        B.coeffRef(i, i) = dx * 6. / 8. + psi_ip12 / dx + psi_im12 / dx;
        B.coeffRef(i, i + 1) = dx * 1. / 8. - psi_ip12 / dx;
        //B.coeffRef(i, i + 2) = 0.0;
    }
    for (int i = 1; i < nx - 1; ++i)
    {
        tmp[i] = dx * (1. / 8. * u_giv[i - 1] + 6. / 8. * u_giv[i] + 1. / 8. * u_giv[i + 1]);
        rhs[i] = tmp[i];
    }
    int i = 0;
    B.coeffRef(i, i) = 1.0;
    B.coeffRef(i, i + 1) = 0.0;
    B.coeffRef(i, i + 2) = 0.0;
    tmp[i] = u_giv[i];
    rhs[i] = tmp[i];
    B.coeffRef(i + 1, i) = 0.0;
    B.coeffRef(i + 1, i + 1) = 1.0;
    B.coeffRef(i + 1, i + 2) = 0.0;
    tmp[i + 1] = u_giv[i + 1];
    rhs[i + 1] = tmp[i + 1];

    i = nx - 1;
    B.coeffRef(i, i - 2) = 0.0;
    B.coeffRef(i, i - 1) = 0.0;
    B.coeffRef(i, i) = 1.0;
    tmp[i] = u_giv[i];
    rhs[i] = tmp[i];
    B.coeffRef(i - 1, i - 2) = 0.0;
    B.coeffRef(i - 1, i - 1) = 1.0;
    B.coeffRef(i - 1, i) = 0.0;
    tmp[i - 1] = u_giv[i - 1];
    rhs[i - 1] = tmp[i - 1];

    //for (int i = 0; i < nx; ++i)
    //{
    //    solution[i] = u0[i];
    //}

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverB;
    solverB.compute(B);
    solution = solverB.solve(rhs);
    //solution = solverB.solveWithGuess(rhs, solution);

    for (int i = 0; i < nx; ++i)
    {
        u->push_back(solution[i]);
    }

    return u;
};


std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(int nx, int ny, double dx, double dy, double c_error, 
std::vector<double> u0, std::vector<double> u0_xixi, std::vector<double> u0_etaeta)
{
    int nxny = u0_xixi.size();
    auto err = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nxny, 0.0);

    Eigen::SparseMatrix<double> A(nxny, ny);
    Eigen::VectorXd solution(nxny);           // solution vector u
    Eigen::VectorXd rhs(nxny);                // RHS

    double alpha = 0.125;
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            A.coeffRef(i-1, j - 1) = alpha * alpha;
            A.coeffRef(i  , j - 1) = alpha * (1 - 2. * alpha) - c_error;
            A.coeffRef(i+1, j - 1) = alpha * alpha;
            A.coeffRef(i-1, j) = alpha * (1 - 2. * alpha) - c_error;
            A.coeffRef(i  , j) = (1 - 2. * alpha) * (1 - 2. * alpha) + 4. * c_error;
            A.coeffRef(i+1, j) = alpha * (1 - 2. * alpha) - c_error;
            A.coeffRef(i-1, j + 1) = alpha * alpha;
            A.coeffRef(i  , j + 1) = alpha * (1 - 2. * alpha) - c_error;
            A.coeffRef(i+1, j + 1) = alpha * alpha;
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        A.coeffRef(i  , j    ) = 1.0;
        A.coeffRef(i  , j + 1) = 0.0;
        A.coeffRef(i  , j + 2) = 0.0;
        j = ny - 1;
        A.coeffRef(i  , j    ) = 1.0;
        A.coeffRef(i  , j - 1) = 0.0;
        A.coeffRef(i  , j - 2) = 0.0;
    }
    for (int j = 0; j < ny; ++j)  // vertical direction
    {
        int i = 0;
        A.coeffRef(i    , j) = 1.0;
        A.coeffRef(i + 1, j) = 0.0;
        A.coeffRef(i + 2, j) = 0.0;
        i = nx - 1;
        A.coeffRef(i    , j) = 1.0;
        A.coeffRef(i - 1, j) = 0.0;
        A.coeffRef(i - 2, j) = 0.0;
    }
    
    
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            rhs[p(i, j, nx)] = alpha * alpha * (dx * u0_xixi[p(i - 1, j - 1, nx)] + dy * u0_etaeta[p(i - 1, j - 1, nx)]);
            rhs[p(i, j, nx)] += alpha * (1 - 2. * alpha) * (dx * u0_xixi[p(i - 1, j, nx)] + dy * u0_etaeta[p(i - 1, j, nx)]);
            rhs[p(i, j, nx)] += alpha * alpha * (dx * u0_xixi[p(i - 1, j + 1, nx)] + dy * u0_etaeta[p(i - 1, j + 1, nx)]);
            rhs[p(i, j, nx)] += alpha * (1 - 2. * alpha) * (dx * u0_xixi[p(i, j - 1, nx)] + dy * u0_etaeta[p(i, j - 1, nx)]);
            rhs[p(i, j, nx)] += (1 - 2. * alpha) * (1 - 2. * alpha) * (dx * u0_xixi[p(i, j, nx)] + dy * u0_etaeta[p(i, j, nx)]);
            rhs[p(i, j, nx)] += alpha * (1 - 2. * alpha) * (dx * u0_xixi[p(i, j + 1, nx)] + dy * u0_etaeta[p(i, j + 1, nx)]);
            rhs[p(i, j, nx)] += alpha * alpha * (dx * u0_xixi[p(i + 1, j - 1, nx)] + dy * u0_etaeta[p(i + 1, j - 1, nx)]);
            rhs[p(i, j, nx)] += alpha * (1 - 2. * alpha) * (dx * u0_xixi[p(i + 1, j, nx)] + dy * u0_etaeta[p(i + 1, j, nx)]);
            rhs[p(i, j, nx)] += alpha * alpha * (dx * u0_xixi[p(i + 1, j + 1, nx)] + dy * u0_etaeta[p(i + 1, j + 1, nx)]);
            rhs[p(i, j, nx)] = std::abs(rhs[p(i, j, nx)]);
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        rhs[p(i, j, nx)] = rhs[p(i, j + 1, nx)];
        j = ny - 1;
        rhs[p(i, j, nx)] = rhs[p(i, j - 1, nx)];
    }
    for (int j = 0; j < ny - 1; ++j)  // vertical direction
    {
        int i = 0;
        rhs[p(i, j, nx)] = rhs[p(i + 1, j, nx)];
        i = nx - 1;
        rhs[p(i, j, nx)] = rhs[p(i - 1, j, nx)];
    }
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            // dummy
        }
    }        

    for (int i = 0; i < nx; ++i)
    {
        solution[i] = u0[i];
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
    solverA.compute(A);
    solution = solverA.solve(rhs);
    //solution = solverA.solveWithGuess(rhs, solution);

    for (int i = 0; i < nx; ++i)
    {
        err->push_back(solution[i]); // / (err_max + 1e-13);  // to prevent division bij zero
    }
    return err;
}


int REGULARIZATION::p(int i, int j, int nx)
{
    return j * nx + i;
}
