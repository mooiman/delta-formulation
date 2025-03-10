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
                int p_0 = p_index(i, j, nx);
                int p_n = p_index(i, j + 1, nx);
                int p_e = p_index(i + 1, j, nx);
                int p_s = p_index(i, j - 1, nx);
                int p_w = p_index(i - 1, j, nx);
                u0_xixi[p_0]   = std::abs((u0[p_e] - 2. * u0[p_0] + u0[p_w]));
                u0_etaeta[p_0] = std::abs((u0[p_n] - 2. * u0[p_0] + u0[p_s]));
                u0xixi_max = std::max(u0xixi_max, std::abs(u0_xixi[i]));
                u0etaeta_max = std::max(u0etaeta_max, std::abs(u0_etaeta[i]));
            }
        }
        for (int i = 0; i < nx; ++i)  // horizontal direction
        {
            int j = 0;
            int p_0 = p_index(i, j, nx);
            int p_n = p_index(i, j + 1, nx);
            int p_nn = p_index(i, j + 2, nx);
            int p_s = p_index(i, j - 1, nx);
            int p_ss = p_index(i, j - 2, nx);
            u0_etaeta[p_0] = u0_etaeta[p_nn];
            u0_etaeta[p_n] = u0_etaeta[p_nn];
            j = nx - 1;
            u0_etaeta[p_0] = u0_etaeta[p_ss];
            u0_etaeta[p_s] = u0_etaeta[p_ss];
        }
        for (int j = 0; j < ny; ++j)  // vertical direction
        {
            int i = 0;
            int p_0 = p_index(i, j, nx);
            int p_e = p_index(i + 1, j, nx);
            int p_ee = p_index(i + 2, j, nx);
            int p_w = p_index(i - 1, j, nx);
            int p_ww = p_index(i - 2, j, nx);
            u0_xixi[p_0] = u0_xixi[p_ee];
            u0_xixi[p_e] = u0_xixi[p_ee];
            i = nx - 1;
            u0_xixi[p_0] = u0_xixi[p_ww];
            u0_xixi[p_w] = u0_xixi[p_ww];
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
                    int p_0 = p_index(i, j, nx)
                    psi[p_0] = c_psi * dx * dx * eq8[p_0];
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

    double alpha = 1. / 8.;                               // Linear (spatial) interpolation coefficient
    std::vector<double> mass(3, 0);
    mass[0] = alpha;
    mass[1] = 1.0 - 2. * alpha;
    mass[2] = alpha;

    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            A.coeffRef(i-1, j - 1) = mass[0] * mass[0];
            A.coeffRef(i  , j - 1) = mass[0] * mass[1] - c_error;
            A.coeffRef(i+1, j - 1) = mass[0] * mass[2];
            A.coeffRef(i-1, j    ) = mass[1] * mass[0] - c_error;
            A.coeffRef(i  , j    ) = mass[1] * mass[1] + 4. * c_error;
            A.coeffRef(i+1, j    ) = mass[1] * mass[2] - c_error;
            A.coeffRef(i-1, j + 1) = mass[2] * mass[0];
            A.coeffRef(i  , j + 1) = mass[2] * mass[1] - c_error;
            A.coeffRef(i+1, j + 1) = mass[2] * mass[2];
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
            int p_0  = p_index(i    , j    , nx);
            int p_n  = p_index(i    , j + 1, nx);
            int p_ne = p_index(i + 1, j + 1, nx);
            int p_e  = p_index(i + 1, j    , nx);
            int p_se = p_index(i + 1, j - 1, nx);
            int p_s  = p_index(i    , j - 1, nx);
            int p_sw = p_index(i - 1, j - 1, nx);
            int p_w  = p_index(i - 1, j    , nx);
            int p_nw = p_index(i - 1, j + 1, nx);
            rhs[p_0]  = mass[0] * mass[0] * (dx * u0_xixi[p_sw] + dy * u0_etaeta[p_sw]);
            rhs[p_0] += mass[0] * mass[1] * (dx * u0_xixi[p_s] + dy * u0_etaeta[p_s]);
            rhs[p_0] += mass[0] * mass[2] * (dx * u0_xixi[p_se] + dy * u0_etaeta[p_se]);
            rhs[p_0] += mass[1] * mass[0] * (dx * u0_xixi[p_w] + dy * u0_etaeta[p_w]);
            rhs[p_0] += mass[1] * mass[1] * (dx * u0_xixi[p_0] + dy * u0_etaeta[p_0]);
            rhs[p_0] += mass[1] * mass[2] * (dx * u0_xixi[p_e] + dy * u0_etaeta[p_e]);
            rhs[p_0] += mass[2] * mass[0] * (dx * u0_xixi[p_nw] + dy * u0_etaeta[p_nw]);
            rhs[p_0] += mass[2] * mass[1] * (dx * u0_xixi[p_n] + dy * u0_etaeta[p_n]);
            rhs[p_0] += mass[2] * mass[2] * (dx * u0_xixi[p_ne] + dy * u0_etaeta[p_ne]);
            rhs[p_0] = std::abs(rhs[p_0]);
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        int p_0 = p_index(i, j, nx);
        int p_n = p_index(i, j + 1, nx);
        rhs[p_0] = rhs[p_n];
        j = ny - 1;
        p_0 = p_index(i, j, nx);
        int p_s = p_index(i, j - 1, nx);
        rhs[p_0] = rhs[p_s];
    }
    for (int j = 0; j < ny - 1; ++j)  // vertical direction
    {
        int i = 0;
        int p_0 = p_index(i, j, nx);
        int p_e = p_index(i + 1, j, nx);
        rhs[p_0] = rhs[p_e];
        i = nx - 1;
        p_0 = p_index(i, j, nx);
        int p_w = p_index(i - 1, j, nx);
        rhs[p_0] = rhs[p_w];
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


int REGULARIZATION::p_index(int i, int j, int nx)
{
    return j * nx + i;
}
