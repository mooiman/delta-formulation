//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2025-03-15
//---------------------------------------------------------------
//   DESCRIPTION
//
//   Regularization of given function
//

#include "regularization.h"
#define LOGFILE 0

REGULARIZATION::REGULARIZATION()
{
    m_iter_max = 100;
    m_g = 10.0;
    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    m_u0_xixi_smooth = 0.0;
}
REGULARIZATION::REGULARIZATION(int iter_max, double g) :
    m_iter_max(iter_max),
    m_g(g)
{
    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    m_u0_xixi_smooth = 0.0;
}

void REGULARIZATION::given_function(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& eq8, std::vector<double>& u_giv_in,
    double dx, double c_psi, bool use_eq8)
{
    int nx = (int)u_giv_in.size();
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    std::vector<double> u_giv(nx, 0.);
    std::vector<double> u0(nx, 0.);
    std::vector<double> u1(nx, 0.);
    std::vector<double> u0_xixi(nx, 0.);
    std::vector<double> tmp(nx, 0.);

    std::ofstream log_file;
#if LOGFILE == 1
    std::filesystem::create_directory("./output");
    log_file.open("./output/janm.log", std::ios_base::app);
#endif
    const auto [u_giv_in_min, u_giv_in_max] = std::minmax_element(u_giv_in.begin(), u_giv_in.end());
    double min_range = 0.0001;
    double u_giv_range = *u_giv_in_max - *u_giv_in_min > min_range ? *u_giv_in_max - *u_giv_in_min : min_range;
    u_giv_range = 1.0;

    for (int i = 0; i < nx; ++i)
    {
        // scaling of the input array u_giv_in
        u0[i] = u_giv_in[i] / u_giv_range;
        u_giv[i] = u0[i];        
    }
    for (int iter = 0; iter < m_iter_max; ++iter)
    {
        //std::cout << "Iteration: " << iter << std::endl;
        //std::cout << "Compute second derivative" << std::endl;
        double u0_xixi_max = 0.0;
        for (int i = 1; i < nx - 1; ++i)
        {
            u0_xixi[i] = (u0[i - 1] - 2. * u0[i] + u0[i + 1]);
            u0_xixi_max = std::max(u0_xixi_max, std::abs(u0_xixi[i]));
        }
        int i = 0;
        u0_xixi[i] = u0_xixi[i + 1];
        //u0_xixi[i+1] = u0_xixi[i + 2];
        i = nx - 1;
        u0_xixi[i] = u0_xixi[i -1];
        //u0_xixi[i-1] = u0_xixi[i - 2];
        if (u0_xixi_max < 1.001 * m_u0_xixi_smooth) { return; }

//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(c_psi, u0, u0_xixi));
//------------------------------------------------------------------------------
#if LOGFILE == 1
        log_file << "--- eq8 -----------------------------------------------" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << eq8[i] << std::endl;
        }
#endif

        if (use_eq8)
        {
            for (int i = 0; i < nx; ++i)
            {
                psi[i] = c_psi * dx * dx * eq8[i];
            }
        }
        else
        {
            for (int i = 0; i < nx; ++i)
            {
                psi[i] = c_psi * dx * dx;
            }
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(dx, psi, u_giv));
//------------------------------------------------------------------------------
#if LOGFILE == 1
        log_file << "--- eq7 -----------------------------------------------" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << u0[i] << std::endl;
        }
#endif


        diff_max1 = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(u0[i] - u_giv[i]));
        }
        if (std::abs(diff_max1 - diff_max0) > 1e-12)
        {
            //log_file << std::setprecision(8) << std::scientific << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            diff_max0 = diff_max1;
        }
        else
        {
            //log_file << std::setprecision(8) << std::scientific  << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            m_u0_xixi_smooth = u0_xixi_max;
            break;
        }
        //log_file << std::endl;
    }
#if LOGFILE == 1
    log_file.close();
#endif    
    for (int i = 0; i < nx; ++i)
    {
        u_out[i] = u0[i] * u_giv_range;
    }
}

void REGULARIZATION::artificial_viscosity(std::vector<double>& psi, std::vector<double>& h, std::vector<double>& q, 
    std::vector<double>& zb, double c_psi_in, double dx)
{
    int nx = (int)h.size();
    std::vector<double> h_xixi(nx, 0.);  // second derivative of total depth in computational space
    std::vector<double> q_xixi(nx, 0.);  // second derivative of flow flux in computational space
    std::vector<double> s_xixi(nx, 0.);  // second derivative of water level: s = h + zb in computational space
    std::vector<double> Err_psi(nx, 0.);  //

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);               // solution vector 
    Eigen::VectorXd rhs(nx);                // RHS vector

    double c_psi = c_psi_in;

    for (int i = 0; i < nx; ++i)
    {
        A.coeffRef(i, i) = 1.0;
        rhs[i] = 0.0;
    }

    for (int i = 1; i < nx-1; ++i)
    {
        h_xixi[i] = (h[i - 1] - 2. * h[i] + h[i + 1]);
        q_xixi[i] = (q[i - 1] - 2. * q[i] + q[i + 1]);
        s_xixi[i] = ((h[i-1] + zb[i-1]) - 2. * (h[i] + zb[i]) + (h[i + 1] + zb[i + 1]));
    }
    int i = 0;
    h_xixi[i] = 2. * h_xixi[i + 1] - h_xixi[i + 2];
    q_xixi[i] = 2. * q_xixi[i + 1] - q_xixi[i + 2];
    s_xixi[i] = 2. * s_xixi[i + 1] - s_xixi[i + 2];
    i = nx - 1;
    h_xixi[i] = 2. * h_xixi[i - 1] - h_xixi[i - 2];
    q_xixi[i] = 2. * q_xixi[i - 1] - q_xixi[i - 2];
    s_xixi[i] = 2. * s_xixi[i - 1] - s_xixi[i - 2];
    //

    // eq. 18
    double hbar_im14;
    double hbar_ip14;
    double qbar_im14;
    double qbar_ip14;

    double c_error = 2. * c_psi; //same value as for regularization of given function
    for (int i = 1; i < nx - 1; ++i)
    {
        A.coeffRef(i, i - 1) = m_mass[0] - c_error;
        A.coeffRef(i, i    ) = m_mass[1] + 2. * c_error;
        A.coeffRef(i, i + 1) = m_mass[2] - c_error;

        hbar_im14 = 0.25 * (h[i - 1] + 3. * h[i]);
        hbar_ip14 = 0.25 * (h[i + 1] + 3. * h[i]);
        qbar_im14 = 0.25 * (q[i - 1] + 3. * q[i]);
        qbar_ip14 = 0.25 * (q[i + 1] + 3. * q[i]);

        rhs[i] = 4.0 * c_psi * dx * (
            0.5 * std::sqrt(m_g / hbar_im14) * std::abs(s_xixi[i]) +
            0.5 * std::sqrt(2.) * std::abs(q_xixi[i] / hbar_im14 - qbar_im14 * h_xixi[i] / (hbar_im14 * hbar_im14))
            + 0.5 * std::sqrt(m_g / hbar_ip14) * std::abs(s_xixi[i]) +
            0.5 * std::sqrt(2.) * std::abs(q_xixi[i] / hbar_ip14 - qbar_ip14 * h_xixi[i] / (hbar_ip14 * hbar_ip14))
            );
    }
    // eq. 19
    i = 0;
    A.coeffRef(i, i) = 1.; 
    A.coeffRef(i, i + 1) = -2.0;
    A.coeffRef(i, i + 2) = 1.;
    rhs[i] = 0.0;
    i = 1;
    A.coeffRef(i, i - 1) = 0;
    A.coeffRef(i, i) = 1.;
    A.coeffRef(i, i + 1) = 0;
    rhs[i] = rhs[i];
    i = nx - 1;
    A.coeffRef(i, i - 2) = 1.;
    A.coeffRef(i, i - 1) = -2.0;
    A.coeffRef(i, i) = 1.;
    rhs[i] = 0.0;
    i = nx - 2;
    A.coeffRef(i, i - 1) = 0.0;
    A.coeffRef(i, i) = 1.0;
    A.coeffRef(i, i + 1) = 0.0;
    rhs[i] = rhs[i];

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);
    for (int i = 0; i < nx; ++i)
    {
        psi[i] = solution[i];
    }
}

void REGULARIZATION::first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx)
{
    int nx = (int)eps.size();
    double peclet = 0.0;
    double peclet_threshold = 1.9;

    for (int i = 0; i < nx; ++i)
    {
        peclet = std::abs(u[i] * dx / eps[i]);
        psi[i] = 1.0;
        if (peclet > peclet_threshold)
        { 
            psi[i] = peclet/ peclet_threshold;
        }
    }
}

std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(double dx, std::vector<double> psi, std::vector<double> u_giv)
{
    int nx = (int) psi.size();
    auto u = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nx, 0.0);

    Eigen::SparseMatrix<double> B(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    for (int i = 1; i < nx - 2; ++i)
    {
        double psi_im12 = 0.5 * (psi[i - 1] + psi[i]);
        double psi_ip12 = 0.5 * (psi[i] + psi[i + 1]);
        //psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
        //psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

        //B.coeffRef(i, i - 2) = 0.0;
        B.coeffRef(i, i - 1) = dx * m_mass[0] - psi_im12 / dx;
        B.coeffRef(i, i) = dx * m_mass[1] + psi_ip12 / dx + psi_im12 / dx;
        B.coeffRef(i, i + 1) = dx * m_mass[2] - psi_ip12 / dx;
        //B.coeffRef(i, i + 2) = 0.0;
    }
    for (int i = 1; i < nx - 1; ++i)
    {
        tmp[i] = dx * (m_mass[0] * u_giv[i - 1] + m_mass[1] * u_giv[i] + m_mass[2] * u_giv[i + 1]);
        rhs[i] = tmp[i];
    }
    int i = 0;
    B.coeffRef(i, i) = 1.0;
    B.coeffRef(i, i + 1) = -2.0;
    B.coeffRef(i, i + 2) = 1.0;
    rhs[i] = 0.0;
    i = 1;
    B.coeffRef(i, i - 1) = 0.0;
    B.coeffRef(i, i    ) = 1.0;
    B.coeffRef(i, i + 1) = 0.0;
    tmp[i] = u_giv[i];
    rhs[i] = tmp[i];

    i = nx - 1;
    B.coeffRef(i, i - 2) = 1.0;
    B.coeffRef(i, i - 1) = -2.0;
    B.coeffRef(i, i) = 1.0;
    rhs[i] = 0.0;
    i = nx - 2;
    B.coeffRef(i, i - 1) = 0.0;
    B.coeffRef(i, i    ) = 1.0;
    B.coeffRef(i, i + 1) = 0.0;
    tmp[i] = u_giv[i];
    rhs[i] = tmp[i];

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

std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(double c_error, std::vector<double> u0, std::vector<double> u0_xixi)
{
    int nx = (int) u0_xixi.size();
    auto err = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nx, 0.0);

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    for (int i = 1; i < nx - 1; ++i)
    {
        //A.coeffRef(i, i - 2) = 0.0;
        A.coeffRef(i, i - 1) = m_mass[0] - c_error;
        A.coeffRef(i, i) = m_mass[1] + 2. * c_error;
        A.coeffRef(i, i + 1) = m_mass[2] - c_error;
        //A.coeffRef(i, i + 2) = 0.0;
    }
    for (int i = 1; i < nx - 1; ++i)
    {
        tmp[i] = std::abs(u0_xixi[i]);
        rhs[i] = tmp[i];
    }
    int i = 0;
    A.coeffRef(i, i    ) = 1.;
    A.coeffRef(i, i + 1) = -2.;
    A.coeffRef(i, i + 2) = 1.;
    rhs[i] = 0.0;
    i = 1;
    A.coeffRef(i, i - 1) = 0.0;
    A.coeffRef(i, i    ) = 1.0;
    A.coeffRef(i, i + 1) = 0.0;
    rhs[i] = rhs[i];

    i = nx - 1;
    A.coeffRef(i, i - 2) = 1.0;
    A.coeffRef(i, i - 1) = -2.0;
    A.coeffRef(i, i    ) = 1.0;
    rhs[i] = 0.0;
    i = nx - 2;
    A.coeffRef(i, i - 1) = 0.0;
    A.coeffRef(i, i    ) = 1.0;
    A.coeffRef(i, i + 1) = 0.0;
    rhs[i] = rhs[i];

    for (int i = 0; i < nx; ++i)
    {
        solution[i] = u0[i];
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
    solverA.compute(A);
    solution = solverA.solve(rhs);
    //solution = solverA.solveWithGuess(rhs, solution);

    //const auto [min, max] = std::minmax_element(solution.begin(), solution.end());
    //double err_max = std::abs(*min) < std::abs(*max) ? std::abs(*max) : std::abs(*min);

    for (int i = 0; i < nx; ++i)
    {
        err->push_back(solution[i]); // / (err_max + 1e-13);  // to prevent division bij zero
    }
    return err;
}


