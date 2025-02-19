//---------------------------------------------------------------
//   programmer: J. Mooiman
//   date:       2024-07-28
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
}
REGULARIZATION::REGULARIZATION(int iter_max) :
    m_iter_max(iter_max)
{
}


void REGULARIZATION::old_version(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& err, std::vector<double>& u_giv,
    double dx, double c_psi, bool use_eq8)
{
    int nx = (int) u_giv.size();
    int iter_max = m_iter_max;
    double c_error;
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    std::vector<double> u0(nx, 0.);
    std::vector<double> u1(nx, 0.);
    std::vector<double> u0_xixi(nx, 0.);
    std::vector<double> tmp(nx, 0.);

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::SparseMatrix<double> B(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    std::ofstream log_file;
#if LOGFILE == 1
    std::filesystem::create_directory("./output");
    log_file.open("./output/janm.log", std::ios_base::app);
#endif

    for (int i = 0; i < nx; ++i)
    {
        u0[i] = u_giv[i];
    }
    for (int i = 0; i < nx; ++i)
    {
        A.coeffRef(i, i) = 1.0;
        B.coeffRef(i, i) = 1.0;
        rhs[i] = 0.0;
    }
    for (int iter = 0; iter < iter_max; ++iter)
    {
        //std::cout << "Iteration: " << iter << std::endl;
        //std::cout << "Compute second derivative" << std::endl;
        double u0xixi_max = 0.0;
        for (int i = 1; i < nx - 1; ++i)
        {
            u0_xixi[i] = (u0[i - 1] - 2. * u0[i] + u0[i + 1]);
            u0xixi_max = std::max(u0xixi_max, std::abs(u0_xixi[i]));
        }
        int i = 0;
        u0_xixi[i] = u0_xixi[i + 2];
        u0_xixi[i+1] = u0_xixi[i + 2];
        i = nx - 1;
        u0_xixi[i] = u0_xixi[i - 2];
        u0_xixi[i-1] = u0_xixi[i - 2];

        //std::cout << "Initialization matrix A and rhs" << std::endl;
        c_error = c_psi;

        for (int i = 1; i < nx - 1; ++i)
        {
            //A.coeffRef(i, i - 2) = 0.0;
            A.coeffRef(i, i - 1) = 1. / 8. - c_error;
            A.coeffRef(i, i) = 6. / 8. + 2. * c_error;
            A.coeffRef(i, i + 1) = 1. / 8. - c_error;
            //A.coeffRef(i, i + 2) = 0.0;
        }
        for (int i = 1; i < nx - 1; ++i)
        {
            tmp[i] = std::abs(u0_xixi[i]);
            rhs[i] = tmp[i];
        }
        i = 0;
        A.coeffRef(i, i) = -1.;
        A.coeffRef(i, i + 1) = 1.;
        A.coeffRef(i, i + 2) = 0.;
        rhs[i] = 0.0;

        i = nx - 1;
        A.coeffRef(i, i - 2) = 0.;
        A.coeffRef(i, i - 1) = -1.;
        A.coeffRef(i, i) = 1.;
        rhs[i] = 0.0;

#if LOGFILE == 1
        if (true)
        {
            log_file << "=== Matrix A ==========================================" << std::endl;
            log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A) << std::endl;
            log_file << "=== RHS A =============================================" << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                tmp[i] = rhs[i];
                log_file << std::setprecision(8) << std::scientific << tmp[i] << std::endl;
            }
        }
#endif
        //std::cout << "Copy previous iteration value u0 to solution vector for initial guess" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            solution[i] = u0[i];
        }

        //std::cout << "Initialize matrix for BiCGStab" << std::endl;
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
        //std::cout << "Set solver" << std::endl;
        solverA.compute(A);
        //std::cout << "Solve matrix A with BiCGStab" << std::endl;
        solution = solverA.solve(rhs);
        //solution = solverA.solveWithGuess(rhs, solution);

        for (int i = 0; i < nx; ++i)
        {
            err[i] = solution[i];  // to prevent division bij zero
        }
#if LOGFILE == 1
        log_file << "--- eq8 -----------------------------------------------" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << err[i] << std::endl;
        }
#endif

        if (!use_eq8)
        {
            for (int i = 0; i < nx; ++i)
            {
                psi[i] = c_psi * dx * dx;
            }
        }
        else
        {
            for (int i = 0; i < nx; ++i)
            {
                psi[i] = c_psi * dx * dx * err[i];
                //psi[i] = c_psi * (u_const * dx) * err[i];
            }
        }

        //std::cout << "Defining matrix B and rhs" << std::endl;

        for (int i = 1; i < nx - 2; ++i)
        {
            double psi_im12 = 0.5 * (psi[i] + psi[i - 1]);
            double psi_ip12 = 0.5 * (psi[i + 1] + psi[i]);
            //psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
            //psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

            //B.coeffRef(i, i - 2) = 0.0;
            B.coeffRef(i, i - 1) = dx * 1. / 8. - psi_im12 / dx;
            B.coeffRef(i, i) = dx * 6./8. + psi_ip12 / dx + psi_im12 / dx;
            B.coeffRef(i, i + 1) = dx * 1. / 8. - psi_ip12 / dx;
            //B.coeffRef(i, i + 2) = 0.0;
        }
        for (int i = 1; i < nx - 1; ++i)
        {
            tmp[i] = dx * (1. / 8. * u_giv[i - 1] + 6. / 8. * u_giv[i] + 1. / 8. * u_giv[i + 1]);
            rhs[i] = tmp[i];
        }
        //std::cout << std::endl;
        i = 0;
        B.coeffRef(i, i) = 1.0;
        B.coeffRef(i, i + 1) = 0.0;
        B.coeffRef(i, i + 2) = 0.0;
        tmp[i] = u_giv[i];
        rhs[i] = tmp[i];
        B.coeffRef(i + 1, i) = 0.0;
        B.coeffRef(i + 1, i + 1) = 1.0;
        B.coeffRef(i + 1, i + 2) = 0.0;
        tmp[i+1] = u_giv[i+1];
        rhs[i+1] = tmp[i+1];

        i = nx - 1;
        B.coeffRef(i, i - 2) = 0.0;
        B.coeffRef(i, i - 1) = 0.0;
        B.coeffRef(i, i) = 1.0;
        tmp[i] = u_giv[i];
        rhs[i] = tmp[i];
        B.coeffRef(i - 1, i - 2) = 0.0;
        B.coeffRef(i - 1, i - 1) = 1.0;
        B.coeffRef(i - 1, i    ) = 0.0;
        tmp[i - 1] = u_giv[i - 1];
        rhs[i - 1] = tmp[i - 1];

#if LOGFILE == 1
        if (true)
        {
            //log_file << "=== Matrix B ==========================================" << std::endl;
            //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(B) << std::endl;
            //log_file << "=== RHS B =============================================" << std::endl;
            //log_file << std::setprecision(8) << std::scientific << rhs << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                tmp[i] = rhs[i];
            }
            log_file << "--- eq7 -----------------------------------------------" << std::endl;
            for (int i = 0; i < nx; ++i)
            {
                log_file << std::setprecision(8) << std::scientific << tmp[i] << std::endl;
            }
        }
#endif
        //std::cout << "Copy previous iteration value u0 to solution vector for initial guess" << std::endl;
        for (int i = 0; i < nx; ++i)
        {
            solution[i] = u0[i];
        }

        //std::cout << "Initialize matrix for BiCGStab" << std::endl;
        Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverB;
        //std::cout << "Set solver" << std::endl;
        solverB.compute(B);
        //std::cout << "Solve matrix B with BiCGStab" << std::endl;
        //solution = solverB.solve(rhs);
        solution = solverB.solveWithGuess(rhs, solution);
        //std::cout << "Copy solution to u1-vector" << std::endl;

        diff_max1 = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(solution[i] - u_giv[i]));
            u0[i] = solution[i];
        }
        if (std::abs(diff_max1 - diff_max0) > 1e-12)
        {
            //log_file << std::setprecision(8) << std::scientific << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            diff_max0 = diff_max1;
        }
        else
        {
            //log_file << std::setprecision(8) << std::scientific  << " --- Convergence: " << std::abs(diff_max1 - diff_max0) << std::endl;
            break;
        }
        //log_file << std::endl;
    }
#if LOGFILE == 1
    log_file.close();
#endif    
    for (int i = 0; i < nx; ++i)
    {
        u_out[i] = u0[i];
    }
    return;
}

void REGULARIZATION::given_function(std::vector<double>& u_out, std::vector<double>& psi, std::vector<double>& eq8, std::vector<double>& u_giv_in,
    double dx, double c_psi, bool use_eq8)
{
    int nx = (int) u_giv_in.size();
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
        u0_xixi[i] = u0_xixi[i + 2];
        u0_xixi[i+1] = u0_xixi[i + 2];
        i = nx - 1;
        u0_xixi[i] = u0_xixi[i - 2];
        u0_xixi[i-1] = u0_xixi[i - 2];

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

std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(double c_error, std::vector<double> u0, std::vector<double> u0_xixi)
{
    int nx = u0_xixi.size();
    auto err = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nx, 0.0);

    Eigen::SparseMatrix<double> A(nx, nx);
    Eigen::VectorXd solution(nx);           // solution vector u
    Eigen::VectorXd rhs(nx);                // RHS

    for (int i = 1; i < nx - 1; ++i)
    {
        //A.coeffRef(i, i - 2) = 0.0;
        A.coeffRef(i, i - 1) = 1. / 8. - c_error;
        A.coeffRef(i, i) = 6. / 8. + 2. * c_error;
        A.coeffRef(i, i + 1) = 1. / 8. - c_error;
        //A.coeffRef(i, i + 2) = 0.0;
    }
    for (int i = 1; i < nx - 1; ++i)
    {
        tmp[i] = std::abs(u0_xixi[i]);
        rhs[i] = tmp[i];
    }
    int i = 0;
    A.coeffRef(i, i) = 1.;
    A.coeffRef(i, i + 1) = 0.;
    A.coeffRef(i, i + 2) = 0.;
    rhs[i] = rhs[i + 1];
    //i = 1;
    //A.coeffRef(i, i - 1) = 1. / 8. - c_error;
    //A.coeffRef(i, i) = 6. / 8. + 2. * c_error;
    //A.coeffRef(i, i + 1) = 1. / 8. - c_error;
    //A.coeffRef(i, i + 2) = 0.0;
    //rhs[i] = rhs[i];

    i = nx - 1;
    A.coeffRef(i, i - 2) = 0.;
    A.coeffRef(i, i - 1) = 0.;
    A.coeffRef(i, i) = 1.;
    rhs[i] = rhs[i - 1];
    //i = nx - 2;
    //A.coeffRef(i, i - 2) = 0.0;
    //A.coeffRef(i, i - 1) = 1. / 8. - c_error;
    //A.coeffRef(i, i) = 6. / 8. + 2. * c_error;
    //A.coeffRef(i, i + 1) = 1. / 8. - c_error;
    //rhs[i] = rhs[i];

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


