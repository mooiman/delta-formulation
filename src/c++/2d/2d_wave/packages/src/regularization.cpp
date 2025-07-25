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
    m_g = 10.0;

    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    eps_smooth = 1e-12;
    m_u0_is_smooth = 1.e-10;
}
REGULARIZATION::REGULARIZATION(int iter_max, double g) :
    m_iter_max(iter_max),
    m_g(g)
{
    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    eps_smooth = 1e-12;
    m_u0_is_smooth = 1.e-10;
}

void REGULARIZATION::given_function( 
    std::vector<double>& u_out, std::vector<double>& psi, 
    std::vector<double>& u_giv_in, 
    int nx, int ny, double dx, double dy, double c_psi, std::ofstream& log_file)
{
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    int nxny = nx * ny;
    std::vector<double> u_giv(nxny, 0.);
    std::vector<double> u0(nxny, 0.);
    std::vector<double> u1(nxny, 0.);
    std::vector<double> u0_xixi(nxny, 0.);
    std::vector<double> u0_etaeta(nxny, 0.);
    std::vector<double> eq8(nxny, 0.);
    std::vector<double> tmp(nxny, 0.);

    const auto [u_giv_in_min, u_giv_in_max] = std::minmax_element(u_giv_in.begin(), u_giv_in.end());
    double min_range = 0.0001;
    double u_giv_range = *u_giv_in_max - *u_giv_in_min > min_range ? *u_giv_in_max - *u_giv_in_min : min_range;
    u_giv_range = 1.0;
    
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
        double u0_xixi_max = 0.0;
        double u0_etaeta_max = 0.0;
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                int p_0 = p_index(i, j    , ny);
                int p_n = p_index(i, j + 1, ny);
                int p_e = p_index(i + 1, j, ny);
                int p_s = p_index(i, j - 1, ny);
                int p_w = p_index(i - 1, j, ny);
                u0_xixi[p_0]   = std::abs((u0[p_e] - 2. * u0[p_0] + u0[p_w]));
                u0_etaeta[p_0] = std::abs((u0[p_n] - 2. * u0[p_0] + u0[p_s]));
                u0_xixi_max = std::max(u0_xixi_max, std::abs(u0_xixi[p_0]));
                u0_etaeta_max = std::max(u0_etaeta_max, std::abs(u0_etaeta[p_0]));
            }
        }
        double smooth= std::max(u0_xixi_max, u0_etaeta_max);
        if (smooth < m_u0_is_smooth)
        {
            for (int i = 0; i < nxny; ++i)
            {
                //u_out[i] = u_giv_in[i];
            }
            //return;
        }

        for (int i = 0; i < nx; ++i)  // horizontal direction
        {
            int j = 0;
            int p_0  = p_index(i, j    , ny);
            int p_n  = p_index(i, j + 1, ny);
            int p_nn = p_index(i, j + 2, ny);
            u0_etaeta[p_0] = u0_etaeta[p_nn];
            u0_etaeta[p_n] = u0_etaeta[p_nn];

            j = ny - 1;
            p_0      = p_index(i, j    , ny);
            int p_s  = p_index(i, j - 1, ny);
            int p_ss = p_index(i, j - 2, ny);
            u0_etaeta[p_0] = u0_etaeta[p_ss];
            u0_etaeta[p_s] = u0_etaeta[p_ss];
        }
        for (int j = 0; j < ny; ++j)  // vertical direction
        {
            int i = 0;
            int p_0  = p_index(i    , j, ny);
            int p_e  = p_index(i + 1, j, ny);
            int p_ee = p_index(i + 2, j, ny);
            u0_xixi[p_0] = u0_xixi[p_ee];
            u0_xixi[p_e] = u0_xixi[p_ee];

            i = nx - 1;
            p_0      = p_index(i    , j, ny);
            int p_w  = p_index(i - 1, j, ny);
            int p_ww = p_index(i - 2, j, ny);
            u0_xixi[p_0] = u0_xixi[p_ww];
            u0_xixi[p_w] = u0_xixi[p_ww];
        }

//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(nx, ny, dx, dy, c_psi, u0, u0_xixi, u0_etaeta, log_file));
//------------------------------------------------------------------------------
        double wortel = std::sqrt(dx + dy);
        for (int i = 0; i < nxny; ++i)
        {
            psi[i] = c_psi * (dx * dx + dy * dy) * eq8[i]/wortel;
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(nx, ny, dx, dy, psi, u_giv, log_file));
//------------------------------------------------------------------------------

        diff_max1 = 0.0;
        for (int i = 0; i < nxny; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(u0[i] - u_giv[i]));
        }
        if (std::abs(diff_max1 - diff_max0) > eps_smooth)
        {
            diff_max0 = diff_max1;
        }
        else
        {
            m_u0_is_smooth = std::max(u0_xixi_max, u0_etaeta_max);
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
}
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(int nx, int ny, double dx, double dy, 
    std::vector<double> psi, std::vector<double> u_giv, std::ofstream& log_file)
{
    int nxny = nx * ny;
    auto u = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nxny, 0.0);

    Eigen::SparseMatrix<double> B(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector u
    Eigen::VectorXd rhs(nxny);                // RHS

    for (int i = 2; i < nx - 2; ++i)
    {
        for (int j = 2; j < ny - 2; ++j)
        {
            int p_0  = p_index(i    , j    , ny);
            int p_n  = p_index(i    , j + 1, ny);
            int p_ne = p_index(i + 1, j + 1, ny);
            int p_e  = p_index(i + 1, j    , ny);
            int p_se = p_index(i + 1, j - 1, ny);
            int p_s  = p_index(i    , j - 1, ny);
            int p_sw = p_index(i - 1, j - 1, ny);
            int p_w  = p_index(i - 1, j    , ny);
            int p_nw = p_index(i - 1, j + 1, ny);
            
            B.coeffRef(p_0, p_sw) = dx * dy * m_mass[0] * m_mass[0];
            B.coeffRef(p_0, p_s ) = dx * dy * m_mass[0] * m_mass[1];
            B.coeffRef(p_0, p_se) = dx * dy * m_mass[0] * m_mass[2];
            
            B.coeffRef(p_0, p_w) = dx * dy * m_mass[1] * m_mass[0];
            B.coeffRef(p_0, p_0) = dx * dy * m_mass[1] * m_mass[1];
            B.coeffRef(p_0, p_e) = dx * dy * m_mass[1] * m_mass[2];
            
            B.coeffRef(p_0, p_nw) = dx * dy * m_mass[2] * m_mass[0];
            B.coeffRef(p_0, p_n ) = dx * dy * m_mass[2] * m_mass[1];
            B.coeffRef(p_0, p_ne) = dx * dy * m_mass[2] * m_mass[2];

            // psi should be computed on the 8 interfaces of the control volume
            double psi_n = 0.5 * (psi[p_n] + psi[p_0]);
            double psi_e = 0.5 * (psi[p_e] + psi[p_0]);
            double psi_s = 0.5 * (psi[p_s] + psi[p_0]);
            double psi_w = 0.5 * (psi[p_w] + psi[p_0]);
            //psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
            //psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

            B.coeffRef(p_0, p_s) += -psi_s * dx / dy;
            B.coeffRef(p_0, p_w) += -psi_w * dy / dx;
            B.coeffRef(p_0, p_0) +=  psi_s * dx / dy + psi_w * dy / dx + psi_n * dy / dx + psi_e * dx / dy;
            B.coeffRef(p_0, p_e) += -psi_e * dy / dx;
            B.coeffRef(p_0, p_n) += -psi_n * dx / dy;
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        int p_0  = p_index(i, j, ny);
        int p_n  = p_index(i, j + 1, ny);
        int p_nn = p_index(i, j + 2, ny);
        B.coeffRef(p_0 , p_0 ) =  1.0;
        B.coeffRef(p_0 , p_n ) = -2.0;
        B.coeffRef(p_0 , p_nn) =  1.0;

        j = 1;
        int p_s = p_index(i, j - 1, ny);
        p_0     = p_index(i, j    , ny);
        p_n     = p_index(i, j + 1, ny);
        B.coeffRef(p_0 , p_s) = 0.0;
        B.coeffRef(p_0 , p_0) = 1.0;
        B.coeffRef(p_0 , p_n) = 0.0;

        j = ny - 1;
        p_0      = p_index(i, j, ny);
        p_s      = p_index(i, j - 1, ny);
        int p_ss = p_index(i, j - 2, ny);
        B.coeffRef(p_0, p_0 ) =  1.0;
        B.coeffRef(p_0, p_s ) = -2.0;
        B.coeffRef(p_0, p_ss) =  1.0;

        j = ny - 2;
        p_s = p_index(i, j - 1, ny);
        p_0 = p_index(i, j    , ny);
        p_n = p_index(i, j + 1, ny);
        B.coeffRef(p_0, p_0 ) = 1.0;
        B.coeffRef(p_0, p_s ) = 0.0;
        B.coeffRef(p_0, p_ss) = 0.0;
    }
    for (int j = 1; j < ny - 1; ++j)  // vertical direction
    {
        int i = 0;
        int p_0  = p_index(i, j, ny);
        int p_e  = p_index(i + 1, j, ny);
        int p_ee = p_index(i + 2, j, ny);
        B.coeffRef(p_0 , p_0 ) = 1.0;
        B.coeffRef(p_0 , p_e ) = -2.0;
        B.coeffRef(p_0 , p_ee) = 1.0;

        i = 1;
        int p_w = p_index(i - 1, j, ny);
        p_0 = p_index(i    , j, ny);
        p_e = p_index(i + 1, j, ny);
        B.coeffRef(p_0 , p_w) = 0.0;
        B.coeffRef(p_0 , p_0) = 1.0;
        B.coeffRef(p_0 , p_e) = 0.0;


        i = nx - 1;
        p_0      = p_index(i    , j, ny);
        p_w      = p_index(i - 1, j, ny);
        int p_ww = p_index(i - 2, j, ny);
        B.coeffRef(p_0 , p_0 ) = 1.0;
        B.coeffRef(p_0 , p_w ) = -2.0;
        B.coeffRef(p_0 , p_ww) = 1.0;

        i = nx - 2;
        p_w = p_index(i - 1, j, ny);
        p_0 = p_index(i    , j, ny);
        p_e = p_index(i + 1, j, ny);
        B.coeffRef(p_0 , p_w) = 0.0;
        B.coeffRef(p_0 , p_0) = 1.0;
        B.coeffRef(p_0 , p_e) = 0.0;
    }
    // corner are set after setting the right hand side
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            int p_0  = p_index(i    , j    , ny);
            int p_n  = p_index(i    , j + 1, ny);
            int p_ne = p_index(i + 1, j + 1, ny);
            int p_e  = p_index(i + 1, j    , ny);
            int p_se = p_index(i + 1, j - 1, ny);
            int p_s  = p_index(i    , j - 1, ny);
            int p_sw = p_index(i - 1, j - 1, ny);
            int p_w  = p_index(i - 1, j    , ny);
            int p_nw = p_index(i - 1, j + 1, ny);
            
            rhs[p_0]  = dx * dy * m_mass[0] * m_mass[0] * u_giv[p_sw];
            rhs[p_0] += dx * dy * m_mass[0] * m_mass[1] * u_giv[p_s ];
            rhs[p_0] += dx * dy * m_mass[0] * m_mass[2] * u_giv[p_se];

            rhs[p_0] += dx * dy * m_mass[1] * m_mass[0] * u_giv[p_w];
            rhs[p_0] += dx * dy * m_mass[1] * m_mass[1] * u_giv[p_0];
            rhs[p_0] += dx * dy * m_mass[1] * m_mass[2] * u_giv[p_e];

            rhs[p_0] += dx * dy * m_mass[2] * m_mass[0] * u_giv[p_nw];
            rhs[p_0] += dx * dy * m_mass[2] * m_mass[1] * u_giv[p_n ];
            rhs[p_0] += dx * dy * m_mass[2] * m_mass[2] * u_giv[p_ne];
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        int p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        j = 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = u_giv[p_0];
        
        j = ny - 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        j = ny - 2;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = u_giv[p_0];
    }
    for (int j = 1; j < ny - 1; ++j)  // vertical direction
    {
        int i = 0;
        int p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        i = 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = u_giv[p_0];
        
        i = nx - 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        i = nx - 2;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = u_giv[p_0];
    }
    // corner points
    int p_0 = p_index(0, 0, ny);
    B.coeffRef(p_0, p_0) = 1.0;
    rhs[p_0] = u_giv[p_0];
    p_0 = p_index(nx - 1, 0, ny);
    B.coeffRef(p_0, p_0) = 1.0;
    rhs[p_0] = u_giv[p_0];
    p_0 = p_index(nx - 1, ny - 1, ny);
    B.coeffRef(p_0, p_0) = 1.0;
    rhs[p_0] = u_giv[p_0];
    p_0 = p_index(0, ny - 1, ny);
    B.coeffRef(p_0, p_0) = 1.0;
    rhs[p_0] = u_giv[p_0];

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverB;
    solverB.compute(B);
    solution = solverB.solve(rhs);
    //solution = solverB.solveWithGuess(rhs, solution);

    for (int i = 0; i < nxny; ++i)
    {
        u->push_back(solution[i]);
    }

    return u;
};


std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(int nx, int ny, double dx, double dy, double c_error, 
std::vector<double> u0, std::vector<double> u0_xixi, std::vector<double> u0_etaeta, std::ofstream& log_file)
{
    int nxny = nx*ny;
    auto err = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nxny, NAN);

    Eigen::SparseMatrix<double> A(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector u
    Eigen::VectorXd rhs(nxny);                // RHS

    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            int p_0  = p_index(i    , j    , ny);
            int p_n  = p_index(i    , j + 1, ny);
            int p_ne = p_index(i + 1, j + 1, ny);
            int p_e  = p_index(i + 1, j    , ny);
            int p_se = p_index(i + 1, j - 1, ny);
            int p_s  = p_index(i    , j - 1, ny);
            int p_sw = p_index(i - 1, j - 1, ny);
            int p_w  = p_index(i - 1, j    , ny);
            int p_nw = p_index(i - 1, j + 1, ny);

            A.coeffRef(p_0, p_sw) = m_mass[0] * m_mass[0];
            A.coeffRef(p_0, p_s ) = m_mass[0] * m_mass[1] - c_error;
            A.coeffRef(p_0, p_se) = m_mass[0] * m_mass[2];

            A.coeffRef(p_0, p_w) = m_mass[1] * m_mass[0] - c_error;
            A.coeffRef(p_0, p_0) = m_mass[1] * m_mass[1] + 4. * c_error;
            A.coeffRef(p_0, p_e) = m_mass[1] * m_mass[2] - c_error;

            A.coeffRef(p_0, p_nw) = m_mass[2] * m_mass[0];
            A.coeffRef(p_0, p_n ) = m_mass[2] * m_mass[1] - c_error;
            A.coeffRef(p_0, p_ne) = m_mass[2] * m_mass[2];
        }
    }
    for (int i = 1; i < nx - 1; ++i)  // horizontal direction
    {
        int j = 0;
        int p_0  = p_index(i, j    , ny);
        int p_n  = p_index(i, j + 1, ny);
        int p_nn = p_index(i, j + 2, ny);
        A.coeffRef(p_0, p_0 ) =  1.0;
        A.coeffRef(p_0, p_n ) = -2.0;
        A.coeffRef(p_0, p_nn) =  1.0;

        j = 1;
        int p_s  = p_index(i, j - 1, ny);
        p_0      = p_index(i, j    , ny);
        p_n      = p_index(i, j + 1, ny);
        A.coeffRef(p_0, p_s) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_n) = 0.0;

        j = ny - 1;
        p_0      = p_index(i, j    , ny);
        p_s      = p_index(i, j - 1, ny);
        int p_ss = p_index(i, j - 2, ny);
        A.coeffRef(p_0, p_0 ) =  1.0;
        A.coeffRef(p_0, p_s ) = -2.0;
        A.coeffRef(p_0, p_ss) =  1.0;

        j = ny - 2;
        p_s  = p_index(i, j - 1, ny);
        p_0  = p_index(i, j    , ny);
        p_n  = p_index(i, j + 1, ny);
        A.coeffRef(p_0, p_s) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_n) = 0.0;
    }
    for (int j = 0; j < ny; ++j)  // vertical direction
    {
        int i = 0;
        int p_0  = p_index(i, j    , ny);
        int p_e  = p_index(i + 1, j, ny);
        int p_ee = p_index(i + 2, j, ny);
        A.coeffRef(p_0, p_0 ) =  1.0;
        A.coeffRef(p_0, p_e ) = -2.0;
        A.coeffRef(p_0, p_ee) =  1.0;

        i = 2;
        int p_w = p_index(i - 1, j, ny);
        p_0     = p_index(i    , j, ny);
        p_e     = p_index(i + 1, j, ny);
        A.coeffRef(p_0, p_w) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_e) = 0.0;

        i = nx - 1;
        p_0      = p_index(i    , j, ny);
        p_w      = p_index(i - 1, j, ny);
        int p_ww = p_index(i - 2, j, ny);
        A.coeffRef(p_0, p_0 ) =  1.0;
        A.coeffRef(p_0, p_w)  = -2.0;
        A.coeffRef(p_0, p_ww) =  1.0;

        i = nx - 2;
        p_w = p_index(i - 1, j, ny);
        p_0 = p_index(i    , j, ny);
        p_e = p_index(i + 1, j, ny);
        A.coeffRef(p_0, p_w) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_e) = 0.0;
    }
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            int p_0  = p_index(i    , j    , ny);
            int p_n  = p_index(i    , j + 1, ny);
            int p_ne = p_index(i + 1, j + 1, ny);
            int p_e  = p_index(i + 1, j    , ny);
            int p_se = p_index(i + 1, j - 1, ny);
            int p_s  = p_index(i    , j - 1, ny);
            int p_sw = p_index(i - 1, j - 1, ny);
            int p_w  = p_index(i - 1, j    , ny);
            int p_nw = p_index(i - 1, j + 1, ny);
            rhs[p_0]  = (u0_xixi[p_sw] + u0_etaeta[p_sw]);
            rhs[p_0] += (u0_xixi[p_s ] + u0_etaeta[p_s]);
            rhs[p_0] += (u0_xixi[p_se] + u0_etaeta[p_se]);
            
            rhs[p_0] += (u0_xixi[p_w] + u0_etaeta[p_w]);
            rhs[p_0] += (u0_xixi[p_0] + u0_etaeta[p_0]);
            rhs[p_0] += (u0_xixi[p_e] + u0_etaeta[p_e]);
            
            rhs[p_0] += (u0_xixi[p_nw] + u0_etaeta[p_nw]);
            rhs[p_0] += (u0_xixi[p_n ] + u0_etaeta[p_n ]);
            rhs[p_0] += (u0_xixi[p_ne] + u0_etaeta[p_ne]);
            rhs[p_0] = std::abs(rhs[p_0]);
        }
    }
    for (int i = 0; i < nx; ++i)  // horizontal direction
    {
        int j = 0;
        int p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        j = 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = rhs[p_0];
        
        j = ny - 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = 0;

        j = ny - 2;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = rhs[p_0];
    }
    for (int j = 0; j < ny; ++j)  // vertical direction
    {
        int i = 0;
        int p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        i = 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = rhs[p_0];

        i = nx - 1;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = 0.0;

        i = nx - 2;
        p_0 = p_index(i, j, ny);
        rhs[p_0] = rhs[p_0];
    }
    // TODO: corner points

    for (int i = 0; i < nxny; ++i)
    {
        solution[i] = u0[i];
    }

    if (false)
    {
        log_file << "=== Matrix eq8 ========================================" << std::endl;
        log_file << std::setprecision(4) << std::scientific << Eigen::MatrixXd(A) << std::endl;
        //log_file << "=== Main diagonal eq8 =================================" << std::endl;
        //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).diagonal() << std::endl;
        //log_file << "=== Eigen values eq8 ==================================" << std::endl;
        //log_file << std::setprecision(8) << std::scientific << Eigen::MatrixXd(A).eigenvalues() << std::endl;
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
    solverA.compute(A);
    solverA.setTolerance(1e-08);
    solution = solverA.solve(rhs);
    //solution = solverA.solveWithGuess(rhs, solution);

    if (false)
    {
        log_file << "=== RHS, solution eq8 =================================" << std::endl;
        for (int i = 0; i < nxny; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << rhs[i] << "  " << solution[i] << std::endl;
        }
        log_file << "=======================================================" << std::endl;
    }

    for (int i = 0; i < nxny; ++i)
    {
        err->push_back(solution[i]); // / (err_max + 1e-13);  // to prevent division bij zero
    }
    return err;
}
int REGULARIZATION::p_index(int i, int j, int ny)
{
    return i * ny + j;
}
