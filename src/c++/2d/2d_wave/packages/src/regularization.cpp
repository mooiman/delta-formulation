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
//   u_tilde - d\dx psi d\dx u_tilde - d\dy psi d\dy u_tilde = u_giv
// 
//   Based on article:
//       From M. Borsboom
//       Applied Numerical Mathematics 26 (1998)
//       Borsboom_developerrorminadaptgridmethod_ApplNumerMath1998.pdf
//

#include "build_matrix_pattern_regularization.h"
#include "perf_timer.h"
#include "print_matrix.h"
#include "regularization.h"

REGULARIZATION::REGULARIZATION()
{
    m_iter_max = 100;
    m_g = 10.0;
    m_logging = "";

    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    m_eps_smooth = 1e-12;
    m_u0_is_smooth = 1.e-10;
}
REGULARIZATION::REGULARIZATION(int iter_max, double g, std::string logging) :
    m_iter_max(iter_max),
    m_g(g),
    m_logging(logging)
{   
    m_alpha = 1./8.;
    m_mass.push_back(m_alpha);
    m_mass.push_back(1.0 - 2. * m_alpha);
    m_mass.push_back(m_alpha);
    m_eps_smooth = 1e-12;
    m_u0_is_smooth = 1.e-10;
}
//------------------------------------------------------------------------------
void REGULARIZATION::given_function(
    std::vector<double>& u_tilde, std::vector<double>& psi_11, std::vector<double>& psi_22,
    std::vector<double>& eq8, std::vector<double>& u_giv_in,
    double c_psi, struct _grid_metric & metric, std::ofstream& log_file)
{
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;
    //std::vector<double> eq8(nxny, 0.);
    std::vector<double> u_giv(nxny, 0.);
    std::vector<double> u0(nxny, 0.);
    std::vector<double> u1(nxny, 0.);
    std::vector<double> Du0_xixi(nxny, 0.);
    std::vector<double> Du0_etaeta(nxny, 0.);
    std::vector<double> tmp(nxny, 0.);

    const auto [u_giv_in_min, u_giv_in_max] = std::minmax_element(u_giv_in.begin(), u_giv_in.end());
    double min_range = 0.0001;
    double u_giv_range = *u_giv_in_max - *u_giv_in_min > min_range ? *u_giv_in_max - *u_giv_in_min : min_range;
    u_giv_range = 1.0;

    for (size_t i = 0; i < nxny; ++i)
    {
        // scaling of the input array u_giv_in
        u0[i] = u_giv_in[i] / u_giv_range;
        u_giv[i] = u0[i];
    }
    for (int iter = 0; iter < m_iter_max; ++iter)
    {
        //std::cout << "Iteration: " << iter << std::endl;
        //std::cout << "Compute second derivative" << std::endl;
        double du0_dxi = 0.0;
        double du0_deta = 0.0;
        double Du0_xixi_max = 0.0;
        double Du0_etaeta_max = 0.0;
        for (size_t j = 1; j < ny - 1; ++j)
        {
            for (size_t i = 1; i < nx - 1; ++i)
            {
                size_t p_0 = p_index(i, j    , ny);
                size_t p_n = p_index(i, j + 1, ny);
                size_t p_e = p_index(i + 1, j, ny);
                size_t p_s = p_index(i, j - 1, ny);
                size_t p_w = p_index(i - 1, j, ny);
                du0_dxi  = 0.5 * (u0[p_e] - u0[p_0]) + (u0[p_0] - u0[p_w]);
                du0_deta = 0.5 * (u0[p_n] - u0[p_0]) + (u0[p_0] - u0[p_s]);
                Du0_xixi[p_0]   = std::abs((u0[p_e] - 2. * u0[p_0] + u0[p_w]));  // - metric.ddx_dxi2[i]  / metric.dx_dxi[i]  * du0_dxi );
                Du0_etaeta[p_0] = std::abs((u0[p_n] - 2. * u0[p_0] + u0[p_s]));  // - metric.ddy_deta2[i] / metric.dy_deta[i] * du0_deta);
                Du0_xixi_max = std::max(Du0_xixi_max, std::abs(Du0_xixi[p_0]));
                Du0_etaeta_max = std::max(Du0_etaeta_max, std::abs(Du0_etaeta[p_0]));
            }
        }
        double smooth= std::max(Du0_xixi_max, Du0_etaeta_max);
        if (smooth < m_u0_is_smooth)
        {
            for (size_t i = 0; i < nxny; ++i)
            {
                //u_out[i] = u_giv_in[i];
            }
            //return;
            smooth = smooth + 1.e-12;
        }

        for (size_t i = 0; i < nx; ++i)  // south and north boundary
        {
            size_t j = 0;
            size_t p_0  = p_index(i, j    , ny);
            size_t p_n  = p_index(i, j + 1, ny);
            size_t p_nn = p_index(i, j + 2, ny);
            Du0_etaeta[p_0] = Du0_etaeta[p_nn];
            Du0_etaeta[p_n] = Du0_etaeta[p_nn];

            j = ny - 1;
            p_0         = p_index(i, j    , ny);
            size_t p_s  = p_index(i, j - 1, ny);
            size_t p_ss = p_index(i, j - 2, ny);
            Du0_etaeta[p_0] = Du0_etaeta[p_ss];
            Du0_etaeta[p_s] = Du0_etaeta[p_ss];
        }
        for (size_t j = 0; j < ny; ++j)  // west and east boundary
        {
            size_t i = 0;
            size_t p_0  = p_index(i    , j, ny);
            size_t p_e  = p_index(i + 1, j, ny);
            size_t p_ee = p_index(i + 2, j, ny);
            Du0_xixi[p_0] = Du0_xixi[p_ee];
            Du0_xixi[p_e] = Du0_xixi[p_ee];

            i = nx - 1;
            p_0         = p_index(i    , j, ny);
            size_t p_w  = p_index(i - 1, j, ny);
            size_t p_ww = p_index(i - 2, j, ny);
            Du0_xixi[p_0] = Du0_xixi[p_ww];
            Du0_xixi[p_w] = Du0_xixi[p_ww];
        }

//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(metric, c_psi, u0, Du0_xixi, Du0_etaeta, log_file));  // eq8 is computed in computational space
//------------------------------------------------------------------------------
        for (size_t i = 0; i < nxny; ++i)
        {
            //psi_11[i] = -c_psi * metric.dx_dxi[i]  * metric.dx_dxi[i]  * eq8[i];
            //psi_22[i] = -c_psi * metric.dy_deta[i] * metric.dy_deta[i] * eq8[i];
            psi_11[i] = -c_psi * eq8[i];
            psi_22[i] = -c_psi * eq8[i];
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(metric, psi_11, psi_22, u_giv, log_file));
//------------------------------------------------------------------------------

        diff_max1 = 0.0;
        for (size_t i = 0; i < nxny; ++i)
        {
            diff_max1 = std::max(diff_max1, std::abs(u0[i] - u_giv[i]));
        }
        if (std::abs(diff_max1 - diff_max0) > m_eps_smooth)
        {
            diff_max0 = diff_max1;
        }
        else
        {
            m_u0_is_smooth = std::max(Du0_xixi_max, Du0_etaeta_max);
            break;
        }
    }
    for (size_t i = 0; i < nxny; ++i)
    {
        u_tilde[i] = u0[i] * u_giv_range;
    }
}
//------------------------------------------------------------------------------
void REGULARIZATION::artificial_viscosity(std::vector<double>& psi, 
    std::vector<double>& h, std::vector<double>& q, std::vector<double>& r, std::vector<double>& zb, 
    double c_psi_in, struct _grid_metric & metric, std::ofstream& log_file)
{
    int status = 0;
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector 
    Eigen::VectorXd rhs(nxny);                // RHS vector

    std::vector< Eigen::Triplet<double> > triplets; 
    triplets.reserve(9 * nxny); // 9-points per row, nrow = nxny; large enough to avoid reallocation
    
    status = build_matrix_pattern_regularization(triplets, nx, ny);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed(); // Very important for valuePtr() access

    if (m_logging == "pattern")
    {
        std::string header_text = "=== Matrix pattern regularization =====================";
        print_matrix_pattern(A, 1, nx, ny, header_text, log_file);
    }

    double *values = A.valuePtr();         // pointer to all non-zero values
    const int* outer = A.outerIndexPtr();   // row start pointers

    double c_psi = c_psi_in;

    START_TIMERN(Regularization diffusion);
    double psi_1 = 0.0;
    double psi_2 = 0.0;
    //
    // viscosity
    //
    // For the moment only interior nodes (2025-08-13)
    // Do not clear the rows, but add the viscosity terms

    // south-west corner
    for (size_t row = 0; row < 1; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_south_west_psi( values, row, c_eq, rhs, 
            h, 1.0, nx, ny);
    }
    // west boundary
    for (size_t row = 1; row < 1 * (ny - 1); row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        status = reg_boundary_west_psi( values, row, c_eq, rhs,
            h, psi_1, psi_2, 
            1.0, metric);
    }
    // north-west corner
    for (size_t row = 1 * (ny - 1); row < 1 * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_west_psi( values, row, c_eq, rhs, 
            h, 1.0, nx, ny);
    }

    // interior with south and north boundary
    for (size_t row = 1 * ny; row < 1 * (nx - 1) * ny; row += 1) 
    {
        size_t c_eq = outer[row    ];
        size_t p_0 = c_eq/(9);

        if (row % ny == 0) {
            // south boundary, over write coefficients
            status = reg_boundary_south_psi( values, row, c_eq, rhs,
                h, psi_1, psi_2, 
                1.0, metric);
            continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary, over write coefficients
            status = reg_boundary_north_psi( values, row, c_eq, rhs,
                h, psi_1, psi_2, 
                1.0, metric);
            continue;
        }
        status = reg_interior_matrix_psi( values, row, c_eq,
             c_psi, metric);
        // overwrite the right hand side
        status = reg_interior_rhs_psi(row, c_eq, rhs, 
            h, q, r, c_psi, m_g, metric);
    }
    // south-east corner
    for (size_t row = 1 * (nx - 1) * ny; row < 1 * (nx - 1) * ny + 1; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_south_east_psi( values, row, c_eq, rhs, 
            h, 1.0, nx, ny);
    }
    // east boundary
    for (size_t row = 1 * (nx - 1) * ny + 1; row < 1 * nx * ny - 1; row += 1) 
    {
        size_t c_eq = (size_t) outer[row    ];
        status = 0;
        size_t p_0 = c_eq/(9);

        status = reg_boundary_east_psi( values, row, c_eq, rhs,
            h, psi_1, psi_2, 
            1.0, metric);
    }
    // north-east corner
    for (size_t row = 1 * nx * ny - 1; row < 1 * nx * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_east_psi( values, row, c_eq, rhs, 
            h, 1.0, nx, ny);
    }
    STOP_TIMER(Regularization diffusion);

    if (m_logging == "matrix_psi")
    {
        std::string header_text = "=== Regularization matrix =============================";
        log_file << std::showpos << std::setprecision(3) << std::scientific;
        print_matrix(A, 1, nx, ny, header_text, log_file);
    }
    if (m_logging == "rhs_psi")
    {
        std::string header_text = "=== Regularization rhs ================================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(rhs, 1, nx, ny, header_text, log_file);
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solverA;
    solverA.compute(A);
    solverA.setTolerance(1e-12);
    solution = solverA.solve(rhs);

    if (m_logging == "matrix")
    {
        std::string header_text = "=== Regularization solution ===========================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(solution, 1, nx, ny, header_text, log_file);
    }
    for (size_t i = 0; i < nxny; ++i)
    {
        psi[i] = metric.dx_dxi[i] * metric.dy_deta[i] * std::abs(solution[i]);  // multiplied by an approximation of det(J)
    }

    return;
}
//------------------------------------------------------------------------------
void REGULARIZATION::artificial_viscosity_xi(std::vector<double>& psi, 
    std::vector<double>& h, std::vector<double>& q, std::vector<double>& zb, 
    double c_psi_in, struct _grid_metric & metric, std::ofstream& log_file)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;

    int status = 1; 

    std::vector<double> h_xixi(nxny, 0.);  // second derivative of total depth in computational space
    std::vector<double> q_xixi(nxny, 0.);  // second derivative of flow flux in computational space
    std::vector<double> s_xixi(nxny, 0.);  // second derivative of water level: s = h + zb in computational space

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector 
    Eigen::VectorXd rhs(nxny);                // RHS vector

    double c_psi = c_psi_in;

    for (size_t i = 1; i < nx-1; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            size_t p_w = p_index(i - 1, j, ny);
            size_t p_0 = p_index(i    , j, ny);
            size_t p_e = p_index(i + 1, j, ny);
            h_xixi[p_0] = (h[p_w] - 2. * h[p_0] + h[p_e]);
            q_xixi[p_0] = (q[p_w] - 2. * q[p_0] + q[p_e]);
            s_xixi[p_0] = ((h[p_w] + zb[p_w]) - 2. * (h[p_0] + zb[p_0]) + (h[p_e] + zb[p_e]));
        }
    }
    for (size_t j = 0; j < ny; ++j)  // loop over west and east boundary
    {
        size_t i = 0;
        size_t p_0  = p_index(i    , j, ny);
        size_t p_e  = p_index(i + 1, j, ny);
        size_t p_ee = p_index(i + 2, j, ny);
        h_xixi[p_0] = 2. * h_xixi[p_e] - h_xixi[p_ee];
        q_xixi[p_0] = 2. * q_xixi[p_e] - q_xixi[p_ee];
        s_xixi[p_0] = 2. * s_xixi[p_e] - s_xixi[p_ee];

        i = nx - 1;
        p_0         = p_index(i    , j, ny);
        size_t p_w  = p_index(i - 1, j, ny);
        size_t p_ww = p_index(i - 2, j, ny);
        h_xixi[p_0] = 2. * h_xixi[p_w] - h_xixi[p_ww];
        q_xixi[p_0] = 2. * q_xixi[p_w] - q_xixi[p_ww];
        s_xixi[p_0] = 2. * s_xixi[p_w] - s_xixi[p_ww];
    }

    // eq. 18
    double hbar_im14;
    double hbar_ip14;
    double qbar_im14;
    double qbar_ip14;

    double c_error = 1.0; // same value as for regularization of given function
    for (size_t i = 1; i < nx - 1; ++i)
    {
        for (size_t j = 0; j < ny ; ++j)  // loop also over south and north boundary
        {
            size_t p_w = p_index(i - 1, j, ny);
            size_t p_0 = p_index(i    , j, ny);
            size_t p_e = p_index(i + 1, j, ny);

            A.coeffRef(p_0, p_w) = m_mass[0] - c_error;
            A.coeffRef(p_0, p_0) = m_mass[1] + 2. * c_error;
            A.coeffRef(p_0, p_e) = m_mass[2] - c_error;

            hbar_im14 = 0.25 * (h[p_w] + 3. * h[p_0]);
            hbar_ip14 = 0.25 * (h[p_e] + 3. * h[p_0]);
            qbar_im14 = 0.25 * (q[p_w] + 3. * q[p_0]);
            qbar_ip14 = 0.25 * (q[p_e] + 3. * q[p_0]);

            rhs[p_0] = 32.0 * c_psi * metric.dx_dxi[p_0] * (
                  0.0625 * std::sqrt(m_g / hbar_im14) * std::abs(s_xixi[p_0])
                + 0.0625 * std::sqrt(2.) * (std::abs(q_xixi[p_0] / hbar_im14 - qbar_im14 * h_xixi[p_0] / (hbar_im14 * hbar_im14)))
                + 0.0625 * std::sqrt(m_g / hbar_ip14) * std::abs(s_xixi[p_0])
                + 0.0625 * std::sqrt(2.) * (std::abs(q_xixi[p_0] / hbar_ip14 - qbar_ip14 * h_xixi[p_0] / (hbar_ip14 * hbar_ip14)))
                );
        }
    }
    // eq. 19
    // west and east boundary
    for (size_t j = 0; j < ny; ++j)  // loop also over south and north boundary
    {
        size_t p_ww;  
        size_t p_w;  
        size_t p_0;  
        size_t p_e;  
        size_t p_ee;  

        // west boundary
        size_t i = 0;
        p_0  = p_index(i    , j, ny);
        p_e  = p_index(i + 1, j, ny);
        p_ee = p_index(i + 2, j, ny);
        A.coeffRef(p_0, p_0 ) = 1.0;
        A.coeffRef(p_0, p_e ) = -2.0;
        A.coeffRef(p_0, p_ee) = 1.0;
        rhs[p_0] = 0.0;

        i = 1;
        p_w = p_index(i - 1, j, ny);
        p_0 = p_index(i    , j, ny);
        p_e = p_index(i + 1, j, ny);
        A.coeffRef(p_0, p_w) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_e) = 0.0;
        rhs[p_0] = rhs[p_0];

        // east boundary
        i = nx - 1;
        p_ww = p_index(i - 2, j, ny);
        p_w  = p_index(i - 1, j, ny);
        p_0  = p_index(i    , j, ny);
        A.coeffRef(p_0,  p_ww) = 1.0;
        A.coeffRef(p_0,  p_w ) = -2.0;
        A.coeffRef(p_0,  p_0 ) = 1.0;
        rhs[p_0] = 0.0;

        i = nx - 2;
        p_w = p_index(i - 1, j, ny);
        p_0 = p_index(i    , j, ny);
        p_e = p_index(i + 1, j, ny);
        A.coeffRef(p_0, p_w) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_e) = 0.0;
        rhs[p_0] = rhs[p_0];
    }
    if (m_logging == "matrix_visc")
    {
        std::string header_text = "=== Viscosity xi, matrix  =============================";
        log_file << std::showpos << std::setprecision(3) << std::scientific;
        print_matrix(A, 1, nx, ny, header_text, log_file);
        header_text = "=== Viscosity xi, RHS =================================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(rhs, 1, nx, ny, header_text, log_file);
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);
    if (m_logging == "matrix_visc")
    {
        std::string header_text = "=== Viscosity xi, Solution ============================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(solution, 1, nx, ny, header_text, log_file);
    }

    for (size_t i = 0; i < nxny; ++i)
    {
        psi[i] = solution[i];
    }
}
//------------------------------------------------------------------------------
void REGULARIZATION::artificial_viscosity_eta(std::vector<double>& psi, 
    std::vector<double>& h, std::vector<double>& r, std::vector<double>& zb, 
    double c_psi_in, struct _grid_metric & metric, std::ofstream& log_file)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;

    int status = 1; 

    std::vector<double> h_etaeta(nxny, 0.);  // second derivative of total depth in computational space
    std::vector<double> r_etaeta(nxny, 0.);  // second derivative of flow flux in computational space
    std::vector<double> s_etaeta(nxny, 0.);  // second derivative of water level: s = h + zb in computational space

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector 
    Eigen::VectorXd rhs(nxny);                // RHS vector

    double c_psi = c_psi_in;

    for (size_t i = 0; i < nx; ++i)  // loop also over west and east boundary
    {
        for (size_t j = 1; j < ny -1 ; ++j)
        {
            size_t p_s = p_index(i, j - 1, ny);
            size_t p_0 = p_index(i, j    , ny);
            size_t p_n = p_index(i, j + 1, ny);
            h_etaeta[p_0] = (h[p_s] - 2. * h[p_0] + h[p_n]);
            r_etaeta[p_0] = (r[p_s] - 2. * r[p_0] + r[p_n]);
            s_etaeta[p_0] = ((h[p_s] + zb[p_s]) - 2. * (h[p_0] + zb[p_0]) + (h[p_n] + zb[p_n]));
        }
    }
    for (size_t i = 1; i < nx - 1; ++i)  // loop over south and north boundary
    {
        size_t j = 0;
        size_t p_0  = p_index(i, j    , ny);
        size_t p_n  = p_index(i, j + 1, ny);
        size_t p_nn = p_index(i, j + 2, ny);
        h_etaeta[p_0] = 2. * h_etaeta[p_n] - h_etaeta[p_nn];
        r_etaeta[p_0] = 2. * r_etaeta[p_n] - r_etaeta[p_nn];
        s_etaeta[p_0] = 2. * s_etaeta[p_n] - s_etaeta[p_nn];

        j = ny - 1;
        p_0         = p_index(i, j    , ny);
        size_t p_s  = p_index(i, j - 1, ny);
        size_t p_ss = p_index(i, j - 2, ny);
        h_etaeta[p_0] = 2. * h_etaeta[p_s] - h_etaeta[p_ss];
        r_etaeta[p_0] = 2. * r_etaeta[p_s] - r_etaeta[p_ss];
        s_etaeta[p_0] = 2. * s_etaeta[p_s] - s_etaeta[p_ss];
    }

    // eq. 18
    double hbar_jm14;
    double hbar_jp14;
    double rbar_jm14;
    double rbar_jp14;

    double c_error = c_psi; // same value as for regularization of given function
    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 1; j < ny - 1; ++j)
        {
            size_t p_s = p_index(i, j - 1, ny);
            size_t p_0 = p_index(i, j    , ny);
            size_t p_n = p_index(i, j + 1, ny);
            
            A.coeffRef(p_0, p_s) = m_mass[0] - c_error;
            A.coeffRef(p_0, p_0) = m_mass[1] + 2. * c_error;
            A.coeffRef(p_0, p_n) = m_mass[2] - c_error;

            hbar_jm14 = 0.25 * (h[p_s] + 3. * h[p_0]);
            hbar_jp14 = 0.25 * (h[p_n] + 3. * h[p_0]);
            rbar_jm14 = 0.25 * (r[p_s] + 3. * r[p_0]);
            rbar_jp14 = 0.25 * (r[p_n] + 3. * r[p_0]);

            rhs[p_0] = 32.0 * c_psi * metric.dy_deta[p_0] * (
                  0.0625 * std::sqrt(m_g / hbar_jm14) * std::abs(s_etaeta[p_0])
                + 0.0625 * std::sqrt(2.) * (std::abs(r_etaeta[p_0] / hbar_jm14 - rbar_jm14 * h_etaeta[p_0] / (hbar_jm14 * hbar_jm14)))
                + 0.0625 * std::sqrt(m_g / hbar_jp14) * std::abs(s_etaeta[p_0])
                + 0.0625 * std::sqrt(2.) * (std::abs(r_etaeta[p_0] / hbar_jp14 - rbar_jp14 * h_etaeta[p_0] / (hbar_jp14 * hbar_jp14)))
                );
        }
    }
    // eq. 19
    for (size_t i = 0; i < nx; ++i)  // loop over south and north boundary
    {
        size_t j = 0;
        size_t p_ss;  
        size_t p_s;  
        size_t p_0;  
        size_t p_n;  
        size_t p_nn; 
        
        // south boundary
        j = 0;
        p_0  = p_index(i, j    , ny);
        p_n  = p_index(i, j + 1, ny);
        p_nn = p_index(i, j + 2, ny);
        A.coeffRef(p_0, p_0 ) = 1.0;
        A.coeffRef(p_0, p_n ) = -2.0;
        A.coeffRef(p_0, p_nn) = 1.0;
        rhs[p_0] = 0.0;

        j = 1;
        p_s = p_index(i, j - 1, ny);
        p_0 = p_index(i, j    , ny);
        p_n = p_index(i, j + 1, ny);
        A.coeffRef(p_0, p_s) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_n) = 0.0;
        rhs[p_0] = rhs[p_0];

        j = ny - 1;
        p_0  = p_index(i, j    , ny);
        p_s  = p_index(i, j - 1, ny);
        p_ss = p_index(i, j - 2, ny);
        A.coeffRef(p_0, p_0 ) =  1.0;
        A.coeffRef(p_0, p_s ) = -2.0;
        A.coeffRef(p_0, p_ss) =  1.0;
        rhs[p_0] = 0.0;

        j = ny - 2;
        p_s = p_index(i, j - 1, ny);
        p_0 = p_index(i, j    , ny);
        p_n = p_index(i, j + 1, ny);
        A.coeffRef(p_0, p_s) = 0.0;
        A.coeffRef(p_0, p_0) = 1.0;
        A.coeffRef(p_0, p_n) = 0.0;
        rhs[p_0] = rhs[p_0];
    }

    if (m_logging == "matrix_visc")
    {
        std::string header_text = "=== Viscosity eta, matrix =============================";
        log_file << std::showpos << std::setprecision(3) << std::scientific;
        print_matrix(A, 1, nx, ny, header_text, log_file);
        header_text = "=== Viscosity eta, RHS ================================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(rhs, 1, nx, ny, header_text, log_file);
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    solver.setTolerance(1e-12);
    solution = solver.solve(rhs);
    if (m_logging == "matrix_visc")
    {
        std::string header_text = "=== Viscosity eta, solution ===========================";

        log_file << std::setprecision(8) << std::scientific;
        print_vector(solution, 1, nx, ny, header_text, log_file);
    }

    for (size_t i = 0; i < nxny; ++i)
    {
        psi[i] = solution[i];
    }
}
//------------------------------------------------------------------------------
void REGULARIZATION::first_derivative( std::vector<double> & psi, std::vector<double>& eps,  std::vector<double>& u, struct _grid_metric & metric)
{
    size_t nx = eps.size();
    std::vector<double> pe(nx, 0.0);
    double pe_threshold = 1.8;

    for (size_t i = 0; i < nx; ++i)
    {
        pe[i] = std::abs(u[i] * metric.dx_dxi[i] / eps[i]);
    }
    for (size_t i = 0; i < nx; ++i)
    {
        psi[i] = 1.0;
        if (pe[i] > pe_threshold)
        {
            psi[i] = pe[i] / pe_threshold;
        }
    }
}
//------------------------------------------------------------------------------
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(struct _grid_metric & metric, 
    std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& u_giv, std::ofstream& log_file)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx * ny;
    auto u = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nxny, 0.0);

    Eigen::SparseMatrix<double, Eigen::RowMajor>  B(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector u
    Eigen::VectorXd rhs(nxny);                // RHS
    std::vector< Eigen::Triplet<double> > triplets; 
    triplets.reserve(9 * nxny); // 9-points per row, nrow = nxny; large enough to avoid reallocation
    
    int status = build_matrix_pattern_regularization(triplets, nx, ny);
    B.setFromTriplets(triplets.begin(), triplets.end());
    B.makeCompressed(); // Very important for valuePtr() access

    double* values = B.valuePtr();         // pointer to all non-zero values
    const int* outer = B.outerIndexPtr();   // row start pointers

    //
    // Add the diffusion term
    //
    // row 0, 1, 2; sw-corner
    // row 3,..., 3 * (ny - 1) - 1: west boundary
    // row 3 * (ny - 1), +1, +2; nw-corner
    // row std::fmod(col , 3 * ny), +1, +2; south boundary 
    // row std::fmod(col + 3, 3 * ny), +1, +2; north boundary
    // row 3 * (nx - 1) * ny, +1, +2; se-corner
    // row 3 * (nx - 1) * ny + 3, ..., 3 * nx * ny - 3 - 1: east boundary
    // row 3 * nx * ny - 3, +1, +2 : ne-corner
    START_TIMERN(Regularization solve eq7);
    //
    // viscosity
    //
    // For the moment only interior nodes (2025-08-13)
    // Do not clear the rows, but add the viscosity terms

    // south-west corner
    for (size_t row = 0; row < 1; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_south_west_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, metric);
    }
    // west boundary
    for (size_t row = 1; row < 1 * (ny - 1); row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        status = reg_boundary_west_utilde(values, row, c_eq, rhs,
            u_giv, psi_11[p_0], psi_22[p_0], 
            (double) 1.0, metric);
    }
    // north-west corner
    for (size_t row = 1 * (ny - 1); row < 1 * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_west_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, metric);
    }

    // interior with south and north boundary
    for (size_t row = 1 * ny; row < 1 * (nx - 1) * ny; row += 1) 
    {
        size_t c_eq = outer[row    ];
        size_t p_0 = c_eq/(9);

        if (row % ny == 0) {
            // south boundary, over write coefficients
            status = reg_boundary_south_utilde(values, row, c_eq, rhs, 
                u_giv, psi_11[p_0], psi_22[p_0], 
                (double) 1.0, metric);
            continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary, over write coefficients
            status = reg_boundary_north_utilde(values, row, c_eq, rhs, 
                u_giv, psi_11[p_0], psi_22[p_0], 
                (double) 1.0, metric);
            continue;
        }

        status = utilde_interior_matrix(values, row, c_eq, rhs, 
            psi_11[p_0], psi_22[p_0], metric);
        status = utilde_interior_rhs(values, row, c_eq, rhs,
            u_giv, metric);
    }
    // south-east corner
    for (size_t row = 1 * (nx - 1) * ny; row < 1 * (nx - 1) * ny + 1; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_south_east_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, metric);
    }
    // east boundary
    for (size_t row = 1 * (nx - 1) * ny + 1; row < 1 * nx * ny - 1; row += 1) 
    {
        int c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        status = reg_boundary_east_utilde(values, row, c_eq, rhs, 
            u_giv, psi_11[p_0], psi_22[p_0], 
            (double) 1.0, metric);
    }
    // north-east corner
    for (size_t row = 1 * nx * ny - 1; row < 1 * nx * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_east_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, metric);
    }

    if (m_logging == "pattern")
    {
        std::string header_text = "=== Matrix eq7 ========================================";
        log_file << std::showpos << std::setprecision(3) << std::scientific;
        print_matrix(B, 1, nx, ny, header_text, log_file);
        header_text = "=== RHS eq7 ===========================================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(rhs, 1, nx, ny, header_text, log_file);
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverB;
    solverB.compute(B);
    solution = solverB.solve(rhs);
    //solution = solverB.solveWithGuess(rhs, solution);

    STOP_TIMER(Regularization solve eq7);

    for (size_t i = 0; i < nxny; ++i)
    {
        u->push_back(solution[i]);
    }

    return u;
}
//------------------------------------------------------------------------------
std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(struct _grid_metric & metric,
     double c_error, std::vector<double>& u0, std::vector<double>& u0_xixi, std::vector<double>& u0_etaeta, 
     std::ofstream& log_file)
{
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    size_t nxny = nx*ny;
    auto eq8 = std::make_unique<std::vector<double>> ();
    std::vector<double> tmp(nxny, NAN);

    Eigen::SparseMatrix<double, Eigen::RowMajor>  A(nxny, nxny);
    Eigen::VectorXd solution(nxny);           // solution vector u
    Eigen::VectorXd rhs(nxny);                // RHS
    std::vector< Eigen::Triplet<double> > triplets; 
    triplets.reserve(9 * nxny); // 9-points per row, nrow = nxny; large enough to avoid reallocation
    
    int status = build_matrix_pattern_regularization(triplets, nx, ny);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed(); // Very important for valuePtr() access

    if (m_logging == "pattern")
    {
        std::string header_text = "=== Matrix build matrix pattern (solve eq8) ===========";
        print_matrix_pattern(A, 1, nx, ny, header_text,log_file);
    }

    double u0_xixi_max = 0.0;
    double u0_etaeta_max = 0.0;
    //for (size_t i = 0; i < nxny; ++i)
    //{
    //    u0_xixi_max = std::max(u0_xixi_max, u0_xixi[i]);
    //    u0_etaeta_max = std::max(u0_etaeta_max, u0_etaeta[i]);
    //}
    //if (u0_xixi_max == 0.0) { u0_xixi_max = 1.0; }
    //if (u0_etaeta_max == 0.0) { u0_etaeta_max = 1.0; }
    u0_xixi_max = 1.0;
    u0_etaeta_max = 1.0;

    double* values = A.valuePtr();         // pointer to all non-zero values
    const int* outer = A.outerIndexPtr();   // row start pointers

    size_t c_eq;
    size_t p0;

    // south-west corner
    for (size_t row = 0; row < 1; row += 1)
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_sw

        values[c_eq     ] = -2.0;
        values[c_eq +  1] = 1.0;
        values[c_eq +  2] = 0.0;
        values[c_eq +  3] = 1.0;
        values[c_eq +  4] = 0.0;
        values[c_eq +  5] = 0.0;
        values[c_eq +  6] = 0.0;
        values[c_eq +  7] = 0.0;
        values[c_eq +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // west boundary
    for (size_t row = 1; row < (ny - 1); row += 1)
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_w

        values[c_eq     ] = 0.0;
        values[c_eq +  1] = -1.0;
        values[c_eq +  2] = 0.0;
        values[c_eq +  3] = 0.0;
        values[c_eq +  4] = 1.0;
        values[c_eq +  5] = 0.0;
        values[c_eq +  6] = 0.0;
        values[c_eq +  7] = 0.0;
        values[c_eq +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-west corner
    for (size_t row = ny - 1; row < ny; row += 1)
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_nw

        values[c_eq     ] = 0.0;
        values[c_eq +  1] = 1.0;
        values[c_eq +  2] = -2.0;
        values[c_eq +  3] = 0.0;
        values[c_eq +  4] = 0.0;
        values[c_eq +  5] = 1.0;
        values[c_eq +  6] = 0.0;
        values[c_eq +  7] = 0.0;
        values[c_eq +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // interior with south and north boundary
    for (size_t row = 1 * ny; row < (nx - 1) * ny; row += 1) 
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_0

        if (row % ny == 0) {
            // south boundary
            c_eq = outer[row    ];
            p0 = c_eq/(9);  // p_s

            values[c_eq     ] = 0.0;
            values[c_eq +  1] = 0.0;
            values[c_eq +  2] = 0.0;
            values[c_eq +  3] = -1.0;
            values[c_eq +  4] = 1.0;
            values[c_eq +  5] = 0.0;
            values[c_eq +  6] = 0.0;
            values[c_eq +  7] = 0.0;
            values[c_eq +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary
            c_eq = outer[row    ];
            p0 = c_eq/(9);  // p_n

            values[c_eq     ] = 0.0;
            values[c_eq +  1] = 0.0;
            values[c_eq +  2] = 0.0;
            values[c_eq +  3] = 0.0;
            values[c_eq +  4] = -1.0;
            values[c_eq +  5] = 1.0;
            values[c_eq +  6] = 0.0;
            values[c_eq +  7] = 0.0;
            values[c_eq +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }

        size_t p_sw = p0 - ny - 1;
        size_t p_w  = p0 - ny    ;
        size_t p_nw = p0 - ny + 1;
        size_t p_s  = p0 - 1;
        size_t p_0  = p0    ;
        size_t p_n  = p0 + 1;
        size_t p_se = p0 + ny - 1;
        size_t p_e  = p0 + ny    ;
        size_t p_ne = p0 + ny + 1;

        //values[c_eq     ] = m_mass[0] * m_mass[0];
        //values[c_eq +  1] = m_mass[1] * m_mass[0] - c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[c_eq +  2] = m_mass[0] * m_mass[2];
        //values[c_eq +  3] = m_mass[0] * m_mass[1] - c_error * u0_xixi[p_w] / u0_xixi_max;
        //values[c_eq +  4] = m_mass[1] * m_mass[1] + 2. * c_error * u0_xixi[p_w] / u0_xixi_max + 2. * c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[c_eq +  5] = m_mass[2] * m_mass[1] - c_error * u0_xixi[p_w] / u0_xixi_max;
        //values[c_eq +  6] = m_mass[2] * m_mass[0];
        //values[c_eq +  7] = m_mass[1] * m_mass[2] - c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[c_eq +  8] = m_mass[2] * m_mass[2];
        
        values[c_eq     ] = m_mass[0] * m_mass[0] - 1.0/6.0 * c_error;
        values[c_eq +  1] = m_mass[0] * m_mass[1] - 4.0/6.0 * c_error;
        values[c_eq +  2] = m_mass[0] * m_mass[2] - 1.0/6.0 * c_error;
        values[c_eq +  3] = m_mass[1] * m_mass[0] - 4.0/6.0 * c_error;
        values[c_eq +  4] = m_mass[1] * m_mass[1] + 20.0/6.0 * c_error;
        values[c_eq +  5] = m_mass[1] * m_mass[2] - 4.0/6.0 * c_error;
        values[c_eq +  6] = m_mass[2] * m_mass[0] - 1.0/6.0 * c_error;
        values[c_eq +  7] = m_mass[2] * m_mass[1] - 4.0/6.0 * c_error;
        values[c_eq +  8] = m_mass[2] * m_mass[2] - 1.0/6.0 * c_error;

        rhs[row]  =  1./64. * (u0_xixi[p_sw] + u0_etaeta[p_sw])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] +=  6./64. * (u0_xixi[p_w ] + u0_etaeta[p_w ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] +=  1./64. * (u0_xixi[p_nw] + u0_etaeta[p_nw])/(u0_xixi_max + u0_etaeta_max);
        
        rhs[row] +=  6./64. * (u0_xixi[p_s ] + u0_etaeta[p_s ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += 36./64. * (u0_xixi[p_0 ] + u0_etaeta[p_0 ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] +=  6./64. * (u0_xixi[p_n ] + u0_etaeta[p_n ])/(u0_xixi_max + u0_etaeta_max);
        
        rhs[row] +=  1./64. * (u0_xixi[p_se] + u0_etaeta[p_se])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] +=  6./64. * (u0_xixi[p_e ] + u0_etaeta[p_e ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] +=  1./64. * (u0_xixi[p_ne] + u0_etaeta[p_ne])/(u0_xixi_max + u0_etaeta_max);
        //rhs[row] = std::abs((u0_xixi[p_0 ] + u0_etaeta[p_0 ])/(u0_xixi_max + u0_etaeta_max));
    }
    // south-east corner
    for (size_t row = (nx - 1) * ny; row < (nx - 1) * ny + 1; row += 1)
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_se

        values[c_eq     ] = 0.0;
        values[c_eq +  1] = 0.0;
        values[c_eq +  2] = 0.0;
        values[c_eq +  3] = -1.0;
        values[c_eq +  4] = 0.0;
        values[c_eq +  5] = 0.0;
        values[c_eq +  6] = 2.0;
        values[c_eq +  7] = -1.0;
        values[c_eq +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // east boundary
    for (size_t row = (nx - 1) * ny + 1; row < nx * ny - 1; row += 1) 
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_e

        values[c_eq     ] = 0.0;
        values[c_eq +  1] = 0.0;
        values[c_eq +  2] = 0.0;
        values[c_eq +  3] = 0.0;
        values[c_eq +  4] = -1.0;
        values[c_eq +  5] = 0.0;
        values[c_eq +  6] = 0.0;
        values[c_eq +  7] = 1.0;
        values[c_eq +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-east corner
    for (size_t row = nx * ny - 1; row < nx * ny; row += 1)
    {
        c_eq = outer[row    ];
        p0 = c_eq/(9);  // p_ne

        values[c_eq     ] = 0.0;
        values[c_eq +  1] = 0.0;
        values[c_eq +  2] = 0.0;
        values[c_eq +  3] = 0.0;
        values[c_eq +  4] = 0.0;
        values[c_eq +  5] = -1.0;
        values[c_eq +  6] = 0.0;
        values[c_eq +  7] = -1.0;
        values[c_eq +  8] = 2.0;
        rhs[row    ] = 0.0;
    }

    for (size_t i = 0; i < nxny; ++i)
    {
        solution[i] = u0[i];
    }

    if (m_logging == "matrix")
    {
        std::string header_text = "=== Matrix eq8 ========================================";
        log_file << std::showpos << std::setprecision(3) << std::scientific;
        print_matrix(A, 1, nx, ny, header_text, log_file);
        header_text = "=== RHS eq8 ===========================================";
        log_file << std::setprecision(8) << std::scientific;
        print_vector(rhs, 1, nx, ny, header_text, log_file);
    }

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
    solverA.compute(A);
    solverA.setTolerance(1e-12);
    solution = solverA.solve(rhs);
    //solution = solverA.solveWithGuess(rhs, solution);

    for (size_t i = 0; i < nxny; ++i)
    {
        eq8->push_back(solution[i]); // / (err_max + 1e-13);  // to prevent division bij zero
    }
    return eq8;
}
inline void REGULARIZATION::add_value(double * values, size_t col, double data)
{ 
    values[col] += data; 
}
size_t REGULARIZATION::p_index(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}
