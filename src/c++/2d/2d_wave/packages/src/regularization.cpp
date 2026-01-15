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
    std::vector<double>& u_tilde, std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& eq8,
    std::vector<double>& x, std::vector<double>& y,
    std::vector<double>& u_giv_in,
    size_t nx, size_t ny, double c_psi, std::ofstream& log_file)
{
    double diff_max0 = 0.0;
    double diff_max1 = 0.0;
    size_t nxny = nx * ny;
    std::vector<double> u_giv(nxny, 0.);
    std::vector<double> u0(nxny, 0.);
    std::vector<double> u1(nxny, 0.);
    std::vector<double> u0_xixi(nxny, 0.);
    std::vector<double> u0_etaeta(nxny, 0.);
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
        double u0_xixi_max = 0.0;
        double u0_etaeta_max = 0.0;
        for (size_t j = 1; j < ny - 1; ++j)
        {
            for (size_t i = 1; i < nx - 1; ++i)
            {
                size_t p_0 = p_index(i, j    , ny);
                size_t p_n = p_index(i, j + 1, ny);
                size_t p_e = p_index(i + 1, j, ny);
                size_t p_s = p_index(i, j - 1, ny);
                size_t p_w = p_index(i - 1, j, ny);
                u0_xixi[p_0]   = std::abs((u0[p_e] - 2. * u0[p_0] + u0[p_w]));
                u0_etaeta[p_0] = std::abs((u0[p_n] - 2. * u0[p_0] + u0[p_s]));
                u0_xixi_max = std::max(u0_xixi_max, std::abs(u0_xixi[p_0]));
                u0_etaeta_max = std::max(u0_etaeta_max, std::abs(u0_etaeta[p_0]));
            }
        }
        double smooth= std::max(u0_xixi_max, u0_etaeta_max);
        if (smooth < m_u0_is_smooth)
        {
            for (size_t i = 0; i < nxny; ++i)
            {
                //u_out[i] = u_giv_in[i];
            }
            //return;
            smooth = smooth + 1.e-12;
        }

        for (size_t i = 0; i < nx; ++i)  // horizontal direction
        {
            size_t j = 0;
            size_t p_0  = p_index(i, j    , ny);
            size_t p_n  = p_index(i, j + 1, ny);
            size_t p_nn = p_index(i, j + 2, ny);
            u0_etaeta[p_0] = u0_etaeta[p_nn];
            u0_etaeta[p_n] = u0_etaeta[p_nn];

            j = ny - 1;
            p_0         = p_index(i, j    , ny);
            size_t p_s  = p_index(i, j - 1, ny);
            size_t p_ss = p_index(i, j - 2, ny);
            u0_etaeta[p_0] = u0_etaeta[p_ss];
            u0_etaeta[p_s] = u0_etaeta[p_ss];
        }
        for (size_t j = 0; j < ny; ++j)  // vertical direction
        {
            size_t i = 0;
            size_t p_0  = p_index(i    , j, ny);
            size_t p_e  = p_index(i + 1, j, ny);
            size_t p_ee = p_index(i + 2, j, ny);
            u0_xixi[p_0] = u0_xixi[p_ee];
            u0_xixi[p_e] = u0_xixi[p_ee];

            i = nx - 1;
            p_0      = p_index(i    , j, ny);
            size_t p_w  = p_index(i - 1, j, ny);
            size_t p_ww = p_index(i - 2, j, ny);
            u0_xixi[p_0] = u0_xixi[p_ww];
            u0_xixi[p_w] = u0_xixi[p_ww];
        }

//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(nx, ny, x, y, c_psi, u0, u0_xixi, u0_etaeta, log_file));  // eq8 is computed in computational space
//------------------------------------------------------------------------------
        for (size_t i = 0; i < nxny; ++i)
        {
            double dxi = 1.0;
            double deta = 1.0;
            psi_11[i] = c_psi * (dxi * dxi + deta * deta) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
            psi_22[i] = c_psi * (dxi * dxi + deta * deta) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(nx, ny, x, y, psi_11, psi_22, u_giv, log_file));
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
            m_u0_is_smooth = std::max(u0_xixi_max, u0_etaeta_max);
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
    std::vector<double>& x, std::vector<double>& y, size_t nx, size_t ny, 
    double c_psi_in, std::ofstream& log_file)
{
    int status = 0;
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
            x,  y,
            h, psi_1, psi_2, 
            1.0, nx, ny);
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
                x,  y,
                h, psi_1, psi_2, 
                1.0, nx, ny);
                continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary, over write coefficients
            status = reg_boundary_north_psi( values, row, c_eq, rhs,
                x,  y,
                h, psi_1, psi_2, 
                1.0, nx, ny);
                continue;
        }
        status = reg_interior_matrix_psi( values, row, c_eq,
             c_psi,  x,  y, nx, ny);
        // overwrite the right hand side
        status = reg_interior_rhs_psi(row, c_eq, rhs, 
            h, q, r, x,  y, c_psi, m_g, nx, ny);

        rhs[row] = 
            0.0
            ;
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
        int c_eq = (size_t) outer[row    ];
        status = 0;
        size_t p_0 = c_eq/(9);

        status = reg_boundary_east_psi( values, row, c_eq, rhs,
            x,  y,
            h, psi_1, psi_2, 
            1.0, nx, ny);
            continue;
    }
    // north-east corner
    for (size_t row = 1 * nx * ny - 1; row < 1 * nx * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_east_psi( values, row, c_eq, rhs, 
            h, 1.0, nx, ny);
    }
    STOP_TIMER(Regularization diffusion);

    Eigen::BiCGSTAB< Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverA;
    solverA.compute(A);
    solverA.setTolerance(1e-12);
    solution = solverA.solve(rhs);

    for (size_t i = 1; i < nxny - 1; ++i)
    {
        psi[i] = solution[i];
    }

    return;
}
//------------------------------------------------------------------------------
void REGULARIZATION::first_derivative( std::vector<double> & psi, std::vector<double>& eps,  std::vector<double>& u, double dx)
{
    size_t nx = eps.size();
    std::vector<double> pe(nx, 0.0);
    double pe_threshold = 1.8;

    for (size_t i = 0; i < nx; ++i)
    {
        pe[i] = std::abs(u[i] * dx / eps[i]);
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
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y, 
    std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& u_giv, std::ofstream& log_file)
{
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
    START_TIMERN(Regularization diffusion);
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
           u_giv, (double) 1.0, nx, ny);
    }
    // west boundary
    for (size_t row = 1; row < 1 * (ny - 1); row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = -psi_11[p_0];
        double psi_2 = -psi_22[p_0];
        status = reg_boundary_west_utilde(values, row, c_eq, rhs,
            x, y, u_giv, psi_1, psi_2, 
            (double) 1.0, nx, ny);
    }
    // north-west corner
    for (size_t row = 1 * (ny - 1); row < 1 * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_west_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }

    // interior with south and north boundary
    for (size_t row = 1 * ny; row < 1 * (nx - 1) * ny; row += 1) 
    {
        size_t c_eq = outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = -psi_11[p_0];
        double psi_2 = -psi_22[p_0];

        if (row % ny == 0) {
            // south boundary, over write coefficients
            status = reg_boundary_south_utilde(values, row, c_eq, rhs, 
                x, y, u_giv, psi_1, psi_2, 
                (double) 1.0, nx, ny);
            continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary, over write coefficients
            status = reg_boundary_north_utilde(values, row, c_eq, rhs, 
                x, y, u_giv, psi_1, psi_2, 
                (double) 1.0, nx, ny);
            continue;
        }

        status = reg_interior_utilde(values, row, c_eq, rhs, u_giv, x, y, nx, ny);
        status = diffusion_matrix_and_rhs(values, row, c_eq, rhs,
                    x, y, u_giv, psi_1, psi_2, (double) 1.0, nx, ny);
    }
    // south-east corner
    for (size_t row = 1 * (nx - 1) * ny; row < 1 * (nx - 1) * ny + 1; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_south_east_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }
    // east boundary
    for (size_t row = 1 * (nx - 1) * ny + 1; row < 1 * nx * ny - 1; row += 1) 
    {
        int c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = -psi_11[p_0];
        double psi_2 = -psi_22[p_0];

        status = reg_boundary_east_utilde(values, row, c_eq, rhs, 
            x, y, u_giv, psi_1, psi_2, 
            (double) 1.0, nx, ny);
    }
    // north-east corner
    for (size_t row = 1 * nx * ny - 1; row < 1 * nx * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_east_utilde(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }
    STOP_TIMER(Regularization diffusion);


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

    for (size_t i = 0; i < nxny; ++i)
    {
        u->push_back(solution[i]);
    }

    return u;
}
//------------------------------------------------------------------------------
std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y,
     double c_error, std::vector<double>& u0, std::vector<double>& u0_xixi, std::vector<double>& u0_etaeta, 
    std::ofstream& log_file)
{
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
    for (size_t i = 0; i < nxny; ++i)
    {
        u0_xixi_max = std::max(u0_xixi_max, u0_xixi[i]);
        u0_etaeta_max = std::max(u0_etaeta_max, u0_etaeta[i]);
    }
    if (u0_xixi_max == 0.0) { u0_xixi_max = 1.0; }
    if (u0_etaeta_max == 0.0) { u0_etaeta_max = 1.0; }
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

        rhs[row]  = (u0_xixi[p_sw] + u0_etaeta[p_sw])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_w ] + u0_etaeta[p_w ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_nw] + u0_etaeta[p_nw])/(u0_xixi_max + u0_etaeta_max);
        
        rhs[row] += (u0_xixi[p_s ] + u0_etaeta[p_s ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_0 ] + u0_etaeta[p_0 ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_n ] + u0_etaeta[p_n ])/(u0_xixi_max + u0_etaeta_max);
        
        rhs[row] += (u0_xixi[p_se] + u0_etaeta[p_se])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_e ] + u0_etaeta[p_e ])/(u0_xixi_max + u0_etaeta_max);
        rhs[row] += (u0_xixi[p_ne] + u0_etaeta[p_ne])/(u0_xixi_max + u0_etaeta_max);
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
//------------------------------------------------------------------------------
int REGULARIZATION::reg_interior_rhs_psi( size_t row, size_t c_eq, Eigen::VectorXd& rhs, 
    std::vector<double>& h, std::vector<double>& q, std::vector<double>& r,
    std::vector<double>& x, std::vector<double>& y,
    double c_psi, double g, size_t nx, size_t ny)
{
    size_t nxny = nx * ny;
    std::vector<double> Err_psi(nxny, 0.0);

    std::vector<double> d2h_dxi2(nxny, 0.0);
    std::vector<double> d2h_dxideta(nxny, 0.0);
    std::vector<double> d2h_deta2(nxny, 0.0);
    std::vector<double> d2q_dxi2(nxny, 0.0);
    std::vector<double> d2q_dxideta(nxny, 0.0);
    std::vector<double> d2q_deta2(nxny, 0.0);
    std::vector<double> d2r_dxi2(nxny, 0.0);
    std::vector<double> d2r_dxideta(nxny, 0.0);
    std::vector<double> d2r_deta2(nxny, 0.0);
    std::vector<double> d2s_dxi2(nxny, 0.0);
    std::vector<double> d2s_dxideta(nxny, 0.0);
    std::vector<double> d2s_deta2(nxny, 0.0);
    std::vector<double> s(nxny, 0.0);
    std::vector<size_t> p(9);

    for (size_t i = 1; i < nx - 1; ++i)
    {
        for (size_t j = 1; j < ny - 1; ++j)
        {
            p[0] = p_index(i - 1, j - 1, ny);
            p[1] = p_index(i - 1, j    , ny);
            p[2] = p_index(i - 1, j + 1, ny);
            p[3] = p_index(i    , j - 1, ny);
            p[4] = p_index(i    , j    , ny);
            p[5] = p_index(i    , j + 1, ny);
            p[6] = p_index(i + 1, j - 1, ny);
            p[7] = p_index(i + 1, j    , ny);
            p[8] = p_index(i + 1, j + 1, ny);
            d2h_dxi2[p[4]]    = d2udxi2(h, p);
            d2q_dxi2[p[4]]    = d2udxi2(q, p);
            d2r_dxi2[p[4]]    = d2udxi2(r, p);
            d2s_dxi2[p[4]]    = d2udxi2(s, p);
            d2h_dxideta[p[4]] = d2udxideta(h, p);
            d2q_dxideta[p[4]] = d2udxideta(q, p);
            d2r_dxideta[p[4]] = d2udxideta(r, p);
            d2s_dxideta[p[4]] = d2udxideta(s, p);
            d2h_deta2[p[4]]   = d2udeta2(h, p);
            d2q_deta2[p[4]]   = d2udeta2(q, p);
            d2r_deta2[p[4]]   = d2udeta2(r, p);
            d2s_deta2[p[4]]   = d2udeta2(s, p);
        }
    }

    for (size_t i = 0; i < nx; ++i)  // south and north boundary
    {
        size_t j = 0;  // south boundary
        size_t p_0  = p_index(i, j    , ny);
        size_t p_n  = p_index(i, j + 1, ny);
        size_t p_nn = p_index(i, j + 2, ny);
        d2h_dxi2[p_0]    = 2.0 * d2h_dxi2[p_n]    - d2h_dxi2[p_nn];
        d2q_dxi2[p_0]    = 2.0 * d2q_dxi2[p_n]    - d2q_dxi2[p_nn];
        d2r_dxi2[p_0]    = 2.0 * d2r_dxi2[p_n]    - d2r_dxi2[p_nn];
        d2s_dxi2[p_0]    = 2.0 * d2s_dxi2[p_n]    - d2s_dxi2[p_nn];
        d2h_dxideta[p_0] = 2.0 * d2h_dxideta[p_n] - d2h_dxideta[p_nn];
        d2q_dxideta[p_0] = 2.0 * d2q_dxideta[p_n] - d2q_dxideta[p_nn];
        d2r_dxideta[p_0] = 2.0 * d2r_dxideta[p_n] - d2r_dxideta[p_nn];
        d2s_dxideta[p_0] = 2.0 * d2s_dxideta[p_n] - d2s_dxideta[p_nn];
        d2h_deta2[p_0]   = 2.0 * d2h_deta2[p_n]   - d2h_deta2[p_nn];
        d2q_deta2[p_0]   = 2.0 * d2q_deta2[p_n]   - d2q_deta2[p_nn];
        d2r_deta2[p_0]   = 2.0 * d2r_deta2[p_n]   - d2r_deta2[p_nn]; 
        d2s_deta2[p_0]   = 2.0 * d2s_deta2[p_n]   - d2s_deta2[p_nn]; 

        j = ny - 1;  // north boundary
        p_0         = p_index(i, j    , ny);
        size_t p_s  = p_index(i, j - 1, ny);
        size_t p_ss = p_index(i, j - 2, ny);
        d2h_dxi2[p_0]    = 2.0 * d2h_dxi2[p_s]    - d2h_dxi2[p_ss];
        d2q_dxi2[p_0]    = 2.0 * d2q_dxi2[p_s]    - d2q_dxi2[p_ss];
        d2r_dxi2[p_0]    = 2.0 * d2r_dxi2[p_s]    - d2r_dxi2[p_ss];
        d2s_dxi2[p_0]    = 2.0 * d2s_dxi2[p_s]    - d2s_dxi2[p_ss];
        d2h_dxideta[p_0] = 2.0 * d2h_dxideta[p_s] - d2h_dxideta[p_ss];
        d2q_dxideta[p_0] = 2.0 * d2q_dxideta[p_s] - d2q_dxideta[p_ss];
        d2r_dxideta[p_0] = 2.0 * d2r_dxideta[p_s] - d2r_dxideta[p_ss];
        d2s_dxideta[p_0] = 2.0 * d2s_dxideta[p_s] - d2s_dxideta[p_ss];
        d2h_deta2[p_0]   = 2.0 * d2h_deta2[p_s]   - d2h_deta2[p_ss];
        d2q_deta2[p_0]   = 2.0 * d2q_deta2[p_s]   - d2q_deta2[p_ss];
        d2r_deta2[p_0]   = 2.0 * d2r_deta2[p_s]   - d2r_deta2[p_ss]; 
        d2s_deta2[p_0]   = 2.0 * d2s_deta2[p_s]   - d2s_deta2[p_ss]; 
    }
    for (size_t j = 0; j < ny; ++j)  // west and east boundary
    {
        size_t i = 0;  // west boundary
        size_t p_0  = p_index(i    , j, ny);
        size_t p_e  = p_index(i + 1, j, ny);
        size_t p_ee = p_index(i + 2, j, ny);
        d2h_dxi2[p_0]    = 2.0 * d2h_dxi2[p_e]    - d2h_dxi2[p_ee];
        d2q_dxi2[p_0]    = 2.0 * d2q_dxi2[p_e]    - d2q_dxi2[p_ee];
        d2r_dxi2[p_0]    = 2.0 * d2r_dxi2[p_e]    - d2r_dxi2[p_ee];
        d2s_dxi2[p_0]    = 2.0 * d2s_dxi2[p_e]    - d2s_dxi2[p_ee];
        d2h_dxideta[p_0] = 2.0 * d2h_dxideta[p_e] - d2h_dxideta[p_ee];
        d2q_dxideta[p_0] = 2.0 * d2q_dxideta[p_e] - d2q_dxideta[p_ee];
        d2r_dxideta[p_0] = 2.0 * d2r_dxideta[p_e] - d2r_dxideta[p_ee];
        d2s_dxideta[p_0] = 2.0 * d2s_dxideta[p_e] - d2s_dxideta[p_ee];
        d2h_deta2[p_0]   = 2.0 * d2h_deta2[p_e]   - d2h_deta2[p_ee];
        d2q_deta2[p_0]   = 2.0 * d2q_deta2[p_e]   - d2q_deta2[p_ee];
        d2r_deta2[p_0]   = 2.0 * d2r_deta2[p_e]   - d2r_deta2[p_ee]; 
        d2s_deta2[p_0]   = 2.0 * d2s_deta2[p_e]   - d2s_deta2[p_ee]; 

        i = nx - 1;  // east boundary
        p_0      = p_index(i    , j, ny);
        size_t p_w  = p_index(i - 1, j, ny);
        size_t p_ww = p_index(i - 2, j, ny);
        d2h_dxi2[p_0]    =  2.0 * d2h_dxi2[p_w]    - d2h_dxi2[p_ww];
        d2q_dxi2[p_0]    =  2.0 * d2q_dxi2[p_w]    - d2q_dxi2[p_ww];
        d2r_dxi2[p_0]    =  2.0 * d2r_dxi2[p_w]    - d2r_dxi2[p_ww];
        d2s_dxi2[p_0]    =  2.0 * d2s_dxi2[p_w]    - d2s_dxi2[p_ww];
        d2h_dxideta[p_0] =  2.0 * d2h_dxideta[p_w] - d2h_dxideta[p_ww];
        d2q_dxideta[p_0] =  2.0 * d2q_dxideta[p_w] - d2q_dxideta[p_ww];
        d2r_dxideta[p_0] =  2.0 * d2r_dxideta[p_w] - d2r_dxideta[p_ww];
        d2s_dxideta[p_0] =  2.0 * d2s_dxideta[p_w] - d2s_dxideta[p_ww];
        d2h_deta2[p_0]   =  2.0 * d2h_deta2[p_w]   - d2h_deta2[p_ww];
        d2q_deta2[p_0]   =  2.0 * d2q_deta2[p_w]   - d2q_deta2[p_ww];
        d2r_deta2[p_0]   =  2.0 * d2r_deta2[p_w]   - d2r_deta2[p_ww]; 
        d2s_deta2[p_0]   =  2.0 * d2s_deta2[p_w]   - d2s_deta2[p_ww]; 
    }

    rhs.setZero();

    // Erro based on potential energy
    //
    //   \sqrt{g \widehat{h}}
    // 
    for (size_t i = 1; i < nx - 1; ++i)
    {
        for (size_t j = 1; j < ny - 1; ++j)
        {
            p[0] = p_index(i - 1, j - 1, ny);
            p[1] = p_index(i - 1, j    , ny);
            p[2] = p_index(i - 1, j + 1, ny);
            p[3] = p_index(i    , j - 1, ny);
            p[4] = p_index(i    , j    , ny);
            p[5] = p_index(i    , j + 1, ny);
            p[6] = p_index(i + 1, j - 1, ny);
            p[7] = p_index(i + 1, j    , ny);
            p[8] = p_index(i + 1, j + 1, ny);

            double f1 = F1(h, p, x, y, nx, ny);
            double f2 = F2(h, p, x, y, nx, ny);
            double f3 = F3(h, p, x, y, nx, ny);
            rhs[p[4]] = std::sqrt( g/h[p[4]] ) * (
                1.0/16.0 * f1 + 1.0/8.0 * f2 + 1.0/16.0 * f3
                );
        }
    }
    return 0;
}

inline double REGULARIZATION::d2udxi2(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    double dxi = 1.0;
    double deta = 1.0;
    double retval = 0.0;

    retval = 1.0 * u[p[0]] +
             4.0 * u[p[1]] +
             1.0 * u[p[2]] +
             4.0 * u[p[3]] +
            -20.0 * u[p[4]] +
             4.0 * u[p[5]] +
             1.0 * u[p[6]] +
             4.0 * u[p[7]] +
             1.0 * u[p[8]];

    return retval/6.0;
}
inline double REGULARIZATION::d2udxideta(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    double dxi = 1.0;
    double deta = 1.0;
    double retval = 0.0;

    retval = 1.0 * u[p[0]] +
             0.0 * u[p[1]] +
            -1.0 * u[p[2]] +
             0.0 * u[p[3]] +
             0.0 * u[p[4]] +
             0.0 * u[p[5]] +
            -1.0 * u[p[6]] +
             0.0 * u[p[7]] +
             1.0 * u[p[8]];

    return retval/(4.0 * dxi * deta);
}
inline double REGULARIZATION::d2udeta2(std::vector<double> & u, std::vector<size_t>& p)
{
    // Computational space
    return d2udxi2(u, p);
}
inline double REGULARIZATION::F1(std::vector<double> & u, std::vector<size_t>& p, 
    std::vector<double> & x, std::vector<double> &y, 
    size_t nx, size_t ny )
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dx_deta = 1.0;
    double dxi_dx = 1.0;
    double deta_dx = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid
    double d2xi_dx2 = 0.0;
    double d2eta_dx2 = 0.0;

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dx_dxi * dx_dxi * (
          dxi_dx * dxi_dx * d2u_dxi2 
        + 2.0 * dxi_dx * deta_dx * d2u_dxideta 
        + deta_dx * deta_dx * d2u_deta2 
        + d2xi_dx2 * du_dxi 
        + d2eta_dx2 * du_deta
        );

    return retval;
}
inline double REGULARIZATION::F2(std::vector<double> & u, std::vector<size_t>& p,
    std::vector<double> & x, std::vector<double> &y, 
    size_t nx, size_t ny )
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dy_deta = 1.0;
    double dxi_dx = 1.0;
    double dxi_dy = 1.0;
    double deta_dx = 1.0;
    double deta_dy = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dx_dxi * dy_deta * (
        dxi_dx * dxi_dy * d2u_dxi2 
        + (dxi_dx * deta_dy + deta_dx * dxi_dy) * d2u_dxideta 
        + deta_dx * deta_dy * d2u_deta2
        + d2xi_dxdy * du_dxi
        + d2eta_dxdy * du_deta
        );

        return 0.0;
}
inline double REGULARIZATION::F3(std::vector<double> & u, std::vector<size_t>& p, 
    std::vector<double> & x, std::vector<double> &y, 
    size_t nx, size_t ny )
{
    double retval = 0.0;

    double dx_dxi = 1.0;
    double dy_deta = 1.0;
    double dxi_dx = 1.0;
    double dxi_dy = 1.0;
    double deta_dx = 1.0;
    double deta_dy = 1.0;
    double d2xi_dxdy = 0.0;  // Assume no cuvature in grid
    double d2eta_dxdy = 0.0; // Assume no cuvature in grid
    double d2xi_dy2 = 0.0;  // Assume no cuvature in grid
    double d2eta_dy2 = 0.0; // Assume no cuvature in grid

    double du_dxi = 1.0;
    double du_deta = 1.0;
    double d2u_dxi2 = d2udxi2(u, p);
    double d2u_dxideta = d2udxideta(u, p);
    double d2u_deta2 = d2udeta2(u, p);

    retval = dy_deta * dy_deta * (
          dxi_dy * dxi_dy * d2u_dxi2 
        + 2.0 * dxi_dy * deta_dy * d2u_dxideta 
        + deta_dy * deta_dy * d2u_deta2 
        + d2xi_dy2 * du_dxi 
        + d2eta_dy2 * du_deta
        );

    return retval;
}
size_t REGULARIZATION::p_index(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}
