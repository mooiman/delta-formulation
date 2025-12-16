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

#include "build_matrix_pattern_regularization.h"
#include "perf_timer.h"
#include "regularization.h"

REGULARIZATION::REGULARIZATION()
{
    m_iter_max = 100;
    m_g = 10.0;

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

void REGULARIZATION::given_function( 
    std::vector<double>& u_out, std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& eq8, 
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

        double dxi = 1.0;
        double deta = 1.0;
//------------------------------------------------------------------------------
        eq8 = *(this->solve_eq8(nx, ny, dxi, deta, c_psi, u0, u0_xixi, u0_etaeta, log_file));  // eq8 is computed in computational space
//------------------------------------------------------------------------------
        for (size_t i = 0; i < nxny; ++i)
        {
            psi_11[i] = -c_psi * (dxi * dxi + deta * deta) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
            psi_22[i] = -c_psi * (dxi * dxi + deta * deta) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
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
        u_out[i] = u0[i] * u_giv_range;
    }
}

void REGULARIZATION::first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx)
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
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(size_t nx, size_t ny, std::vector<double>& x, std::vector<double>& y, 
    std::vector<double> psi_11, std::vector<double> psi_22, std::vector<double> u_giv, std::ofstream& log_file)
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
        status = reg_corner_south_west(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }
    // west boundary
    for (size_t row = 1; row < 1 * (ny - 1); row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = psi_11[p_0];
        double psi_2 = psi_22[p_0];
        status = reg_boundary_west(values, row, c_eq, rhs,
            x, y, u_giv, psi_1, psi_2, 
            (double) 1.0, nx, ny);
    }
    // north-west corner
    for (size_t row = 1 * (ny - 1); row < 1 * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_west(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }

    // interior with south and north boundary
    for (size_t row = 1 * ny; row < 1 * (nx - 1) * ny; row += 1) 
    {
        size_t c_eq = outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = psi_11[p_0];
        double psi_2 = psi_22[p_0];

        if (row % ny == 0) {
            // south boundary, over write coefficients
            status = reg_boundary_south(values, row, c_eq, rhs, 
                x, y, u_giv, psi_1, psi_2, 
                (double) 1.0, nx, ny);
            continue;
        }
        if ((row + 1) % ny == 0) {
            // north boundary, over write coefficients
            status = reg_boundary_north(values, row, c_eq, rhs, 
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
        status = reg_corner_south_east(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }
    // east boundary
    for (size_t row = 1 * (nx - 1) * ny + 1; row < 1 * nx * ny - 1; row += 1) 
    {
        int c_eq = (size_t) outer[row    ];
        size_t p_0 = c_eq/(9);

        double psi_1 = psi_11[p_0];
        double psi_2 = psi_22[p_0];

        status = reg_boundary_east(values, row, c_eq, rhs, 
            x, y, u_giv, psi_1, psi_2, 
            (double) 1.0, nx, ny);
    }
    // north-east corner
    for (size_t row = 1 * nx * ny - 1; row < 1 * nx * ny; row += 1)
    {
        size_t c_eq = (size_t) outer[row    ];
        status = reg_corner_north_east(values, row, c_eq, rhs, 
           u_giv, (double) 1.0, nx, ny);
    }
    STOP_TIMER(Regularization diffusion);


    if (m_logging == "matrix")
    {
        log_file << "=== Matrix eq7 ========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            for (size_t j = 0; j < nxny; ++j)
            {
                log_file << std::showpos << std::setprecision(3) << std::scientific << B.coeff(i, j) << " ";
                if ((j+1) % ny == 0) { log_file << "| "; }
            }
            log_file << std::endl;
            if ((i+1) % ny == 0) { log_file << std::endl; }
        }
        log_file << "=== RHS eq7 ===========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
            if ((i+1) % ny == 0) { log_file << std::endl; }
        }
        log_file << "=======================================================" << std::endl;
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
};


std::unique_ptr<std::vector<double>>  REGULARIZATION::solve_eq8(size_t nx, size_t ny, double dx, double dy, double c_error, 
std::vector<double> u0, std::vector<double> u0_xixi, std::vector<double> u0_etaeta, std::ofstream& log_file)
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
        log_file << "=== Matrix build matrix pattern =======================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            for (size_t j = 0; j < nxny; ++j)
            {
                if (A.coeff(i, j) != 0.0)
                {
                    log_file << "* ";
                }
                else
                {
                    log_file << "- ";
                }
            }
            log_file << std::endl;
            if ((i+1) % ny == 0) { log_file << std::endl; }
        }
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
        
        values[c_eq     ] = m_mass[0] * m_mass[0];
        values[c_eq +  1] = m_mass[0] * m_mass[1] - c_error;
        values[c_eq +  2] = m_mass[0] * m_mass[2];
        values[c_eq +  3] = m_mass[1] * m_mass[0] - c_error;
        values[c_eq +  4] = m_mass[1] * m_mass[1] + 4. * c_error;
        values[c_eq +  5] = m_mass[1] * m_mass[2] - c_error;
        values[c_eq +  6] = m_mass[2] * m_mass[0];
        values[c_eq +  7] = m_mass[2] * m_mass[1] - c_error;
        values[c_eq +  8] = m_mass[2] * m_mass[2];

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
        log_file << "=== Matrix eq8 ========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            for (size_t j = 0; j < nxny; ++j)
            {
                log_file << std::showpos << std::setprecision(3) << std::scientific << A.coeff(i, j) << " ";
                if ((j+1) % ny == 0) { log_file << "| "; }
            }
            log_file << std::endl;
            if ((i+1) % ny == 0) { log_file << std::endl; }
        }
        log_file << "=== RHS eq8 ===========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
            if ((i+1) % ny == 0) { log_file << std::endl; }
        }
        log_file << "=======================================================" << std::endl;
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
size_t REGULARIZATION::p_index(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}
