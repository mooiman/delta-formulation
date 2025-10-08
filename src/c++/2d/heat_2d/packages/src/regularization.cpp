//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formulation and Modified Newton iteration 
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
#include "build_matrix_pattern_regularization.h"

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
    std::vector<double>& u_out, std::vector<double>& psi_11, std::vector<double>& psi_22, std::vector<double>& eq8, 
    std::vector<double>& u_giv_in,
    size_t nx, size_t ny, double dx, double dy, double c_psi, std::ofstream& log_file)
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
        eq8 = *(this->solve_eq8(nx, ny, dx, dy, c_psi, u0, u0_xixi, u0_etaeta, log_file));
//------------------------------------------------------------------------------
        for (size_t i = 0; i < nxny; ++i)
        {
            psi_11[i] = c_psi * (dx * dx + dy * dy) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
            psi_22[i] = c_psi * (dx * dx + dy * dy) * eq8[i]/2.0;  // divide by 2: then is equal to 1D if dx=dy
        }
//------------------------------------------------------------------------------
        u0 = *(this->solve_eq7(nx, ny, dx, dy, psi_11, psi_22, u_giv, log_file));
//------------------------------------------------------------------------------

        diff_max1 = 0.0;
        for (size_t i = 0; i < nxny; ++i)
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
    for (size_t i = 0; i < nxny; ++i)
    {
        u_out[i] = u0[i] * u_giv_range;
    }
}

void REGULARIZATION::first_derivative(std::vector<double>& psi, std::vector<double>& eps, std::vector<double>& u, double dx)
{
    size_t nx = (int)eps.size();
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
std::unique_ptr<std::vector<double>> REGULARIZATION::solve_eq7(size_t nx, size_t ny, double dx, double dy, 
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

    size_t col;
    size_t p_0;

    // south-west corner
    for (size_t row = 0; row < 1; row += 1)
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_sw

        values[col     ] = 1.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = u_giv[p_0];
    }
    // west boundary
    for (size_t row = 1; row < (ny - 1); row += 1)
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_w

        values[col     ] = 0.0;
        values[col +  1] = -1.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-west corner
    for (size_t row = ny - 1; row < ny; row += 1)
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_nw

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 1.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = u_giv[p_0];
    }
    //
    // first internal west boundary column
    for (size_t row = ny; row < 2 * ny; row += 1) 
    {
        if (std::fmod(row , ny) == 0) {
            // south boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_s

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 1.0;
            values[col +  4] = -1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if (std::fmod(row + 1, ny) == 0) {
            // north boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_n

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = -1.0;
            values[col +  5] = 1.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        col = outer[row    ];
        p_0 = col/(9);  // p_nw

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = u_giv[p_0];
    }

    // interior with south and north boundary
    for (size_t row = 2 * ny; row < (nx - 2) * ny; row += 1) 
    {
        if (std::fmod(row , ny) == 0) {
            // south boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_s

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 1.0;
            values[col +  4] = -1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if (std::fmod(row - 1, ny) == 0) {
            // south boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_s

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = 1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = u_giv[p_0];
            continue;
        }

        if (std::fmod(row + 1, ny) == 0) {
            // north boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_n

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = -1.0;
            values[col +  5] = 1.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if (std::fmod(row + 2, ny) == 0) {
            // north boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_n

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = 1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = u_giv[p_0];
            continue;
        }

        col = outer[row    ];
        size_t p0 = col/(9);  // p_s

        size_t p_sw = p0 - ny - 1;
        size_t p_w  = p0 - ny    ;
        size_t p_nw = p0 - ny + 1;
        size_t p_s  = p0 - 1;
        size_t p_0  = p0    ;
        size_t p_n  = p0 + 1;
        size_t p_se = p0 + ny - 1;
        size_t p_e  = p0 + ny    ;
        size_t p_ne = p0 + ny + 1;
        
        values[col     ] = dx * dy * m_mass[0] * m_mass[0];
        values[col +  1] = dx * dy * m_mass[0] * m_mass[1];
        values[col +  2] = dx * dy * m_mass[0] * m_mass[2];
        values[col +  3] = dx * dy * m_mass[1] * m_mass[0];
        values[col +  4] = dx * dy * m_mass[1] * m_mass[1];
        values[col +  5] = dx * dy * m_mass[1] * m_mass[2];
        values[col +  6] = dx * dy * m_mass[2] * m_mass[0];
        values[col +  7] = dx * dy * m_mass[2] * m_mass[1];
        values[col +  8] = dx * dy * m_mass[2] * m_mass[2];

        // psi should be computed on the 8 interfaces of the control volume
        double psi_n = 0.5 * (psi_22[p_n] + psi_22[p_0]);
        double psi_e = 0.5 * (psi_11[p_e] + psi_11[p_0]);
        double psi_s = 0.5 * (psi_22[p_s] + psi_22[p_0]);
        double psi_w = 0.5 * (psi_11[p_w] + psi_11[p_0]);
        //psi_im12 = 2. * psi[i] * psi[i - 1] / (psi[i] + psi[i - 1]);
        //psi_ip12 = 2. * psi[i + 1] * psi[i] / (psi[i + 1] + psi[i]);

        values[col    ] += 0.0;
        values[col + 1] += -psi_s * dx / dy;
        values[col + 2] += 0.0;
        values[col + 3] += -psi_w * dy / dx;
        values[col + 4] +=  psi_s * dx / dy + psi_w * dy / dx + psi_n * dy / dx + psi_e * dx / dy;
        values[col + 5] += -psi_e * dy / dx;
        values[col + 6] += 0.0;
        values[col + 7] += -psi_n * dx / dy;
        values[col + 8] += 0.0;
            
        rhs[row]  = dx * dy * m_mass[0] * m_mass[0] * u_giv[p_sw];
        rhs[row] += dx * dy * m_mass[0] * m_mass[1] * u_giv[p_s ];
        rhs[row] += dx * dy * m_mass[0] * m_mass[2] * u_giv[p_se];

        rhs[row] += dx * dy * m_mass[1] * m_mass[0] * u_giv[p_w ];
        rhs[row] += dx * dy * m_mass[1] * m_mass[1] * u_giv[p_0 ];
        rhs[row] += dx * dy * m_mass[1] * m_mass[2] * u_giv[p_e ];

        rhs[row] += dx * dy * m_mass[2] * m_mass[0] * u_giv[p_nw];
        rhs[row] += dx * dy * m_mass[2] * m_mass[1] * u_giv[p_n ];
        rhs[row] += dx * dy * m_mass[2] * m_mass[2] * u_giv[p_ne];
    }
    // first internal east boundary column
    for (size_t row = (nx - 2) * ny; row < (nx - 1) * ny; row += 1) 
    {
        if (std::fmod(row , ny) == 0) {
            // south boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_s

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 1.0;
            values[col +  4] = -1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if (std::fmod(row + 1, ny) == 0) {
            // north boundary
            col = outer[row    ];
            p_0 = col/(9);  // p_n

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = -1.0;
            values[col +  5] = 1.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        col = outer[row    ];
        p_0 = col/(9);  // p_nw

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = u_giv[p_0];
    }
    // south-east corner
    for (size_t row = (nx - 1) * ny; row < (nx - 1) * ny + 1; row += 1)
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_se

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = 1.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = u_giv[p_0];
}
    // east boundary
    for (size_t row = (nx - 1) * ny + 1; row < nx * ny - 1; row += 1) 
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_e

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = -1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 1.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-east corner
    for (size_t row = nx * ny - 1; row < nx * ny; row += 1)
    {
        col = outer[row    ];
        p_0 = col/(9);  // p_ne

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 1.0;
        rhs[row    ] = u_giv[p_0];
    }

    if (false)
    {
        log_file << "=== Matrix eq7 ========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            for (size_t j = 0; j < nxny; ++j)
            {
                log_file << std::showpos << std::setprecision(3) << std::scientific << B.coeff(i, j) << " ";
                if (std::fmod(j+1,ny) == 0) { log_file << "| "; }
            }
            log_file << std::endl;
            if (std::fmod(i+1,ny) == 0) { log_file << std::endl; }
        }
        log_file << "=== RHS eq7 ===========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
            if (std::fmod(i+1,ny) == 0) { log_file << std::endl; }
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

    if (false) 
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
            if (std::fmod(i+1,ny) == 0) { log_file << std::endl; }
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

    double* values = A.valuePtr();         // pointer to all non-zero values
    const int* outer = A.outerIndexPtr();   // row start pointers

    size_t col;
    size_t p0;

    // south-west corner
    for (size_t row = 0; row < 1; row += 1)
    {
        col = outer[row    ];
        p0 = col/(9);  // p_sw

        values[col     ] = -2.0;
        values[col +  1] = 1.0;
        values[col +  2] = 0.0;
        values[col +  3] = 1.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // west boundary
    for (size_t row = 1; row < (ny - 1); row += 1)
    {
        col = outer[row    ];
        p0 = col/(9);  // p_w

        values[col     ] = 0.0;
        values[col +  1] = -1.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-west corner
    for (size_t row = ny - 1; row < ny; row += 1)
    {
        col = outer[row    ];
        p0 = col/(9);  // p_nw

        values[col     ] = 0.0;
        values[col +  1] = 1.0;
        values[col +  2] = -2.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 1.0;
        values[col +  6] = 0.0;
        values[col +  7] = 0.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // interior with south and north boundary
    for (size_t row = ny; row < (nx - 1) * ny; row += 1) 
    {
        col = outer[row    ];
        p0 = col/(9);  // p_0

        if (std::fmod(row , ny) == 0) {
            // south boundary
            col = outer[row    ];
            p0 = col/(9);  // p_s

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = -1.0;
            values[col +  4] = 1.0;
            values[col +  5] = 0.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
            rhs[row    ] = 0.0;
            continue;
        }
        if (std::fmod(row + 1, ny) == 0) {
            // north boundary
            col = outer[row    ];
            p0 = col/(9);  // p_n

            values[col     ] = 0.0;
            values[col +  1] = 0.0;
            values[col +  2] = 0.0;
            values[col +  3] = 0.0;
            values[col +  4] = 1.0;
            values[col +  5] = -1.0;
            values[col +  6] = 0.0;
            values[col +  7] = 0.0;
            values[col +  8] = 0.0;
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

        //values[col     ] = m_mass[0] * m_mass[0];
        //values[col +  1] = m_mass[1] * m_mass[0] - c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[col +  2] = m_mass[0] * m_mass[2];
        //values[col +  3] = m_mass[0] * m_mass[1] - c_error * u0_xixi[p_w] / u0_xixi_max;
        //values[col +  4] = m_mass[1] * m_mass[1] + 2. * c_error * u0_xixi[p_w] / u0_xixi_max + 2. * c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[col +  5] = m_mass[2] * m_mass[1] - c_error * u0_xixi[p_w] / u0_xixi_max;
        //values[col +  6] = m_mass[2] * m_mass[0];
        //values[col +  7] = m_mass[1] * m_mass[2] - c_error * u0_etaeta[p_w] / u0_etaeta_max;
        //values[col +  8] = m_mass[2] * m_mass[2];
        
        values[col     ] = m_mass[0] * m_mass[0];
        values[col +  1] = m_mass[0] * m_mass[1] - c_error;
        values[col +  2] = m_mass[0] * m_mass[2];
        values[col +  3] = m_mass[1] * m_mass[0] - c_error;
        values[col +  4] = m_mass[1] * m_mass[1] + 4. * c_error;
        values[col +  5] = m_mass[1] * m_mass[2] - c_error;
        values[col +  6] = m_mass[2] * m_mass[0];
        values[col +  7] = m_mass[2] * m_mass[1] - c_error;
        values[col +  8] = m_mass[2] * m_mass[2];

        rhs[row]  = (u0_xixi[p_sw] + u0_etaeta[p_sw]);
        rhs[row] += (u0_xixi[p_w ] + u0_etaeta[p_w ]);
        rhs[row] += (u0_xixi[p_nw] + u0_etaeta[p_nw]);
        
        rhs[row] += (u0_xixi[p_s ] + u0_etaeta[p_s ]);
        rhs[row] += (u0_xixi[p_0 ] + u0_etaeta[p_0 ]);
        rhs[row] += (u0_xixi[p_n ] + u0_etaeta[p_n ]);
        
        rhs[row] += (u0_xixi[p_se] + u0_etaeta[p_se]);
        rhs[row] += (u0_xixi[p_e ] + u0_etaeta[p_e ]);
        rhs[row] += (u0_xixi[p_ne] + u0_etaeta[p_ne]);
        rhs[row] = std::abs(rhs[p_0]);
    }
    // south-east corner
    for (size_t row = (nx - 1) * ny; row < (nx - 1) * ny + 1; row += 1)
    {
        col = outer[row    ];
        p0 = col/(9);  // p_se

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 1.0;
        values[col +  4] = 0.0;
        values[col +  5] = 0.0;
        values[col +  6] = -2.0;
        values[col +  7] = 1.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
}
    // east boundary
    for (size_t row = (nx - 1) * ny + 1; row < nx * ny - 1; row += 1) 
    {
        col = outer[row    ];
        p0 = col/(9);  // p_e

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = -1.0;
        values[col +  5] = 0.0;
        values[col +  6] = 0.0;
        values[col +  7] = 1.0;
        values[col +  8] = 0.0;
        rhs[row    ] = 0.0;
    }
    // north-east corner
    for (size_t row = nx * ny - 1; row < nx * ny; row += 1)
    {
        col = outer[row    ];
        p0 = col/(9);  // p_ne

        values[col     ] = 0.0;
        values[col +  1] = 0.0;
        values[col +  2] = 0.0;
        values[col +  3] = 0.0;
        values[col +  4] = 0.0;
        values[col +  5] = 1.0;
        values[col +  6] = 0.0;
        values[col +  7] = 1.0;
        values[col +  8] = -2.0;
        rhs[row    ] = 0.0;
    }

    for (size_t i = 0; i < nxny; ++i)
    {
        solution[i] = u0[i];
    }

    if (false)
    {
        log_file << "=== Matrix eq8 ========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            for (size_t j = 0; j < nxny; ++j)
            {
                log_file << std::showpos << std::setprecision(3) << std::scientific << A.coeff(i, j) << " ";
                if (std::fmod(j+1,ny) == 0) { log_file << "| "; }
            }
            log_file << std::endl;
            if (std::fmod(i+1,ny) == 0) { log_file << std::endl; }
        }
        log_file << "=== RHS eq8 ===========================================" << std::endl;
        for (size_t i = 0; i < nxny; ++i)
        {
            log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
            if (std::fmod(i+1,ny) == 0) { log_file << std::endl; }
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
size_t REGULARIZATION::p_index(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}
