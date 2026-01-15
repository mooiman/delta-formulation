//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 2D shallow water equations, fully implicit with delta-formuation and Modified Newton iteration 
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

#include "print_matrix.h"

void print_matrix_pattern(Eigen::SparseMatrix<double, Eigen::RowMajor> A, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file)
{
    size_t nxny = nx * ny;

    log_file << header_text << std::endl;
    for (size_t i = 0; i < n_eq * nxny; ++i)
    {
        for (size_t j = 0; j < n_eq * nxny; ++j)
        {
            if (A.coeff(i, j) != 0.0)
            {
                log_file << "* ";
            }
            else
            {
                log_file << "- ";
            }
            if ((j+1) % n_eq == 0) { if (n_eq != 1) { log_file << "| ";} }
        }
        log_file << std::endl;
        if ((i+1) % n_eq == 0) 
        { 
            if (n_eq != 1) {log_file << std::endl; } 
        }
    }
}
void print_matrix(Eigen::SparseMatrix<double, Eigen::RowMajor> A, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file)
{
    size_t nxny = nx * ny;

    log_file << header_text << std::endl;
    for (size_t i = 0; i < n_eq * nxny; ++i)
    {
        log_file << std::showpos << std::setprecision(3) << i << "   ";
        for (size_t j = 0; j < n_eq * nxny; ++j)
        {
            log_file << std::showpos << std::setprecision(3) << std::scientific << A.coeff(i, j) << " ";
            if ((j+1) % n_eq == 0) { 
                if (n_eq != 1) { log_file << "| "; }
            }
        }
        log_file << std::endl;
        if ((i+1) % n_eq == 0) 
        { 
            if (n_eq != 1) {log_file << std::endl; } 
        }
    }
    //log_file << "=== Diagonal dominant == diag, off_diag, -, + =========" << std::endl;
    //for (size_t i = 0; i < 3 * nxny; ++i)
    //{
    //    double off_diag = 0.0;
    //    double diag = 0.0;
    //    log_file << std::showpos << std::setprecision(3) << i << "   ";
    //    for (int j = 0; j < 3*nxny; ++j)
    //    {
    //        if (i != j ) { off_diag += std::abs(A.coeff(i,j)); }
    //        else { diag = std::abs(A.coeff(i,j)); }
    //    }
    //    log_file << std::showpos << std::setprecision(5) << std::scientific << diag << " " << off_diag << " " 
    //        <<  diag - off_diag << " " <<  diag + off_diag << " ";
    //    log_file << std::endl;
    //    if (std::fmod(i+1,3) == 0) { log_file << std::endl; }
    //}
}
void print_vector(Eigen::VectorXd& rhs, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file)
{
    size_t nxny = nx * ny;

    log_file << header_text << std::endl;
    for (size_t i = 0; i < n_eq * nxny; ++i)
    {
        log_file << std::showpos << std::setprecision(3) << i << "   ";
        log_file << std::setprecision(8) << std::scientific << rhs[i] << std::endl;
        if ((i+1) % n_eq == 0) 
        { 
            if (n_eq != 1) {log_file << std::endl; } 
        }
    }
}
