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
#pragma once

#include <fstream>
#include <iomanip>      // std::setprecision
#include <iostream>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

void print_matrix_pattern(Eigen::SparseMatrix<double, Eigen::RowMajor> A, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file);
void print_matrix(Eigen::SparseMatrix<double, Eigen::RowMajor> A, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file);
void print_vector(Eigen::VectorXd& rhs, size_t n_eq, size_t nx, size_t ny, std::string header_text, std::ofstream& log_file);
