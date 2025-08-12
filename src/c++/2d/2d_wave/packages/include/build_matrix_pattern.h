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

#ifndef __BUILD_MATRIX_PATTERN_H__
#define __BUILD_MATRIX_PATTERN_H__

#include <vector>
#include <Eigen/Sparse>

    int build_matrix_pattern(std::vector< Eigen::Triplet<double> >& triplets, int nx, int ny);
    inline int bmp_index(int i, int j, int nx);

#endif  // __BUILD_MATRIX_PATTERN_H__
