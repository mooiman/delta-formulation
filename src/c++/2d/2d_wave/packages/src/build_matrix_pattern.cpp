//
// programmer: Jan Mooiman
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

#include "build_matrix_pattern.h"

//------------------------------------------------------------------------------
int build_matrix_pattern(std::vector< Eigen::Triplet<double> >& triplets,
int nx, int ny)
{
    //    x - - - - - - - x - - - - - - - x 
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    x - - - - - - - x - - - - - - - x  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    |               |               |  
    //    x - - - - - - - x - - - - - - - x 

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int ph_sw = bmp_index(i - 1, j - 1, ny);  
            int ph_s  = bmp_index(i    , j - 1, ny);  
            int ph_se = bmp_index(i + 1, j - 1, ny);  
            int ph_w  = bmp_index(i - 1, j    , ny);  
            int ph_0  = bmp_index(i    , j    , ny); // central point of control volume
            int ph_e  = bmp_index(i + 1, j    , ny);  
            int ph_nw = bmp_index(i - 1, j + 1, ny);  
            int ph_n  = bmp_index(i    , j + 1, ny);  
            int ph_ne = bmp_index(i + 1, j + 1, ny);  

            int c_eq = 3 * ph_0;
            int q_eq = 3 * ph_0 + 1;
            int r_eq = 3 * ph_0 + 2;
            
            triplets.emplace_back(c_eq, 3 * ph_0     , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_sw    , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_s     , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_se    , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_w     , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_e     , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_nw    , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_n     , 1.0);
            triplets.emplace_back(c_eq, 3 * ph_ne    , 1.0);

            triplets.emplace_back(c_eq, 3 * ph_0  + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_sw + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_s  + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_se + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_w  + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_e  + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_nw + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_n  + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_ne + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * ph_0  + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_sw + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_s  + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_se + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_w  + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_e  + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_nw + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_n  + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * ph_ne + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * ph_0     , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_sw    , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_s     , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_se    , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_w     , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_e     , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_nw    , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_n     , 1.0);
            triplets.emplace_back(q_eq, 3 * ph_ne    , 1.0);

            triplets.emplace_back(q_eq, 3 * ph_0  + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_sw + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_s  + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_se + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_w  + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_e  + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_nw + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_n  + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_ne + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * ph_0  + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_sw + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_s  + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_se + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_w  + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_e  + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_nw + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_n  + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * ph_ne + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * ph_0     , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_sw    , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_s     , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_se    , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_w     , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_e     , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_nw    , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_n     , 1.0);
            triplets.emplace_back(r_eq, 3 * ph_ne    , 1.0);

            triplets.emplace_back(r_eq, 3 * ph_0  + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_sw + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_s  + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_se + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_w  + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_e  + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_nw + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_n  + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_ne + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * ph_0  + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_sw + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_s  + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_se + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_w  + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_e  + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_nw + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_n  + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * ph_ne + 2, 1.0);
        }
    }
    //
    // Boundary
    //
    if (true)  // nnorth boundary
    {
        int j = ny - 1;
        for (int i = 1; i < nx - 1; ++i)
        {
            int p0 = bmp_index(i - 1, j    , ny);
            int p1 = bmp_index(i    , j    , ny);
            int p2 = bmp_index(i + 1, j    , ny);
            int p3 = bmp_index(i - 1, j - 1, ny);
            int p4 = bmp_index(i    , j - 1, ny);
            int p5 = bmp_index(i + 1, j - 1, ny);
            int p6 = bmp_index(i - 1, j - 2, ny);
            int p7 = bmp_index(i    , j - 2, ny);
            int p8 = bmp_index(i + 1, j - 2, ny);

            int c_eq = 3 * p1;    // continuity equation
            int q_eq = c_eq + 1;  // q-momentum equation
            int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
        }
    }
    if (true)  // eeast boundary
    {
        int i = nx - 1;
        for (int j = 1; j < ny - 1; ++j)
        {
            int p0 = bmp_index(i    , j - 1, ny);
            int p1 = bmp_index(i - 1, j - 1, ny);
            int p2 = bmp_index(i - 2, j - 1, ny);
            int p3 = bmp_index(i    , j    , ny);
            int p4 = bmp_index(i - 1, j    , ny);
            int p5 = bmp_index(i - 2, j    , ny);
            int p6 = bmp_index(i    , j + 1, ny);
            int p7 = bmp_index(i - 1, j + 1, ny);
            int p8 = bmp_index(i - 2, j + 1, ny);

            int c_eq = 3 * p3;    // continuity equation
            int q_eq = c_eq + 1;  // q-momentum equation
            int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
        }
    }
    if (true)  // ssouth boundary
    {
        int j = 0;
        for (int i = 1; i < nx - 1; ++i)
        {
            int p0 = bmp_index(i - 1, j    , ny);
            int p1 = bmp_index(i    , j    , ny);
            int p2 = bmp_index(i + 1, j    , ny);
            int p3 = bmp_index(i - 1, j + 1, ny);
            int p4 = bmp_index(i    , j + 1, ny);
            int p5 = bmp_index(i + 1, j + 1, ny);
            int p6 = bmp_index(i - 1, j + 2, ny);
            int p7 = bmp_index(i    , j + 2, ny);
            int p8 = bmp_index(i + 1, j + 2, ny);

            int c_eq = 3 * p1;    // continuity equation
            int q_eq = c_eq + 1;  // q-momentum equation
            int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
        }
    }
    if (true)  // wwest boundary(2D)
    {
        int i = 0;
        for (int j = 1; j < ny - 1; ++j)
        {            
            int p0 = bmp_index(i    , j - 1, ny);
            int p1 = bmp_index(i + 1, j - 1, ny);
            int p2 = bmp_index(i + 2, j - 1, ny);
            int p3 = bmp_index(i    , j    , ny);
            int p4 = bmp_index(i + 1, j    , ny);
            int p5 = bmp_index(i + 2, j    , ny);
            int p6 = bmp_index(i    , j + 1, ny);
            int p7 = bmp_index(i + 1, j + 1, ny);
            int p8 = bmp_index(i + 2, j + 1, ny);

            int c_eq = 3 * p3;     // continuity equation
            int q_eq = c_eq + 1;   // q-momentum equation
            int r_eq = c_eq + 2;   // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
        }
    }
    //
    //corner nodes
    //
    if (true)  // NE-corner
    {
        int i = nx - 1;
        int j = ny - 1;

        int p0 = bmp_index(i    , j    , ny);
        int p1 = bmp_index(i - 1, j    , ny);
        int p2 = bmp_index(i - 2, j    , ny);
        int p3 = bmp_index(i    , j - 1, ny);
        int p4 = bmp_index(i - 1, j - 1, ny);
        int p5 = bmp_index(i - 2, j - 1, ny);
        int p6 = bmp_index(i    , j - 2, ny);
        int p7 = bmp_index(i - 1, j - 2, ny);
        int p8 = bmp_index(i - 2, j - 2, ny);

        int c_eq = 3 * p0;    // continuity equation
        int q_eq = c_eq + 1;  // q-momentum equation
        int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
    }
    if (true)  // SE-corner
    {
        int i = nx - 1;
        int j = 0;

        int p0 = bmp_index(i    , j   , ny);
        int p1 = bmp_index(i - 1, j   , ny);
        int p2 = bmp_index(i - 2, j   , ny);
        int p3 = bmp_index(i    , j + 1, ny);
        int p4 = bmp_index(i - 1, j + 1, ny);
        int p5 = bmp_index(i - 2, j + 1, ny);
        int p6 = bmp_index(i    , j + 2, ny);
        int p7 = bmp_index(i - 1, j + 2, ny);
        int p8 = bmp_index(i - 2, j + 2, ny);

        int c_eq = 3 * p0;    // continuity equation
        int q_eq = c_eq + 1;  // q-momentum equation
        int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
    }
    if (true)  // SW-corner
    {
        int i = 0;
        int j = 0;

        int p0 = bmp_index(i    , j   , ny);
        int p1 = bmp_index(i + 1, j   , ny);
        int p2 = bmp_index(i + 2, j   , ny);
        int p3 = bmp_index(i    , j + 1, ny);
        int p4 = bmp_index(i + 1, j + 1, ny);
        int p5 = bmp_index(i + 2, j + 1, ny);
        int p6 = bmp_index(i    , j + 2, ny);
        int p7 = bmp_index(i + 1, j + 2, ny);
        int p8 = bmp_index(i + 2, j + 2, ny);

        int c_eq = 3 * p0;    // continuity equation
        int q_eq = c_eq + 1;  // q-momentum equation
        int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
    }
    if (true)  // NW-corner
    {
        int i = 0;
        int j = ny - 1;
        int p0 = bmp_index(i    , j   , ny);
        int p1 = bmp_index(i + 1, j   , ny);
        int p2 = bmp_index(i + 2, j   , ny);
        int p3 = bmp_index(i    , j - 1, ny);
        int p4 = bmp_index(i + 1, j - 1, ny);
        int p5 = bmp_index(i + 2, j - 1, ny);
        int p6 = bmp_index(i    , j - 2, ny);
        int p7 = bmp_index(i + 1, j - 2, ny);
        int p8 = bmp_index(i + 2, j - 2, ny);

        int c_eq = 3 * p0;    // continuity equation
        int q_eq = c_eq + 1;  // q-momentum equation
        int r_eq = c_eq + 2;  // r-momentum equation

            triplets.emplace_back(c_eq, 3 * p0    , 1.0);
            triplets.emplace_back(c_eq, 3 * p1    , 1.0);
            triplets.emplace_back(c_eq, 3 * p2    , 1.0);
            triplets.emplace_back(c_eq, 3 * p3    , 1.0);
            triplets.emplace_back(c_eq, 3 * p4    , 1.0);
            triplets.emplace_back(c_eq, 3 * p5    , 1.0);
            triplets.emplace_back(c_eq, 3 * p6    , 1.0);
            triplets.emplace_back(c_eq, 3 * p7    , 1.0);
            triplets.emplace_back(c_eq, 3 * p8    , 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(c_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(c_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(q_eq, 3 * p0    , 1.0);
            triplets.emplace_back(q_eq, 3 * p1    , 1.0);
            triplets.emplace_back(q_eq, 3 * p2    , 1.0);
            triplets.emplace_back(q_eq, 3 * p3    , 1.0);
            triplets.emplace_back(q_eq, 3 * p4    , 1.0);
            triplets.emplace_back(q_eq, 3 * p5    , 1.0);
            triplets.emplace_back(q_eq, 3 * p6    , 1.0);
            triplets.emplace_back(q_eq, 3 * p7    , 1.0);
            triplets.emplace_back(q_eq, 3 * p8    , 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(q_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(q_eq, 3 * p8 + 2, 1.0);

            triplets.emplace_back(r_eq, 3 * p0    , 1.0);
            triplets.emplace_back(r_eq, 3 * p1    , 1.0);
            triplets.emplace_back(r_eq, 3 * p2    , 1.0);
            triplets.emplace_back(r_eq, 3 * p3    , 1.0);
            triplets.emplace_back(r_eq, 3 * p4    , 1.0);
            triplets.emplace_back(r_eq, 3 * p5    , 1.0);
            triplets.emplace_back(r_eq, 3 * p6    , 1.0);
            triplets.emplace_back(r_eq, 3 * p7    , 1.0);
            triplets.emplace_back(r_eq, 3 * p8    , 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 1, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 1, 1.0);

            triplets.emplace_back(r_eq, 3 * p0 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p1 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p2 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p3 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p4 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p5 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p6 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p7 + 2, 1.0);
            triplets.emplace_back(r_eq, 3 * p8 + 2, 1.0);
    }
        
    return 0;
}
inline int bmp_index(int i, int j, int ny_in)
{
    return i * ny_in + j;
}

