//
// programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the HEAT-equation in 2 dimensions, fully implicit with delta-formuation and Modified Newton iteration 
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
size_t nx, size_t ny)
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

    for (size_t i = 1; i < nx - 1; ++i)
    {
        for (size_t j = 1; j < ny - 1; ++j)
        {
            size_t ph_sw = bmp_idx(i - 1, j - 1, ny);  
            size_t ph_s  = bmp_idx(i    , j - 1, ny);  
            size_t ph_se = bmp_idx(i + 1, j - 1, ny);  
            size_t ph_w  = bmp_idx(i - 1, j    , ny);  
            size_t ph_0  = bmp_idx(i    , j    , ny); // central point of control volume
            size_t ph_e  = bmp_idx(i + 1, j    , ny);  
            size_t ph_nw = bmp_idx(i - 1, j + 1, ny);  
            size_t ph_n  = bmp_idx(i    , j + 1, ny);  
            size_t ph_ne = bmp_idx(i + 1, j + 1, ny);  

            size_t c_eq = ph_0;
            
            triplets.emplace_back(c_eq, ph_0     , 1.0);
            triplets.emplace_back(c_eq, ph_sw    , 1.0);
            triplets.emplace_back(c_eq, ph_s     , 1.0);
            triplets.emplace_back(c_eq, ph_se    , 1.0);
            triplets.emplace_back(c_eq, ph_w     , 1.0);
            triplets.emplace_back(c_eq, ph_e     , 1.0);
            triplets.emplace_back(c_eq, ph_nw    , 1.0);
            triplets.emplace_back(c_eq, ph_n     , 1.0);
            triplets.emplace_back(c_eq, ph_ne    , 1.0);
        }
    }
    //
    // Boundary
    //
    if (true)  // nnorth boundary
    {
        size_t j = ny - 1;
        for (size_t i = 1; i < nx - 1; ++i)
        {
            size_t p0 = bmp_idx(i - 1, j    , ny);
            size_t p1 = bmp_idx(i    , j    , ny);
            size_t p2 = bmp_idx(i + 1, j    , ny);
            size_t p3 = bmp_idx(i - 1, j - 1, ny);
            size_t p4 = bmp_idx(i    , j - 1, ny);
            size_t p5 = bmp_idx(i + 1, j - 1, ny);
            size_t p6 = bmp_idx(i - 1, j - 2, ny);
            size_t p7 = bmp_idx(i    , j - 2, ny);
            size_t p8 = bmp_idx(i + 1, j - 2, ny);

            size_t c_eq = p1;    // continuity equation

            triplets.emplace_back(c_eq, p0    , 1.0);
            triplets.emplace_back(c_eq, p1    , 1.0);
            triplets.emplace_back(c_eq, p2    , 1.0);
            triplets.emplace_back(c_eq, p3    , 1.0);
            triplets.emplace_back(c_eq, p4    , 1.0);
            triplets.emplace_back(c_eq, p5    , 1.0);
            triplets.emplace_back(c_eq, p6    , 1.0);
            triplets.emplace_back(c_eq, p7    , 1.0);
            triplets.emplace_back(c_eq, p8    , 1.0);

        }
    }
    if (true)  // eeast boundary
    {
        size_t i = nx - 1;
        for (size_t j = 1; j < ny - 1; ++j)
        {
            size_t p0 = bmp_idx(i    , j - 1, ny);
            size_t p1 = bmp_idx(i - 1, j - 1, ny);
            size_t p2 = bmp_idx(i - 2, j - 1, ny);
            size_t p3 = bmp_idx(i    , j    , ny);
            size_t p4 = bmp_idx(i - 1, j    , ny);
            size_t p5 = bmp_idx(i - 2, j    , ny);
            size_t p6 = bmp_idx(i    , j + 1, ny);
            size_t p7 = bmp_idx(i - 1, j + 1, ny);
            size_t p8 = bmp_idx(i - 2, j + 1, ny);

            size_t c_eq = p3;    // continuity equation

            triplets.emplace_back(c_eq, p0    , 1.0);
            triplets.emplace_back(c_eq, p1    , 1.0);
            triplets.emplace_back(c_eq, p2    , 1.0);
            triplets.emplace_back(c_eq, p3    , 1.0);
            triplets.emplace_back(c_eq, p4    , 1.0);
            triplets.emplace_back(c_eq, p5    , 1.0);
            triplets.emplace_back(c_eq, p6    , 1.0);
            triplets.emplace_back(c_eq, p7    , 1.0);
            triplets.emplace_back(c_eq, p8    , 1.0);

        }
    }
    if (true)  // ssouth boundary
    {
        size_t j = 0;
        for (size_t i = 1; i < nx - 1; ++i)
        {
            size_t p0 = bmp_idx(i - 1, j    , ny);
            size_t p1 = bmp_idx(i    , j    , ny);
            size_t p2 = bmp_idx(i + 1, j    , ny);
            size_t p3 = bmp_idx(i - 1, j + 1, ny);
            size_t p4 = bmp_idx(i    , j + 1, ny);
            size_t p5 = bmp_idx(i + 1, j + 1, ny);
            size_t p6 = bmp_idx(i - 1, j + 2, ny);
            size_t p7 = bmp_idx(i    , j + 2, ny);
            size_t p8 = bmp_idx(i + 1, j + 2, ny);

            size_t c_eq = p1;    // continuity equation

            triplets.emplace_back(c_eq, p0    , 1.0);
            triplets.emplace_back(c_eq, p1    , 1.0);
            triplets.emplace_back(c_eq, p2    , 1.0);
            triplets.emplace_back(c_eq, p3    , 1.0);
            triplets.emplace_back(c_eq, p4    , 1.0);
            triplets.emplace_back(c_eq, p5    , 1.0);
            triplets.emplace_back(c_eq, p6    , 1.0);
            triplets.emplace_back(c_eq, p7    , 1.0);
            triplets.emplace_back(c_eq, p8    , 1.0);

        }
    }
    if (true)  // wwest boundary(2D)
    {
        size_t i = 0;
        for (size_t j = 1; j < ny - 1; ++j)
        {            
            size_t p0 = bmp_idx(i    , j - 1, ny);
            size_t p1 = bmp_idx(i + 1, j - 1, ny);
            size_t p2 = bmp_idx(i + 2, j - 1, ny);
            size_t p3 = bmp_idx(i    , j    , ny);
            size_t p4 = bmp_idx(i + 1, j    , ny);
            size_t p5 = bmp_idx(i + 2, j    , ny);
            size_t p6 = bmp_idx(i    , j + 1, ny);
            size_t p7 = bmp_idx(i + 1, j + 1, ny);
            size_t p8 = bmp_idx(i + 2, j + 1, ny);

            size_t c_eq = p3;     // continuity equation

            triplets.emplace_back(c_eq, p0    , 1.0);
            triplets.emplace_back(c_eq, p1    , 1.0);
            triplets.emplace_back(c_eq, p2    , 1.0);
            triplets.emplace_back(c_eq, p3    , 1.0);
            triplets.emplace_back(c_eq, p4    , 1.0);
            triplets.emplace_back(c_eq, p5    , 1.0);
            triplets.emplace_back(c_eq, p6    , 1.0);
            triplets.emplace_back(c_eq, p7    , 1.0);
            triplets.emplace_back(c_eq, p8    , 1.0);

        }
    }
    //
    //corner nodes
    //
    if (true)  // NE-corner
    {
        size_t i = nx - 1;
        size_t j = ny - 1;

        size_t p0 = bmp_idx(i    , j    , ny);
        size_t p1 = bmp_idx(i - 1, j    , ny);
        size_t p2 = bmp_idx(i - 2, j    , ny);
        size_t p3 = bmp_idx(i    , j - 1, ny);
        size_t p4 = bmp_idx(i - 1, j - 1, ny);
        size_t p5 = bmp_idx(i - 2, j - 1, ny);
        size_t p6 = bmp_idx(i    , j - 2, ny);
        size_t p7 = bmp_idx(i - 1, j - 2, ny);
        size_t p8 = bmp_idx(i - 2, j - 2, ny);

        size_t c_eq = p0;    // continuity equation

        triplets.emplace_back(c_eq, p0    , 1.0);
        triplets.emplace_back(c_eq, p1    , 1.0);
        triplets.emplace_back(c_eq, p2    , 1.0);
        triplets.emplace_back(c_eq, p3    , 1.0);
        triplets.emplace_back(c_eq, p4    , 1.0);
        triplets.emplace_back(c_eq, p5    , 1.0);
        triplets.emplace_back(c_eq, p6    , 1.0);
        triplets.emplace_back(c_eq, p7    , 1.0);
        triplets.emplace_back(c_eq, p8    , 1.0);

    }
    if (true)  // SE-corner
    {
        size_t i = nx - 1;
        size_t j = 0;

        size_t p0 = bmp_idx(i    , j   , ny);
        size_t p1 = bmp_idx(i - 1, j   , ny);
        size_t p2 = bmp_idx(i - 2, j   , ny);
        size_t p3 = bmp_idx(i    , j + 1, ny);
        size_t p4 = bmp_idx(i - 1, j + 1, ny);
        size_t p5 = bmp_idx(i - 2, j + 1, ny);
        size_t p6 = bmp_idx(i    , j + 2, ny);
        size_t p7 = bmp_idx(i - 1, j + 2, ny);
        size_t p8 = bmp_idx(i - 2, j + 2, ny);

        size_t c_eq = p0;    // continuity equation

        triplets.emplace_back(c_eq, p0    , 1.0);
        triplets.emplace_back(c_eq, p1    , 1.0);
        triplets.emplace_back(c_eq, p2    , 1.0);
        triplets.emplace_back(c_eq, p3    , 1.0);
        triplets.emplace_back(c_eq, p4    , 1.0);
        triplets.emplace_back(c_eq, p5    , 1.0);
        triplets.emplace_back(c_eq, p6    , 1.0);
        triplets.emplace_back(c_eq, p7    , 1.0);
        triplets.emplace_back(c_eq, p8    , 1.0);
    }
    if (true)  // SW-corner
    {
        size_t i = 0;
        size_t j = 0;

        size_t p0 = bmp_idx(i    , j   , ny);
        size_t p1 = bmp_idx(i + 1, j   , ny);
        size_t p2 = bmp_idx(i + 2, j   , ny);
        size_t p3 = bmp_idx(i    , j + 1, ny);
        size_t p4 = bmp_idx(i + 1, j + 1, ny);
        size_t p5 = bmp_idx(i + 2, j + 1, ny);
        size_t p6 = bmp_idx(i    , j + 2, ny);
        size_t p7 = bmp_idx(i + 1, j + 2, ny);
        size_t p8 = bmp_idx(i + 2, j + 2, ny);

        size_t c_eq = p0;    // continuity equation

        triplets.emplace_back(c_eq, p0    , 1.0);
        triplets.emplace_back(c_eq, p1    , 1.0);
        triplets.emplace_back(c_eq, p2    , 1.0);
        triplets.emplace_back(c_eq, p3    , 1.0);
        triplets.emplace_back(c_eq, p4    , 1.0);
        triplets.emplace_back(c_eq, p5    , 1.0);
        triplets.emplace_back(c_eq, p6    , 1.0);
        triplets.emplace_back(c_eq, p7    , 1.0);
        triplets.emplace_back(c_eq, p8    , 1.0);

    }
    if (true)  // NW-corner
    {
        size_t i = 0;
        size_t j = ny - 1;

        size_t p0 = bmp_idx(i    , j   , ny);
        size_t p1 = bmp_idx(i + 1, j   , ny);
        size_t p2 = bmp_idx(i + 2, j   , ny);
        size_t p3 = bmp_idx(i    , j - 1, ny);
        size_t p4 = bmp_idx(i + 1, j - 1, ny);
        size_t p5 = bmp_idx(i + 2, j - 1, ny);
        size_t p6 = bmp_idx(i    , j - 2, ny);
        size_t p7 = bmp_idx(i + 1, j - 2, ny);
        size_t p8 = bmp_idx(i + 2, j - 2, ny);

        size_t c_eq = p0;    // continuity equation

        triplets.emplace_back(c_eq, p0    , 1.0);
        triplets.emplace_back(c_eq, p1    , 1.0);
        triplets.emplace_back(c_eq, p2    , 1.0);
        triplets.emplace_back(c_eq, p3    , 1.0);
        triplets.emplace_back(c_eq, p4    , 1.0);
        triplets.emplace_back(c_eq, p5    , 1.0);
        triplets.emplace_back(c_eq, p6    , 1.0);
        triplets.emplace_back(c_eq, p7    , 1.0);
        triplets.emplace_back(c_eq, p8    , 1.0);
    }
        
    return 0;
}
inline size_t bmp_idx(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}

