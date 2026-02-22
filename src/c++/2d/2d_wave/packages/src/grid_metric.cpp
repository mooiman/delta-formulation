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

#include "grid_metric.h"

//------------------------------------------------------------------------------
int grid_metric(struct _grid_metric & metric)
{    
    //    3 - - - - - - - 2
    //    |       |       | 
    //    |       |       | 
    //    |       |       | 
    //    | - - - - - - - | 
    //    |       |       | 
    //    |       |       | 
    //    |       |       | 
    //    0 - - - - - - - 1 

    //
    // The terms are added to the matrix coefficients and rhs, they already contain contributions from other terms in momentum equation
    // 
    // loop of over elements
    //
    std::vector<double>& x = metric.x;
    std::vector<double>& y = metric.y; 
    size_t nx = metric.nx;
    size_t ny = metric.ny;
    std::vector<double>& dx_dxi    = metric.dx_dxi ; 
    std::vector<double>& dy_dxi    = metric.dy_dxi ; 
    std::vector<double>& dx_deta   = metric.dx_deta; 
    std::vector<double>& dy_deta   = metric.dy_deta;
    std::vector<double>& ddx_dxi2  = metric.ddx_dxi2; 
    std::vector<double>& ddy_dxi2  = metric.ddy_dxi2; 
    std::vector<double>& ddx_deta2 = metric.ddx_deta2; 
    std::vector<double>& ddy_deta2 = metric.ddy_deta2;

    size_t node0, node1, node2, node3;
    for (size_t i = 0; i < nx - 1; ++i)
    {
        for (size_t j = 0; j < ny - 1; ++j)
        {
            node0 = idx (i    , j    , ny);
            node1 = idx (i + 1, j    , ny);
            node2 = idx (i + 1, j + 1, ny);
            node3 = idx (i    , j + 1, ny);
            double dx = (x[node1] - x[node0]);
            double dy = (y[node1] - y[node0]);
            dx_dxi[node0]   += 0.25 * dx;
            dx_dxi[node1]   += 0.25 * dx;
            dy_dxi[node0]   += 0.25 * dy;
            dy_dxi[node1]   += 0.25 * dy;
            ddx_dxi2[node1] += dx;
            ddx_dxi2[node0] -= dx;
            ddy_dxi2[node1] += dy;
            ddy_dxi2[node0] -= dy;

            dx = (x[node2] - x[node3]);
            dy = (y[node2] - y[node3]);
            dx_dxi[node2]   += 0.25 * dx;
            dx_dxi[node3]   += 0.25 * dx;
            dy_dxi[node2]   += 0.25 * dy;
            dy_dxi[node3]   += 0.25 * dy;
            ddx_dxi2[node2] += dx;
            ddx_dxi2[node3] -= dx;
            ddy_dxi2[node2] += dy;
            ddy_dxi2[node3] -= dy;

            dx = (x[node3] - x[node0]);
            dy = (y[node3] - y[node0]);
            dx_deta[node0]   += 0.25 * dx;
            dx_deta[node3]   += 0.25 * dx;
            dy_deta[node0]   += 0.25 * dy;
            dy_deta[node3]   += 0.25 * dy;
            ddx_deta2[node3] += dx;
            ddx_deta2[node0] -= dx;
            ddy_deta2[node3] += dy;
            ddy_deta2[node0] -= dy;

            dx = (x[node2] - x[node1]);
            dy = (y[node2] - y[node1]);
            dx_deta[node2]   += 0.25 * dx;
            dx_deta[node1]   += 0.25 * dx;
            dy_deta[node2]   += 0.25 * dy;
            dy_deta[node1]   += 0.25 * dy;
            ddx_deta2[node2] += dx;
            ddx_deta2[node1] -= dx;
            ddy_deta2[node2] += dy;
            ddy_deta2[node1] -= dy;
        }
    }
    // west boundary
    size_t i = 0;
    size_t j = 0;
    i = 0;  // for (size_t i = 0; i < 1; ++i)
    {
        for (size_t j = 1; j < ny - 2; ++j)
        {
            node0 = idx (i    , j    , ny);
            node1 = idx (i + 1, j    , ny);
            node2 = idx (i + 1, j + 1, ny);
            node3 = idx (i    , j + 1, ny);

            dx_dxi[node0]  = dx_dxi[node1];
            dx_dxi[node3]  = dx_dxi[node2];
            dy_dxi[node0]  = dy_dxi[node1];
            dy_dxi[node3]  = dy_dxi[node2];

            dx_deta[node0]  = dx_deta[node1];
            dx_deta[node3]  = dx_deta[node2];
            dy_deta[node0]  = dy_deta[node1];
            dy_deta[node3]  = dy_deta[node2];

            ddx_dxi2[node0] = ddx_dxi2[node1];
            ddx_dxi2[node3] = ddx_dxi2[node2];
            ddy_dxi2[node0] = ddy_dxi2[node1];
            ddy_dxi2[node3] = ddy_dxi2[node2];

            ddx_deta2[node0] = ddx_deta2[node1];
            ddx_deta2[node3] = ddx_deta2[node2];
            ddy_deta2[node0] = ddy_deta2[node1];
            ddy_deta2[node3] = ddy_deta2[node2];
        }
    }
    // south boundary
    for (size_t i = 1; i < nx - 2; ++i)
    {
        j = 0;  // for (size_t j = 0; j < 1; ++j)
        {
            node0 = idx (i    , j    , ny);
            node1 = idx (i + 1, j    , ny);
            node2 = idx (i + 1, j + 1, ny);
            node3 = idx (i    , j + 1, ny);

            dx_dxi[node0]  = dx_dxi[node3];
            dx_dxi[node1]  = dx_dxi[node2];
            dy_dxi[node0]  = dy_dxi[node3];
            dy_dxi[node1]  = dy_dxi[node2];

            dx_deta[node0]  = dx_deta[node3];
            dx_deta[node1]  = dx_deta[node2];
            dy_deta[node0]  = dy_deta[node3];
            dy_deta[node1]  = dy_deta[node2];

            ddx_dxi2[node0] = ddx_dxi2[node3];
            ddx_dxi2[node1] = ddx_dxi2[node2];
            ddy_dxi2[node0] = ddy_dxi2[node3];
            ddy_dxi2[node1] = ddy_dxi2[node2];

            ddx_deta2[node0] = ddx_deta2[node3];
            ddx_deta2[node1] = ddx_deta2[node2];
            ddy_deta2[node0] = ddy_deta2[node3];
            ddy_deta2[node1] = ddy_deta2[node2];
        }
    }
    // north boundary
    for (size_t i = 1; i < nx - 2; ++i)
    {
        j = ny - 2;
        {
            node0 = idx (i    , j    , ny);
            node1 = idx (i + 1, j    , ny);
            node2 = idx (i + 1, j + 1, ny);
            node3 = idx (i    , j + 1, ny);

            dx_dxi[node2]  = dx_dxi[node1];
            dx_dxi[node3]  = dx_dxi[node0];
            dy_dxi[node2]  = dy_dxi[node1];
            dy_dxi[node3]  = dy_dxi[node0];

            dx_deta[node2]  = dx_deta[node1];
            dx_deta[node3]  = dx_deta[node0];
            dy_deta[node2]  = dy_deta[node1];
            dy_deta[node3]  = dy_deta[node0];

            ddx_dxi2[node3] = ddx_dxi2[node0];
            ddx_dxi2[node2] = ddx_dxi2[node1];
            ddy_dxi2[node3] = ddy_dxi2[node0];
            ddy_dxi2[node2] = ddy_dxi2[node1];

            ddx_deta2[node3] = ddx_deta2[node0];
            ddx_deta2[node2] = ddx_deta2[node1];
            ddy_deta2[node3] = ddy_deta2[node0];
            ddy_deta2[node2] = ddy_deta2[node1];
        }
    }
    // east boundary
    i = nx - 2;  // for (size_t i = nx - 2; i < nx - 1; ++i)
    {
        for (size_t j = 1; j < ny - 2; ++j)
        {
            node0 = idx (i    , j    , ny);
            node1 = idx (i + 1, j    , ny);
            node2 = idx (i + 1, j + 1, ny);
            node3 = idx (i    , j + 1, ny);

            dx_dxi[node1]  = dx_dxi[node0];
            dx_dxi[node2]  = dx_dxi[node3];
            dy_dxi[node1]  = dy_dxi[node0];
            dy_dxi[node2]  = dy_dxi[node3];

            dx_deta[node1]  = dx_deta[node0];
            dx_deta[node2]  = dx_deta[node3];
            dy_deta[node1]  = dy_deta[node0];
            dy_deta[node2]  = dy_deta[node3];

            ddx_dxi2[node1] = ddx_dxi2[node0];
            ddx_dxi2[node2] = ddx_dxi2[node3];
            ddy_dxi2[node1] = ddy_dxi2[node0];
            ddy_dxi2[node2] = ddy_dxi2[node3];

            ddx_deta2[node1] = ddx_deta2[node0];
            ddx_deta2[node2] = ddx_deta2[node3];
            ddy_deta2[node1] = ddy_deta2[node0];
            ddy_deta2[node2] = ddy_deta2[node3];
        }
    }
    // south-west corner
    i = 0;
    j = 0;
    node0 = idx(i, j, ny);
    dx_dxi[node0] = 4. * dx_dxi[node0];
    dy_dxi[node0] = 4. * dy_dxi[node0];
    dx_deta[node0] = 4. * dx_deta[node0];
    dy_deta[node0] = 4. * dy_deta[node0];

    ddx_dxi2 [node0] = 0.0;
    ddy_dxi2 [node0] = 0.0;
    ddx_deta2[node0] = 0.0;
    ddy_deta2[node0] = 0.0;

    // north-west corner
    i = 0;
    j = ny - 2;
    node3 = idx(i, j + 1, ny);
    dx_dxi [node3] = 4. * dx_dxi[node3];
    dy_dxi [node3] = 4. * dy_dxi[node3];
    dx_deta[node3] = 4. * dx_deta[node3];
    dy_deta[node3] = 4. * dy_deta[node3];

    ddx_dxi2 [node3] = 0.0;
    ddy_dxi2 [node3] = 0.0;
    ddx_deta2[node3] = 0.0;
    ddy_deta2[node3] = 0.0;

    // north-east corner
    i = nx - 2;
    j = ny - 2;
    node2 = idx(i + 1, j + 1, ny);
    dx_dxi [node2] = 4. * dx_dxi[node2];
    dy_dxi [node2] = 4. * dy_dxi[node2];
    dx_deta[node2] = 4. * dx_deta[node2];
    dy_deta[node2] = 4. * dy_deta[node2];

    ddx_dxi2 [node2] = 0.0;
    ddy_dxi2 [node2] = 0.0;
    ddx_deta2[node2] = 0.0;
    ddy_deta2[node2] = 0.0;

    // south-east corner
    i = nx - 2;
    j = 0;
    node1 = idx(i + 1, j, ny);
    dx_dxi [node1] = 4. * dx_dxi[node1];
    dy_dxi [node1] = 4. * dy_dxi[node1];
    dx_deta[node1] = 4. * dx_deta[node1];
    dy_deta[node1] = 4. * dy_deta[node1];

    ddx_dxi2 [node1] = 0.0;
    ddy_dxi2 [node1] = 0.0;
    ddx_deta2[node1] = 0.0;
    ddy_deta2[node1] = 0.0;

    return 0;
}
inline size_t idx(size_t i, size_t j, size_t ny_in)
{
    return i * ny_in + j;
}
