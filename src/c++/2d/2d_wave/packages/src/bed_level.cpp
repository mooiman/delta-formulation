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

#include "bed_level.h"

//------------------------------------------------------------------------------
BED_LEVEL::BED_LEVEL()
{
}
//------------------------------------------------------------------------------
BED_LEVEL::~BED_LEVEL()
{
}
//------------------------------------------------------------------------------
long BED_LEVEL::open(std::string filename)
{
    long status = 1;
    m_fname.open(filename);
    if (m_fname.is_open()) {
        status = 0;
    }
    return status;
}
//------------------------------------------------------------------------------
long BED_LEVEL::read(int nx, int ny)
{
    long status = -1;
    std::istringstream iss;
    double missing_value = -999.0;  // used by quickin

    m_bed_given.resize(nx*ny);

    int k = -1;
    std::string token;
    std::string line;
    std::getline(m_fname, line);
    iss.clear();
    iss.str(line);
    for (;;)
    {
        while (iss >> token)
        {
            if (missing_value != stod(token))
            {
                ++k;
                m_bed_given[k] = stod(token);
            }
        }
        std::getline(m_fname, line);
        iss.clear();
        iss.str(line);
        if (std::fmod(k + 1, nx*ny) == 0) { break; }
    }
    if (k + 1 == nx * ny)
    {
        status = 0;
    }
    transpose(m_bed_given, nx, ny);
    return status;

}
//------------------------------------------------------------------------------
std::vector<double> BED_LEVEL::get_bed_level()
{
     return m_bed_given;
}

//
////  PRIVATE  ////////////////////////////////////////////////////////
//
void BED_LEVEL::transpose(std::vector<double>& x, int nx, int ny)
{
    std::vector<double> xt;
    xt.resize(nx * ny);
    int k1;
    int k2;
    double tmp;
    for (int i = 0; i < nx; ++i) 
    {
        for (int j = 0; j < ny; ++j) 
        {
            k1 = i * ny + j;
            k2 = j * nx + i;
            xt[k1] = x[k2];
        }
    }
    for (int i = 0; i < nx * ny; ++i)
    {
        x[i] = xt[i];
    }
}

