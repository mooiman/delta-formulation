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

#ifndef __BED_LEVEL_H__
#define __BED_LEVEL_H__

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

enum class BED_LEVEL_ENUM
{
    NONE = 0,
    FLAT,
    SHOAL,
    SLOPED,
    WAVY,
    WAVY_SLOPED,
    WEIR,
    NR_BED_LEVELS
};

class BED_LEVEL
{
public:
    BED_LEVEL();
    ~BED_LEVEL();
    int set_bed_level_type(std::string geometry_type, BED_LEVEL_ENUM& bed_level_type);
    int initialize_bed_level(BED_LEVEL_ENUM& bed_type, std::vector<double>& x, std::vector<double>& zb, std::string & model_title, double depth);
//    long open(std::string filename);
//    long read(size_t nx, size_t ny);
//    std::vector<double> get_bed_level();

//private:
//    void transpose(std::vector<double>& x, size_t nx, size_t ny);
//
//    std::ifstream m_fname;
//    std::vector<double> m_bed_given;
};

#endif  // __BED_LEVEL_H__
