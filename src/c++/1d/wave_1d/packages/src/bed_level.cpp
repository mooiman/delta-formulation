//
// Programmer: Jan Mooiman
// Email     : jan.mooiman@outlook.com
//
//    Solving the 1D shallow water equations, fully implicit with delta-formulation and Modified Newton iteration 
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
// set bed level type
size_t BED_LEVEL::set_bed_level_type(std::string geometry_type, BED_LEVEL_ENUM& bed_level_type)
{
    if (geometry_type == "WavyBedGeometry") { bed_level_type = BED_LEVEL_ENUM::WAVY; }
    else if (geometry_type == "SlopedBedGeometry") { bed_level_type = BED_LEVEL_ENUM::SLOPED; }
    else if (geometry_type == "SlopedWavyBedGeometry") { bed_level_type = BED_LEVEL_ENUM::WAVY_SLOPED; }
    else if (geometry_type == "UniformGeometry") {bed_level_type = BED_LEVEL_ENUM::FLAT; }
    else if (geometry_type == "WeirBorsboom2000Geometry") { bed_level_type = BED_LEVEL_ENUM::WEIR; }
    else if (geometry_type == "Shoal") { bed_level_type = BED_LEVEL_ENUM::SHOAL; }
    else { bed_level_type = BED_LEVEL_ENUM::NONE; }

    return 0;
}

int BED_LEVEL::initialize_bed_level(BED_LEVEL_ENUM& bed_type, std::vector<double>& x, std::vector<double>& zb, std::string & model_title, double depth)
{
    int status = 0;
    if (bed_type == BED_LEVEL_ENUM::NONE)
    {
        std::cout << "Bed level type BED_LEVEL::NONE not supported" << std::endl;
        status = 1;
    }
    else if (bed_type == BED_LEVEL_ENUM::FLAT)
    {
        std::stringstream depth_strm;
        depth_strm << std::fixed << std::setprecision(2) << depth;
        model_title += "\nFlat bed level at " + depth_strm.str() + " [m]";
        for (int i = 0; i < x.size(); ++i)
        {
            zb[i] = -depth;
        }
    }
    else if (bed_type == BED_LEVEL_ENUM::SHOAL)
    {
        model_title = "Shoaling bedlevel from -10 [m] to -2.5 [m]";
        double slope_begin = 0.0;
        double slope_end = 1000.;
        for (int i = 0; i < x.size(); ++i)
        {
            if (x[i] < slope_begin)
                zb[i] = -10.0;
            else if(x[i] < slope_end)
                zb[i] = -10. + (x[i] - slope_begin) / (slope_end - slope_begin) * 7.5;
            else
                zb[i] = -2.5;
        }
    }
    else if (bed_type == BED_LEVEL_ENUM::SLOPED)
    {
        model_title = "Sloped bed level 10 [cm] per [km] (par 5.1, Platzek2019)";
        // From -3.9 [m] to -4.0 [m]
        double ib = -1e-04;  // 0.1 mm per meter
        double slope_begin = x[1];
        double x0 = x[1];
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] - x0;
            zb[i] = -3.9 + (x[i] - slope_begin) * ib;
        }
    }
    else if (bed_type == BED_LEVEL_ENUM::WEIR)
    {
        // Borsboom_development1Derrorminmovingadaptgridmethod_AdaptMethodLinesCRC2001.pdf
        model_title += "\nWeir: from -12 [m] to -5 [m] and from -5 [m] to -10 [m] (Borsboom_presCASATUE2023)";
        double x0 = x[1];  // x[1] = -250. 
        double zb_def = -12.0;

        double slope_up_begin = 200.;
        double slope_up_end = 250.;
        double slope_down_begin = 350.;
        double slope_down_end = 450.;
        double zb_begin = zb_def;
        double zb_weir = zb_def + 7.0;
        double zb_end = zb_def + 2.0;
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] - x0;
            if (x[i] < slope_up_begin)
                zb[i] = zb_begin;
            else if (x[i] < slope_up_end)
                zb[i] = zb_begin + (x[i] - slope_up_begin) / (slope_up_end - slope_up_begin) * (zb_weir - zb_begin);
            else if (x[i] < slope_down_begin)
                zb[i] = zb_weir;
            else if (x[i] < slope_down_end)
                zb[i] = zb_weir + (x[i] - slope_down_begin) / (slope_down_end - slope_down_begin) * (zb_end - zb_weir);
            else if (x[i] >= slope_down_end)
                zb[i] = zb_end;
        }
    }
    else if (bed_type == BED_LEVEL_ENUM::WAVY)
    {
        model_title = "Wavy bed level (par 5.2, Platzek2019)";
        int nx = (int) x.size();
        double Lx = x[nx - 2] - x[1];  // 2 virtual nodes
        double dx = Lx / (nx - 3);
        double zb_ref = -4.0;  // Mean depth for the wave bed level
        double Ab = 0.3;    // Amplitude of the bed forms
        double Nb = 25.0;  // Number of waves in bed level in the domain
        double fb = 2.0 * M_PI * Nb / ((nx-3) * dx);
        double x0 = x[1];
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] - x0;
            zb[i] = zb_ref - Ab * cos(fb * (x[i]- x[1]));
        }
        zb[0] = zb[1];
        zb[nx - 1] = zb[nx - 2];
    }
    else if (bed_type == BED_LEVEL_ENUM::WAVY_SLOPED)
    {
        model_title = "Sloped wavy bed level (par 5.3, Platzek2019)";
        int nx = (int)x.size();
        double Lx = x[nx - 2] - x[1];
        double dx = Lx / (nx - 3);
        double zb_ref = -3.9;  // Mean depth for the wave bed level
        double Ab = 0.3;    // Amplitude of the bed forms
        double Nb = 25.0;  // Number of waves in bed level in the domain
        double fb = 2.0 * M_PI * Nb / ((nx - 3) * dx);
        double ib = -1e-04;  // 0.1 mm per meter
        double x0 = x[1];
        for (int i = 0; i < x.size(); ++i)
        {
            x[i] = x[i] - x0;
            zb[i] = zb_ref - Ab * cos(fb * (x[i] - x[1])) + ib * (x[i] - x[1]);
        }
        zb[0] = zb[1];
        zb[nx - 1] = zb[nx - 2];
    }
    else
    {
        std::cout << "----------------------------" << std::endl;
        std::cout << "Bed level type not supported" << std::endl;
        std::cout << "Press Enter to finish";
        std::cin.ignore();
        status = 1; 
    }
    return status;
}
