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

double scvf_xi(double c0, double c1, double c2, double c3);
double scvf_eta(double c0, double c1, double c2, double c3);
double c_scv(double c0, double c1, double c2, double c3);
double dcdx_scv(double c0, double c1, double c2, double c3);
double dcdy_scv(double c0, double c1, double c2, double c3);
double dcdx_scvf_n(double c0, double c1, double c2, double c3);
double dcdx_scvf_t(double c0, double c1, double c2, double c3);
double dcdy_scvf_n(double c0, double c1, double c2, double c3);
double dcdy_scvf_t(double c0, double c1, double c2, double c3);
