//---------------------------------------------------------------
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//
//    Solving the 1D advection/diffusion equation, fully implicit with delta-formuation and Modified Newton iteration 
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

#include "include/netcdf.h"
#include "definition_map_file.h"

UGRID1D * create_map_file(std::string nc_mapfile, std::string model_title, std::vector<double> & x, std::vector<std::string>& map_names)
{

    int status = -1;
    int nr_nodes = (int)x.size();
    int nr_edges = nr_nodes - 1;
    std::vector<double> y(nr_nodes, 0.);                      // y-coordinate

    UGRID1D * map_file = new UGRID1D();
    status = map_file->open(nc_mapfile, model_title);
    status = map_file->network();

    status = map_file->def_dimensions(2, 1, 2, nr_nodes, nr_edges);
    std::vector<std::string> dim_names;

    dim_names.push_back("nNetworkEdges");
    status = map_file->add_variable("network1D_edge_length", dim_names, "-", "Real length of branch geometries", "m");
    std::vector<double> ntw_edge_length{ x[nr_nodes - 1] - x[0] };
    status = map_file->put_variable("network1D_edge_length", ntw_edge_length);

    dim_names.clear();
    dim_names.push_back("nNetworkEdges");
    dim_names.push_back("Two");
    status = map_file->add_ntw_edge_nodes("network1D_edge_nodes", dim_names, "edge_node_connectivity", "Start and end nodes of each branch in the network");
    std::vector<int> ntw_edge_nodes{ 0, 1 };
    status = map_file->put_variable_2("network1D_edge_nodes", ntw_edge_nodes);

    dim_names.clear();
    dim_names.push_back("nNetworkNodes");
    status = map_file->add_variable("network1D_node_x", dim_names, "projection_x_coordinate", "x-coordinates of the network connection nodes", "m");
    status = map_file->add_variable("network1D_node_y", dim_names, "projection_y_coordinate", "y-coordinates of the network connection nodes", "m");
    std::vector<double> x_ntw{ x[0], x[nr_nodes - 1] };
    std::vector<double> y_ntw{ y[0], y[nr_nodes - 1] };
    status = map_file->put_variable("network1D_node_x", x_ntw);
    status = map_file->put_variable("network1D_node_y", y_ntw);


    status = map_file->network_geometry();

    dim_names.clear();
    dim_names.push_back("nGeometryNodes");
    status = map_file->add_geom_coordinate("network1D_geom_x", dim_names, "geometry_x_node", "projection_x_coordinate", "x-coordinates of the branch geometries", "m");
    status = map_file->add_geom_coordinate("network1D_geom_y", dim_names, "geometry_y_node", "projection_y_coordinate", "y-coordinates of the branch geometries", "m");
    std::vector<double> x_geom{ x[0], x[nr_nodes - 1] };
    std::vector<double> y_geom{ y[0], y[nr_nodes - 1] };

    status = map_file->put_variable("network1D_geom_x", x_geom);
    status = map_file->put_variable("network1D_geom_y", y_geom);

    dim_names.clear();
    dim_names.push_back("nNetworkEdges");
    status = map_file->add_node_count("network1D_node_count", dim_names, "Number of geometry nodes per branch (including end points)");
    std::vector<int> ntw_edges{ (int)x_geom.size() };
    status = map_file->put_variable("network1D_node_count", ntw_edges);


    status = map_file->mesh1d();

    dim_names.clear();
    dim_names.push_back("nMesh1DEdges");
    dim_names.push_back("Two");
    status = map_file->add_mesh1d_edge_nodes("mesh1D_edge_nodes", dim_names, "Each edge connects two nodes");
    std::vector<int> mesh1d_edge_nodes;
    for (int i = 0; i < nr_edges; ++i)
    {
        mesh1d_edge_nodes.push_back(i);
        mesh1d_edge_nodes.push_back(i + 1);
    }
    status = map_file->put_variable_2("mesh1D_edge_nodes", mesh1d_edge_nodes);


    dim_names.clear();
    dim_names.push_back("nMesh1DEdges");
    status = map_file->add_mesh1d_point_on_branch("mesh1D_edge_branch", dim_names, "Number of branch on which edge is located");
    std::vector<int> mesh1d_edge_branch;
    for (int i = 0; i < nr_edges; ++i)
    {
        mesh1d_edge_branch.push_back(0);
    }
    status = map_file->put_variable("mesh1D_edge_branch", mesh1d_edge_branch);


    dim_names.clear();
    dim_names.push_back("nMesh1DNodes");
    status = map_file->add_mesh1d_point_on_branch("mesh1D_node_branch", dim_names, "Number of branch on which node is located");
    std::vector<int> mesh1d_node_branch;
    for (int i = 0; i < nr_nodes; ++i)
    {
        mesh1d_node_branch.push_back(0);
    }
    status = map_file->put_variable("mesh1D_node_branch", mesh1d_node_branch);

    dim_names.clear();
    dim_names.push_back("nMesh1DNodes");
    status = map_file->add_mesh1d_offset_on_branch("mesh1D_node_offset", dim_names, "Offset along the branch at which the node is located", "m");
    std::vector<double> mesh1d_offset_on_branch;
    for (int i = 0; i < nr_nodes; ++i)
    {
        mesh1d_offset_on_branch.push_back(x[i] - x[0]);
    }
    status = map_file->put_variable("mesh1D_node_offset", mesh1d_offset_on_branch);

    status = map_file->add_time_series();

    dim_names.clear();
    dim_names.push_back("time");
    dim_names.push_back("nMesh1DNodes");

    status = map_file->add_variable(map_names[0], dim_names, "", "Phi", "-", "mesh1D", "node");
    status = map_file->add_variable(map_names[1], dim_names, "", "Water velocity", "m s-1", "mesh1D", "node");
    status = map_file->add_variable(map_names[2], dim_names, "-", "Viscosity (reg)", "m2 s-1", "mesh1D", "node");
    status = map_file->add_variable(map_names[3], dim_names, "-", "Psi", "m2 s-1", "mesh1D", "node");
    status = map_file->add_variable(map_names[4], dim_names, "-", "Peclet", "-", "mesh1D", "node");
    status = map_file->add_variable(map_names[5], dim_names, "-", "Concentration", "-", "mesh1D", "node");

    return map_file;
}