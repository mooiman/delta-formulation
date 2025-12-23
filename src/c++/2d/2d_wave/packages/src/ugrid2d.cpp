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

#include "ugrid2d.h"
#include "compile_date_and_time.h"
#include "interpolations.h"
#include "main_version.h"
#include "include/netcdf.h"

UGRID2D::UGRID2D()
{
    m_ncid = -1;
    m_times = -1;
    m_time_units = "seconds since 2007-01-01 00:00:00";
}
UGRID2D::~UGRID2D()
{
}
std::string UGRID2D::compileDateTime()
{
    std::string str1(compileYear());
    std::string str2(compileMonth());
    std::string str3(compileDay());
    if (str3.size() == 1)
    {
        str3.resize(2);
        str3[1] = str3[0];
        str3[0] = '0';
    }
    return str1 + "-" + str2 + "-" + str3 + " " + __TIME__;
}
int UGRID2D::open(std::string ncfile, std::string model_title)
{
    int status = nc_create(ncfile.c_str(), NC_NETCDF4, &m_ncid);

    const std::chrono::zoned_time now{ std::chrono::current_zone(), std::chrono::system_clock::now() };
    //auto date_time = std::format("{:%FT%H:%M:%OSZ%Oz (%Z)}", now);
    auto date_time = std::format("{:%F %H:%M:%OS %Oz}", now);

    // Define global attributes
    status = set_global_attribute("Title", model_title);
    status = set_global_attribute("Model", "Delta-formulation 2D, C++");
    status = set_global_attribute("Program created", compileDateTime() );
    status = set_global_attribute("Program version", getversionstring_main() );
    status = set_global_attribute("Program build", getbuildstring_main() );
    status = set_global_attribute("Conventions", "CF-1.8 UGRID-1.0");
    status = set_global_attribute("file_created", date_time);
    status = set_global_attribute("reference", "https://www.github.com/mooiman");

    int var_id;
    status = nc_def_var(m_ncid, "projected_coordinate_system", NC_INT, 0, nullptr, &var_id);
    int epsg = 28992;
    std::string epsg_code = "EPSG:28992";
    status = set_attribute("projected_coordinate_system", "name", "Unknown projected");
    status = set_attribute("projected_coordinate_system", "epsg", epsg);
    status = set_attribute("projected_coordinate_system", "grid_mapping_name", "Unknown projected");
    status = set_attribute("projected_coordinate_system", "longitude_of_prime_meridian", 0.);
    status = set_attribute("projected_coordinate_system", "semi_major_axis", 6378137.);
    status = set_attribute("projected_coordinate_system", "semi_minor_axis", 6356752.314245);
    status = set_attribute("projected_coordinate_system", "inverse_flattening", 298.257223563);
    status = set_attribute("projected_coordinate_system", "EPSG_CODE", epsg_code);
    status = set_attribute("projected_coordinate_system", "value", "value is equal to EPSG code");

    return status;
}
int UGRID2D::close()
{
    int status = nc_close(m_ncid);
    return status;
}

int UGRID2D::mesh2d()
{
    int status = -1;
    int i_var;

    std::string var_name("mesh2D");
    status = nc_def_var(m_ncid, var_name.data(), NC_INT, 0, nullptr, &i_var);
    status = set_attribute(var_name, "cf_role", "mesh_topology");
    status = set_attribute(var_name, "edge_coordinates", "mesh2d_edge_x mesh2d_edge_y");
    status = set_attribute(var_name, "edge_dimension", "mesh2d_nEdges");
    status = set_attribute(var_name, "edge_face_connectivity", "mesh2d_edge_faces");
    status = set_attribute(var_name, "edge_node_connectivity", "mesh2d_edge_nodes");
    status = set_attribute(var_name, "face_coordinates", "mesh2d_face_x mesh2d_face_y");
    status = set_attribute(var_name, "long_name", "Mesh 2D");
    status = set_attribute(var_name, "max_face_nodes_dimension", "mesh2d_nMax_face_nodes");
    status = set_attribute(var_name, "face_dimension", "mesh2d_nFaces");
    status = set_attribute(var_name, "face_node_connectivity", "mesh2d_face_nodes");
    status = set_attribute(var_name, "node_coordinates", "mesh2d_node_x mesh2d_node_y");
    status = set_attribute(var_name, "node_dimension", "mesh2d_nNodes");
    status = set_attribute(var_name, "topology_dimension", 2);
    return status;
}

int UGRID2D::def_dimensions(int nr_mesh_nodes, int nr_mesh_edges, int nr_mesh_faces, int mesh2d_nmax_face_nodes)
{
    int dim_id;
    int status = -1;
    status = nc_def_dim(m_ncid, "mesh2d_nEdges", nr_mesh_edges, &dim_id);
    status = nc_def_dim(m_ncid, "mesh2d_nFaces", nr_mesh_faces, &dim_id);
    status = nc_def_dim(m_ncid, "mesh2d_nNodes", nr_mesh_nodes, &dim_id);
    status = nc_def_dim(m_ncid, "mesh2d_nMax_face_nodes", mesh2d_nmax_face_nodes, &dim_id);
    
    status = nc_def_dim(m_ncid, "Two", 2, &dim_id);
    return status;
}
int UGRID2D::add_edge_nodes(size_t nx, size_t ny)
{
    int status;
    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nEdges");
    dim_names.push_back("Two");
    status = this->add_mesh2d_edge_nodes("mesh2d_edge_nodes", dim_names, "Each edge connects two nodes");

    size_t nr_nodes = nx * ny;
    //int nr_edges = (nx - 1) * ny + nx * (ny - 1);
    std::vector<int> node_mask(nr_nodes, 0);

    size_t p0 = 0;
    size_t p1 = 0;
    for (size_t i = 0; i < nx - 1; ++i)
    {
        for (size_t j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j, ny);
            p1 = idx(i + 1, j, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                m_mesh2d_edge_nodes.push_back(p0);
                m_mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i + 1, j    , ny);
            p1 = idx(i + 1, j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                m_mesh2d_edge_nodes.push_back(p0);
                m_mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i + 1, j + 1, ny);
            p1 = idx(i    , j + 1, ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                m_mesh2d_edge_nodes.push_back(p0);
                m_mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                node_mask[p1] += 1;
            }

            p0 = idx(i    , j + 1, ny);
            p1 = idx(i    , j    , ny);
            if (node_mask[p0] <= 1 || node_mask[p1] <= 1)
            {
                m_mesh2d_edge_nodes.push_back(p0);
                m_mesh2d_edge_nodes.push_back(p1);
                node_mask[p0] += 1;
                //node_mask[p1] += 1;
            }
        }
    }
    status = this->put_variable_2("mesh2d_edge_nodes", m_mesh2d_edge_nodes);

    return status;
}
int UGRID2D::add_face_nodes(std::vector<double> & x, std::vector<double> & y, double fill_value, size_t nx, size_t ny)
{
    int status;

    size_t p0;
    size_t p1;
    size_t p2;
    size_t p3;
    std::vector<int> mesh2d_face_nodes;

    for (int i = 0; i < nx - 1; ++i)
    {
        for (int j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            if ((x[p0] != fill_value || y[p0] != fill_value) &&
                (x[p1] != fill_value || y[p1] != fill_value) &&
                (x[p2] != fill_value || y[p2] != fill_value) &&
                (x[p3] != fill_value || y[p3] != fill_value) )
            {
                mesh2d_face_nodes.push_back(p0);
                mesh2d_face_nodes.push_back(p1);
                mesh2d_face_nodes.push_back(p2);
                mesh2d_face_nodes.push_back(p3);
            }
        }
    }
    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nFaces");
    dim_names.push_back("mesh2d_nMax_face_nodes");

    status = this->add_mesh2d_edge_nodes("mesh2d_face_nodes", dim_names, "Each face contains four nodes");
    status = this->put_variable_4("mesh2d_face_nodes", mesh2d_face_nodes);

    return status;
}
int UGRID2D::add_node_coords(std::vector<double> & x, std::vector<double> & y, double fill_value)
{
    int status = 1;
    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nNodes");
    status = this->add_variable("mesh2d_node_x", dim_names, "projection_x_coordinate", "x", "m");
    status = this->add_attribute("mesh2d_node_x", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_node_x", x);

    status = this->add_variable("mesh2d_node_y", dim_names, "projection_y_coordinate", "y", "m");
    status = this->add_attribute("mesh2d_node_y", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_node_y", y);

    return status;
}
int UGRID2D::add_edge_coords(std::vector<double> & x, std::vector<double> & y, double fill_value)
{
    int status = 1;

    // Compute edges centres
    int p0;
    int p1;
    std::vector<double> edge_x;
    std::vector<double> edge_y;
    for (int i = 0; i < m_mesh2d_edge_nodes.size(); i += 2)
    {
        p0 = m_mesh2d_edge_nodes[i];
        p1 = m_mesh2d_edge_nodes[i+1];
        edge_x.push_back(0.5 * (x[p0] + x[p1]));
        edge_y.push_back(0.5 * (y[p0] + y[p1]));
    }

    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nEdges");
    status = this->add_variable("mesh2d_edge_x", dim_names, "projection_x_coordinate", "x", "m");
    status = this->add_attribute("mesh2d_edge_x", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_edge_x", edge_x);

    status = this->add_variable("mesh2d_edge_y", dim_names, "projection_y_coordinate", "y", "m");
    status = this->add_attribute("mesh2d_edge_y", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_edge_y", edge_y);

    return status;
}
int UGRID2D::add_face_mass_centres(std::vector<double> & x, std::vector<double> & y, double fill_value, size_t nx, size_t ny)
{
    int status = 1;

    // Compute mass centres of faces
    size_t p0;
    size_t p1;
    size_t p2;
    size_t p3;
    std::vector<double> xmc;
    std::vector<double> ymc;
    for (size_t i = 0; i < nx - 1; ++i)
    {
        for (size_t j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            xmc.push_back(0.25 * (x[p0] + x[p1] + x[p2] + x[p3]));
            ymc.push_back(0.25 * (y[p0] + y[p1] + y[p2] + y[p3]));
        }
    }
    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nFaces");
    status = this->add_variable("mesh2d_face_x", dim_names, "projection_x_coordinate", "x", "m");
    status = this->add_attribute("mesh2d_face_x", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_face_x", xmc);

    status = this->add_variable("mesh2d_face_y", dim_names, "projection_y_coordinate", "y", "m");
    status = this->add_attribute("mesh2d_face_y", "_FillValue", fill_value);
    status = this->put_variable("mesh2d_face_y", ymc);

    return status;
}
int UGRID2D::add_face_area(std::vector<double> & x, std::vector<double> & y, double fill_value, size_t nx, size_t ny)
{
    int status = 1;

    // Compute area of faces
    size_t p0;
    size_t p1;
    size_t p2;
    size_t p3;
    std::vector<double> cell_area;
    std::vector<double> x_pol(4);
    std::vector<double> y_pol(4);

    for (size_t i = 0; i < nx - 1; ++i)
    {
        for (size_t j = 0; j < ny - 1; ++j)
        {
            p0 = idx(i    , j    , ny);
            p1 = idx(i + 1, j    , ny);
            p2 = idx(i + 1, j + 1, ny);
            p3 = idx(i    , j + 1, ny);
            //
            // Area: $J = x_\xi y_\eta - y_\xi x_\eta$
            //
            x_pol[0] = x[p0];
            x_pol[1] = x[p1];
            x_pol[2] = x[p2];
            x_pol[3] = x[p3];
            y_pol[0] = y[p0];
            y_pol[1] = y[p1];
            y_pol[2] = y[p2];
            y_pol[3] = y[p3];
            double area = polygon_area(x_pol, y_pol);
            cell_area.push_back(area);
        }
    }

    std::vector<std::string> dim_names;
    dim_names.push_back("mesh2d_nFaces");
    status = this->add_variable("cell_area", dim_names, "cell_area", "-", "m2", "mesh2D", "face", "Area of an element bounded by the vertices");
    status = this->add_attribute("cell_area", "coordinates", "mesh2d_face_x, mesh2d_face_y");
    status = this->put_variable("cell_area", cell_area);

    return status;
}

int UGRID2D::add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}
int UGRID2D::add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, 
    std::string unit, std::string mesh, std::string location, std::string comment)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("mesh"), mesh);
    status = set_attribute(var_name, std::string("location"), location);
    if (location == "node")
    {
        std::string xy("mesh2d_node_x mesh2d_node_y");
        status = set_attribute(var_name, std::string("coordinates"), xy);
    }
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    if (comment.size() != 0) { status = set_attribute(var_name, std::string("comment"), comment); }

    return status;
}
int UGRID2D::add_attribute(std::string var_name, std::string att_name, std::string att_value)
{
    int status = -1;
    status = set_attribute(var_name, att_name, att_value);
    return status;
}
int UGRID2D::add_attribute(std::string var_name, std::string att_name, double att_value)
{
    int status = -1;
    status = set_attribute(var_name, att_name, att_value);
    return status;
}

int UGRID2D::add_ntw_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string long_name)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("cf_role"), cf_role);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}
int UGRID2D::add_geom_coordinate(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string std_name, std::string long_name, std::string unit)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("cf_role"), cf_role);
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}

int UGRID2D::add_node_count(std::string var_name, std::vector<std::string> dim_names, std::string long_name)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}

int UGRID2D::add_mesh2d_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string long_name)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    for (int i = 0; i < dim_names.size(); ++i)
    {
        status = nc_inq_dimid(m_ncid, dim_names[i].data(), &dim_id);
        dimids.push_back(dim_id);
    }

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, (int)dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}

int UGRID2D::add_time_series(void)
{
    int dim_id;
    int i_var;
    int status = nc_def_dim(m_ncid, "time", NC_UNLIMITED, &dim_id);
    std::vector<int> dimids;
    dimids.push_back(dim_id);
    status = nc_def_var(m_ncid, "time", NC_DOUBLE, 1, &dimids[0], &i_var);
    status = set_attribute(std::string("time"), std::string("standard_name"), std::string("time"));
    status = set_attribute(std::string("time"), std::string("long_name"), std::string("Time"));
    status = set_attribute(std::string("time"), std::string("units"), m_time_units);
    return status;
}
int UGRID2D::put_time(const int nst, double time)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if ("time" == tmp_var_name)
        {
            size_t start1[] = { (size_t) nst };
            size_t count1[] = { 1 };

            status = nc_put_vara_double(m_ncid, i_var, start1, count1, (const double*)&time);
            break;
        }
    }
    free(tmp_var_name_c);
    return status;
}
int UGRID2D::put_time_variable(std::string var_name, int i_time, std::vector<double> values)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if (var_name == tmp_var_name)
        {
            size_t start1[] = { (size_t) i_time, 0 };
            size_t count1[] = { 1, values.size() };

            status = nc_put_vara_double(m_ncid, i_var, start1, count1, (const double*)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}
int UGRID2D::put_variable(std::string var_name, std::vector<double> values)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if (var_name == tmp_var_name)
        {
            size_t start1[] = { 0 };
            size_t count1[] = { values.size()};

            status = nc_put_vara_double(m_ncid, i_var, start1, count1, (const double*)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}
int UGRID2D::put_variable(std::string var_name, std::vector<int> values)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if (var_name == tmp_var_name)
        {
            size_t start1[] = { 0 };
            size_t count1[] = { values.size() };

            status = nc_put_vara_int(m_ncid, i_var, start1, count1, (int*)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}
int UGRID2D::put_variable_2(std::string var_name, std::vector<int> values)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if (var_name == tmp_var_name)
        {
            size_t start1[] = { 0, 0 };
            size_t count1[] = { values.size()/2, 2 };

            status = nc_put_vara_int(m_ncid, i_var, start1, count1, (const int*)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}

int UGRID2D::put_variable_4(std::string var_name, std::vector<int> values)
{
    int status = -1;
    int ndims, nvars, natts, nunlimited;

    char* tmp_var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    for (int i_var = 0; i_var < nvars; ++i_var)
    {
        status = nc_inq_varname(this->m_ncid, i_var, tmp_var_name_c);
        std::string tmp_var_name(tmp_var_name_c);

        if (var_name == tmp_var_name)
        {
            size_t start1[] = { 0, 0 };
            size_t count1[] = { values.size() / 4, 4 };

            status = nc_put_vara_int(m_ncid, i_var, start1, count1, (const int*)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}
////////////////////////////////////////////////////////////////////////////////
inline size_t UGRID2D::idx(size_t i, size_t j, size_t ny)
{
    return i * ny + j;
}
int UGRID2D::set_global_attribute(std::string att_name, std::string att_value)
{
    int status = nc_put_att_text(m_ncid, NC_GLOBAL, att_name.data(), att_value.size(), att_value.data());
    return status;
}
int UGRID2D::set_attribute(std::string var_name, std::string att_name, double att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_double(m_ncid, i_var, att_name.data(), NC_DOUBLE, 1, &att_value);
    return status;
}
int UGRID2D::set_attribute(std::string var_name, std::string att_name, int att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_int(m_ncid, i_var, att_name.data(), NC_INT, 1, &att_value);
    return status;
}
int UGRID2D::set_attribute(std::string var_name, std::string att_name, std::string att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_text(m_ncid, i_var, att_name.data(), att_value.size(), att_value.data());
    return status;
}