//
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//

#include "ugrid1d.h"
#include "include/netcdf.h"

UGRID1D::UGRID1D()
{
    m_ncid = -1;
    m_times = -1;
    m_time_units = "seconds since 2023-01-01 00:00:00";
}
UGRID1D::~UGRID1D()
{
}
int UGRID1D::open(std::string ncfile, std::string model_title)
{
    int status = nc_create(ncfile.c_str(), NC_NETCDF4, &m_ncid);

    const std::chrono::zoned_time now{ std::chrono::current_zone(), std::chrono::system_clock::now() };
    //auto date_time = std::format("{:%FT%H:%M:%OSZ%Oz (%Z)}", now);
    auto date_time = std::format("{:%F %H:%M:%OS %Oz}", now);

    // Define global attributes
    status = set_global_attribute("Title", model_title);
    status = set_global_attribute("Model", "Delta-Formulation 1D, C++");
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
int UGRID1D::close()
{
    int status = nc_close(m_ncid);
    return status;
}

int UGRID1D::network()
{
    int status = -1;
    int i_var;

    std::string var_name("network1D");
    status = nc_def_var(m_ncid, var_name.data(), NC_INT, 0, nullptr, &i_var);
    status = set_attribute(var_name, "cf_role", "mesh_topology");
    status = set_attribute(var_name, "edge_dimension", "nNetworkEdges");
    status = set_attribute(var_name, "edge_geometry", "network1D_geometry");
    status = set_attribute(var_name, "edge_length", "network1D_edge_length");
    status = set_attribute(var_name, "edge_node_connectivity", "network1D_edge_nodes");
    status = set_attribute(var_name, "long_name", "Network topology");
    status = set_attribute(var_name, "node_coordinates", "network1D_node_x network1D_node_y");
    status = set_attribute(var_name, "topology_dimension", 1);
    return status;
}

int UGRID1D::network_geometry()
{
    int status = -1;
    int i_var;

    std::string var_name("network1D_geometry");
    status = nc_def_var(m_ncid, var_name.data(), NC_INT, 0, nullptr, &i_var);
    status = set_attribute(var_name, std::string("geometry_type"), std::string("line"));
    status = set_attribute(var_name, std::string("long_name"), std::string("1D Geometry"));
    status = set_attribute(var_name, std::string("node_coordinates"), std::string("network1D_geom_x network1D_geom_y"));
    status = set_attribute(var_name, std::string("node_dimension"), std::string("nGeometryNodes"));
    status = set_attribute(var_name, std::string("node_count"), std::string("network1D_node_count"));
    return status;
}

int UGRID1D::mesh1d()
{
    int status = -1;
    int i_var;

    std::string var_name("mesh1D");
    status = nc_def_var(m_ncid, var_name.data(), NC_INT, 0, nullptr, &i_var);
    status = set_attribute(var_name, "cf_role", "mesh_topology");
    status = set_attribute(var_name, "coordinate_space", "network1D");
    status = set_attribute(var_name, "edge_coordinates", "mesh1D_edge_branch");
    status = set_attribute(var_name, "edge_dimension", "nMesh1DEdges");
    status = set_attribute(var_name, "edge_node_connectivity", "mesh1D_edge_nodes");
    status = set_attribute(var_name, "long_name", "Mesh 1D");
    status = set_attribute(var_name, "node_coordinates", "mesh1D_node_branch mesh1D_node_offset");
    status = set_attribute(var_name, "node_dimension", "nMesh1DNodes");
    status = set_attribute(var_name, "node_edge_exchange", "mesh1D_connect_node_edges");
    status = set_attribute(var_name, "topology_dimension", 1);
    return status;
}

int UGRID1D::def_dimensions(int nr_ntw_nodes, int nr_ntw_edges, int nr_geom_nodes, int nr_mesh_nodes, int nr_mesh_edges)
{
    int dim_id;
    int status = -1;
    status = nc_def_dim(m_ncid, "nNetworkNodes", nr_ntw_nodes, &dim_id);
    status = nc_def_dim(m_ncid, "nNetworkEdges", nr_ntw_edges, &dim_id);
    status = nc_def_dim(m_ncid, "nGeometryNodes", nr_geom_nodes, &dim_id);
    status = nc_def_dim(m_ncid, "nMesh1DNodes", nr_mesh_nodes, &dim_id);
    status = nc_def_dim(m_ncid, "nMesh1DEdges", nr_mesh_edges, &dim_id);
    status = nc_def_dim(m_ncid, "Two", 2, &dim_id);
    return status;
}
int UGRID1D::add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;

}
int UGRID1D::add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit, std::string mesh, std::string location)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("mesh"), mesh);
    status = set_attribute(var_name, std::string("location"), location);
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}
int UGRID1D::add_ntw_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string long_name)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("cf_role"), cf_role);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}
int UGRID1D::add_geom_coordinate(std::string var_name, std::vector<std::string> dim_names, std::string cf_role, std::string std_name, std::string long_name, std::string unit)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("cf_role"), cf_role);
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}

int UGRID1D::add_node_count(std::string var_name, std::vector<std::string> dim_names, std::string long_name)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}

int UGRID1D::add_mesh1d_edge_nodes(std::string var_name, std::vector<std::string> dim_names, std::string long_name)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}

int UGRID1D::add_mesh1d_point_on_branch(std::string var_name, std::vector<std::string> dim_names, std::string long_name)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_UINT, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    return status;
}
int UGRID1D::add_mesh1d_offset_on_branch(std::string var_name, std::vector<std::string> dim_names, std::string long_name, std::string unit)
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

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}
int UGRID1D::add_time_series(void)
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
int UGRID1D::put_time(const int nst, double time)
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
int UGRID1D::put_time_variable(std::string var_name, int i_time, std::vector<double> values)
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
int UGRID1D::put_variable(std::string var_name, std::vector<double> values)
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
int UGRID1D::put_variable(std::string var_name, std::vector<int> values)
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
int UGRID1D::put_variable_2(std::string var_name, std::vector<int> values)
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

////////////////////////////////////////////////////////////////////////////////
int UGRID1D::set_global_attribute(std::string att_name, std::string att_value)
{
    int status = nc_put_att_text(m_ncid, NC_GLOBAL, att_name.data(), att_value.size(), att_value.data());
    return status;
}
int UGRID1D::set_attribute(std::string var_name, std::string att_name, double att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_double(m_ncid, i_var, att_name.data(), NC_DOUBLE, 1, &att_value);
    return status;
}
int UGRID1D::set_attribute(std::string var_name, std::string att_name, int att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_int(m_ncid, i_var, att_name.data(), NC_INT, 1, &att_value);
    return status;
}
int UGRID1D::set_attribute(std::string var_name, std::string att_name, std::string att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_text(m_ncid, i_var, att_name.data(), att_value.size(), att_value.data());
    return status;
}