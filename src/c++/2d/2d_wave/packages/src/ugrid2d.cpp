
#include "ugrid2d.h"
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
int UGRID2D::open(std::string ncfile)
{
    int status = nc_create(ncfile.c_str(), NC_NETCDF4, &m_ncid);

    const auto now = std::chrono::system_clock::now();
    //auto date_time = std::format("{:%FT%H:%M:%OSZ%Oz (%Z)}", now);
    auto date_time = std::format("{:%F %H:%M:%OS %Oz}", now);

    // Define global attributes
    status = set_global_attribute("Conventions", "CF-1.8 UGRID-1.0");
    status = set_global_attribute("institution", "Deltares");
    status = set_global_attribute("file_created", date_time );
    status = set_global_attribute("reference", "https://www.deltares.nl");

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
int UGRID2D::add_variable(std::string var_name, std::vector<std::string> dim_names, std::string std_name, std::string long_name, std::string unit, std::string mesh, std::string location)
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
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
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