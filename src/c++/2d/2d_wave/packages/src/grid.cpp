#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>

#include <direct.h>
#include <windows.h>

#define strdup _strdup
#define strcat strcat_s
#define getcwd _getcwd
#define FILE_NAME_MAX _MAX_PATH

#include "grid.h"
#include "include/netcdf.h"
#include "perf_timer.h"

//------------------------------------------------------------------------------
SGRID::SGRID()
{
    _two = 2;
}
//------------------------------------------------------------------------------
SGRID::~SGRID()
{
    int status;
    status = nc_close(m_ncid);
}
//------------------------------------------------------------------------------
long SGRID::open(std::string filename)
{
    long ret_value = 1;
    m_grid_file_name = filename;

    int status = nc_open(m_grid_file_name.c_str(), NC_NOWRITE, &this->m_ncid);
    if (status != NC_NOERR)
    {
        fprintf(stderr, "SGRID::read()\n    Failed to open file: %s\n", m_grid_file_name.c_str());
        return ret_value;
    }
#ifdef NATIVE_C
    fprintf(stderr, "SGRID::read()\n    Opened: %s\n", m_grid_file_name.c_str());
#endif    

    status = this->read_global_attributes();
    for (int i = 0; i < this->global_attributes->count; ++i)
    {
        if (global_attributes->attribute[i]->cvalue.find("UGRID-1") != std::string::npos)
        {
            m_ftype = FILE_TYPE::UGRID;
            ret_value = 1;
        }
        else if (global_attributes->attribute[i]->cvalue.find("SGRID-0.3") != std::string::npos)
        {
            m_ftype = FILE_TYPE::SGRID;
            ret_value = 0;
        }
        else if (global_attributes->attribute[i]->cvalue.find("CF-1.8") != std::string::npos)
        {
            m_ftype = FILE_TYPE::KISS;
            ret_value = 1;
        }
    }
    ret_value = 0;
    return ret_value;
}
//------------------------------------------------------------------------------
long SGRID::read()
{
    int status = -1;
    if (m_ftype == FILE_TYPE::UGRID)
    {
    }
    if (m_ftype == FILE_TYPE::SGRID ||
        m_ftype == FILE_TYPE::KISS)
    {
        START_TIMERN(Read_sgrid_mesh);
        status = this->read_sgrid_mesh();
        STOP_TIMER(Read_sgrid_mesh);

        //START_TIMER(Read_times);
        //status = this->read_times();
        //STOP_TIMER(Read_times);

        START_TIMER(Read_sgrid_variables);
        status = this->read_sgrid_variables();
        STOP_TIMER(Read_sgrid_variables);
    }
    return status;
}
//------------------------------------------------------------------------------
long SGRID::read_global_attributes()
{
    long status;
    int ndims, nvars, natts, nunlimited;
    nc_type att_type;
    size_t att_length;

#ifdef NATIVE_C
    fprintf(stderr, "SGRID::read_global_attributes()\n");
#endif    

    char * att_name_c = (char *)malloc(sizeof(char) * (NC_MAX_NAME + 1));
    att_name_c[0] = '\0';

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);

    this->global_attributes = (struct _global_attributes *)malloc(sizeof(struct _global_attributes));
    this->global_attributes->count = natts;
    this->global_attributes->attribute = (struct _global_attribute **)malloc(sizeof(struct _global_attribute *) * natts);
    for (long i = 0; i < natts; i++)
    {
        this->global_attributes->attribute[i] = new _global_attribute;
    }

    for (long i = 0; i < natts; i++)
    {
        status = nc_inq_attname(this->m_ncid, NC_GLOBAL, i, att_name_c);
        status = nc_inq_att(this->m_ncid, NC_GLOBAL, att_name_c, &att_type, &att_length);
        this->global_attributes->attribute[i]->name = std::string(strdup(att_name_c));
        this->global_attributes->attribute[i]->type = att_type;
        this->global_attributes->attribute[i]->length = att_length;
        if (att_type == NC_CHAR)
        {
            char * att_value_c = (char *)malloc(sizeof(char) * (att_length + 1));
            att_value_c[0] = '\0';
            status = nc_get_att_text(this->m_ncid, NC_GLOBAL, att_name_c, att_value_c);
            att_value_c[att_length] = '\0';
            this->global_attributes->attribute[i]->cvalue = std::string(strdup(att_value_c));
            free(att_value_c);
            att_value_c = nullptr;
        }
        else
        {
#ifdef NATIVE_C
            fprintf(stderr, "    Attribute nc_type: %d\n", att_type);
#endif    
        }
    }
    free(att_name_c);
    att_name_c = nullptr;

    return status;
}
//------------------------------------------------------------------------------
long SGRID::get_global_attribute_value(std::string att_name, std::string* att_value)
{
    long status = 1;  // error by default
    bool att_name_found = false;
    for (int i = 0; i < global_attributes->count; i++)
    {
        if (att_name_found) { exit; }
        if (global_attributes->attribute[i]->name == att_name)
        {
            att_name_found = true;
            *att_value = global_attributes->attribute[i]->cvalue;
            *att_value = std::string(global_attributes->attribute[i]->cvalue);
            status = 0;
        }
    }
    return status;
}

//------------------------------------------------------------------------------
FILE_TYPE SGRID::get_file_type()
{
    return m_ftype;
}
//------------------------------------------------------------------------------
std::string SGRID::get_filename()
{
    return this->m_grid_file_name;
}
//------------------------------------------------------------------------------
long SGRID::read_sgrid_mesh()
{
    long status = 1;
    int ndims, nvars, natts, nunlimited;
    int unlimid;
    std::string cf_role;
    std::string grid_mapping_name;
    std::string std_name;
    std::vector<std::string> tmp_dim_names;
    size_t length;

    char* var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));
    m_mapping = new _mapping();
    set_grid_mapping_epsg(0, "EPSG:0");

    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);
    m_dimids = (size_t*)malloc(sizeof(size_t) * ndims);
    for (int i = 0; i < ndims; i++)
    {
        status = nc_inq_dimlen(this->m_ncid, i, &m_dimids[i]);
        status = nc_inq_dimname(this->m_ncid, i, var_name_c);
        m_dim_names.push_back(std::string(var_name_c));
        m_map_dim[std::string(var_name_c)] = (long)m_dimids[i];
    }
    status = nc_inq_unlimdim(this->m_ncid, &unlimid);
    for (int i_var = 0; i_var < nvars; i_var++)
    {
        nc_type nc_type;
        nc_inq_varndims(this->m_ncid, i_var, &ndims);
        int* var_dimids = NULL;
        if (ndims > 0) {
            var_dimids = (int*)malloc(sizeof(int) * ndims);
        }
        status = nc_inq_var(this->m_ncid, i_var, var_name_c, &nc_type, &ndims, var_dimids, &natts);
        std::string var_name(var_name_c);
#ifdef NATIVE_C
        fprintf(stderr, "SGRID::read_mesh()\n    Variable name: %d - %s\n", i_var + 1, var_name.c_str());
#endif    

        length = 0;
        status = get_attribute(this->m_ncid, i_var, "cf_role", &cf_role);
        if (status == NC_NOERR)
        {
            status = read_variables_with_cf_role(i_var, var_name, cf_role, ndims, var_dimids);
        }
        else if (get_file_type() == FILE_TYPE::KISS)
        {
            // no grid-variable defined, so presume var_name="grid2d"and cf_role="grid_topology"
            status = read_variables_with_cf_role(i_var, "grid2d", "grid_topology", ndims, var_dimids);
        }

        status = get_attribute(this->m_ncid, i_var, "standard_name", &std_name);
        if (std_name == "ocean_sigma_coordinate")
        {
            // this grid contains sigma layers, find sigma dimension name
            tmp_dim_names = get_dimension_names(this->m_ncid, var_name);
            for (int i = 0; i < tmp_dim_names.size(); i++)
            {
                if (tmp_dim_names[i].find("nterface") != std::string::npos)  // Todo: HACK: dimension name should have the sub-string 'interface' or 'Interface'
                {
                    m_map_dim_name["zs_dim_interface"] = tmp_dim_names[i];
                    m_map_dim_name["zs_name_interface"] = var_name;
                }
                else
                {
                    m_map_dim_name["zs_dim_layer"] = tmp_dim_names[i];
                    m_map_dim_name["zs_name_layer"] = var_name;
                }
            }
        }
        else if (std_name == "altitude")
        {
            // this grid contains z- layers, find z-dimension name
            tmp_dim_names = get_dimension_names(this->m_ncid, var_name);
            for (int i = 0; i < tmp_dim_names.size(); i++)
            {
                if (tmp_dim_names[i].find("nterface") != std::string::npos)
                {
                    m_map_dim_name["zs_dim_interface"] = tmp_dim_names[i];
                    m_map_dim_name["zs_name_interface"] = var_name;
                }
                else if (tmp_dim_names[i].find("Layer") != std::string::npos)
                {
                    m_map_dim_name["zs_dim_layer"] = tmp_dim_names[i];
                    m_map_dim_name["zs_name_layer"] = var_name;
                }
            }
        }
        else
        {
            tmp_dim_names = get_dimension_names(this->m_ncid, var_name);
            for (int i = 0; i < tmp_dim_names.size(); i++)
            {
                if (tmp_dim_names[i].find("nBedLayers") != std::string::npos)
                {
                    m_map_dim_name["zs_dim_bed_layer"] = tmp_dim_names[i];
                    m_map_dim_name["zs_name_bed_layer"] = var_name;
                }
            }
        }

        status = get_attribute(this->m_ncid, i_var, "grid_mapping_name", &grid_mapping_name);
        if (status == NC_NOERR)
        {
            status = read_grid_mapping(i_var, var_name, grid_mapping_name);
        }
        free(var_dimids);
        var_dimids = nullptr;
    }

    free(var_name_c);
    var_name_c = nullptr;

    return status;
}
//------------------------------------------------------------------------------
long SGRID::set_grid_mapping_epsg(long epsg, std::string epsg_code)
{
    long status = 0;
    m_mapping->epsg = epsg;
    m_mapping->epsg_code = epsg_code;
    return status;
}
//------------------------------------------------------------------------------
long SGRID::read_sgrid_variables()
{
    // read the attributes of the variable  and the dimensions
#ifdef NATIVE_C
    fprintf(stderr, "SGRID::read_variables()\n");
#endif    
    int ndims, nvars, natts, nunlimited;
    int status;
    struct _variable* var;
    var = NULL;
    m_nr_mesh_var = 0;

    std::string tmp_string;
    char* var_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));
    status = nc_inq(this->m_ncid, &ndims, &nvars, &natts, &nunlimited);

    for (int i_var = 0; i_var < nvars; i_var++)
    {
        nc_type nc_type;
        nc_inq_varndims(this->m_ncid, i_var, &ndims);
        int* var_dimids = NULL;
        if (ndims > 0) {
            var_dimids = (int*)malloc(sizeof(int) * ndims);
        }
        status = nc_inq_var(this->m_ncid, i_var, var_name_c, &nc_type, &ndims, var_dimids, &natts);
        std::string var_name(var_name_c);
        status = get_attribute(this->m_ncid, i_var, "grid", &tmp_string);  // statement, to detect if this is a variable on a structured grid
        if (status == NC_NOERR) { status = get_attribute(this->m_ncid, i_var, "location", &tmp_string); } // each variable does have a location
        if (status != NC_NOERR && m_ftype == FILE_TYPE::KISS) {
            tmp_string = "node";
            status = NC_NOERR;
        }
        if (status == NC_NOERR)
        {
            // This variable is defined on a mesh and has a dimension
            m_nr_mesh_var += 1;
            if (m_nr_mesh_var == 1)
            {
                m_mesh_vars = new _mesh_variable();
                m_mesh_vars->variable = (struct _variable**)malloc(sizeof(struct _variable*));
            }
            else
            {
                m_mesh_vars->variable = (struct _variable**)realloc(m_mesh_vars->variable, sizeof(struct _variable*) * m_nr_mesh_var);
            }
            m_mesh_vars->variable[m_nr_mesh_var - 1] = new _variable();
            m_mesh_vars->nr_vars = m_nr_mesh_var;

            m_mesh_vars->variable[m_nr_mesh_var - 1]->var_name = var_name;
            m_mesh_vars->variable[m_nr_mesh_var - 1]->nc_type = nc_type;
            m_mesh_vars->variable[m_nr_mesh_var - 1]->read = false;
            status = get_attribute(this->m_ncid, i_var, "location", &m_mesh_vars->variable[m_nr_mesh_var - 1]->location);
            if (status != NC_NOERR && m_ftype == FILE_TYPE::KISS)
            {
                m_mesh_vars->variable[m_nr_mesh_var - 1]->location = "node";
                status = NC_NOERR; 
            }
            status = get_attribute(this->m_ncid, i_var, "grid", &m_mesh_vars->variable[m_nr_mesh_var - 1]->mesh);
            if (status != NC_NOERR && m_ftype == FILE_TYPE::KISS)
            {
                m_mesh_vars->variable[m_nr_mesh_var - 1]->mesh = "grid2d";
                status = NC_NOERR;
            }
            status = get_attribute(this->m_ncid, i_var, "coordinates", &m_mesh_vars->variable[m_nr_mesh_var - 1]->coordinates);
            status = get_attribute(this->m_ncid, i_var, "cell_methods", &m_mesh_vars->variable[m_nr_mesh_var - 1]->cell_methods);
            status = get_attribute(this->m_ncid, i_var, "standard_name", &m_mesh_vars->variable[m_nr_mesh_var - 1]->standard_name);
            status = get_attribute(this->m_ncid, i_var, "long_name", &m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name);
            if (status != NC_NOERR || m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name.size() <= 1)
            {
                m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name = m_mesh_vars->variable[m_nr_mesh_var - 1]->standard_name;
                if (m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name.size() <= 1)
                {
                    m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name = m_mesh_vars->variable[m_nr_mesh_var - 1]->var_name;
                }
            }
            status = get_attribute(this->m_ncid, i_var, "units", &m_mesh_vars->variable[m_nr_mesh_var - 1]->units);
            status = get_attribute(this->m_ncid, i_var, "grid_mapping", &m_mesh_vars->variable[m_nr_mesh_var - 1]->grid_mapping);
            status = get_attribute(this->m_ncid, i_var, const_cast<char*>("_FillValue"), &m_mesh_vars->variable[m_nr_mesh_var - 1]->fill_value);
            status = get_attribute(this->m_ncid, i_var, "comment", &m_mesh_vars->variable[m_nr_mesh_var - 1]->comment);  // UGRID
            status = get_attribute(this->m_ncid, i_var, "description", &m_mesh_vars->variable[m_nr_mesh_var - 1]->description);  // SGRID

            // variable with flag_values and flag_meanings
            std::string flag_meanings;
            status = get_attribute(this->m_ncid, i_var, "flag_meanings", &flag_meanings);
            if (status == NC_NOERR)
            {
                size_t length;
                status = nc_inq_attlen(this->m_ncid, i_var, "flag_values", &length);
                int* flag_values_c = (int*)malloc(sizeof(int) * length);
                status = get_attribute(this->m_ncid, i_var, const_cast<char*>("flag_values"), flag_values_c);
                for (int j = 0; j < length; ++j)
                {
                    m_mesh_vars->variable[m_nr_mesh_var - 1]->flag_values.push_back(flag_values_c[j]);
                }
                m_mesh_vars->variable[m_nr_mesh_var - 1]->flag_meanings = tokenize(flag_meanings, ' ');
            }

            m_mesh_vars->variable[m_nr_mesh_var - 1]->draw = false;
            m_mesh_vars->variable[m_nr_mesh_var - 1]->time_series = false;
            for (int j = 0; j < ndims; j++)
            {
                status = nc_inq_dimname(this->m_ncid, var_dimids[j], var_name_c);

                m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.push_back((long)m_dimids[var_dimids[j]]);  // Todo: HACK typecast: size_t -> long
                m_mesh_vars->variable[m_nr_mesh_var - 1]->dim_names.push_back(m_dim_names[var_dimids[j]]);
            }

            int topo_dim;
            int mesh_id;
            status = nc_inq_varid(this->m_ncid, m_mesh_vars->variable[m_nr_mesh_var - 1]->mesh.c_str(), &mesh_id);
            status = get_attribute(this->m_ncid, mesh_id, const_cast<char*>("topology_dimension"), &topo_dim);
            m_mesh_vars->variable[m_nr_mesh_var - 1]->topology_dimension = topo_dim;
            if (m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 1)
            {
                // nothing to do
            }
            if (m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 2 ||
                m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 3 && m_ftype == FILE_TYPE::SGRID ||
                m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 3 && m_ftype == FILE_TYPE::KISS)
            {
                bool contains_time_dimension = false;
                for (int i = 0; i < m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size(); i++)
                {
                    // check if one of the dimension is the time dimension
                }
                if (contains_time_dimension)
                {
                    // nothing special: 2D (time, xy-space)
                }
                else
                {
                    // special: 2D (nSedTot, xy-space)
                    for (int i = 0; i < m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size(); i++)
                    {
                        if (m_mesh_vars->variable[m_nr_mesh_var - 1]->dim_names[i] == "nSedTot")  // todo: hard coded sediment dimension
                        {
                            std::vector<std::string> sed_name = get_string_var(this->m_ncid, "sedfrac_name");
                            int nr_sedtot = m_mesh_vars->variable[m_nr_mesh_var - 1]->dims[i];
                            for (int j = 0; j < nr_sedtot; j++)
                            {
                                if (j > 0)
                                {
                                    m_nr_mesh_var += 1;
                                    m_mesh_vars->variable = (struct _variable**)realloc(m_mesh_vars->variable, sizeof(struct _variable*) * m_nr_mesh_var);
                                    m_mesh_vars->variable[m_nr_mesh_var - 1] = new _variable();
                                    m_mesh_vars->nr_vars = m_nr_mesh_var;
                                }
                                std::stringstream ss;
                                std::streamsize nsig = int(log10(nr_sedtot)) + 1;
                                ss << std::setfill('0') << std::setw(nsig) << j + 1;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->draw = false;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->time_series = false;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->mesh = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->mesh;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->location = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->location;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->var_name = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->var_name;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->coordinates = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->coordinates;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->topology_dimension = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->topology_dimension;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->nc_type = m_mesh_vars->variable[m_nr_mesh_var - 1 - j]->nc_type;
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->sediment_index = j;
                                std::string janm = var_name + "- " + sed_name[j] + " - " + "nSedTot";
                                m_mesh_vars->variable[m_nr_mesh_var - 1]->long_name = strdup(janm.c_str());
                            }
                        }
                    }
                }
            }
            else if (m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 3 )
            {
                bool contains_time_dimension = false;
                for (int i = 0; i < m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size(); i++)
                {
                    // check if one of the dimension is the time dimension
                }
            }
            else if (m_mesh_vars->variable[m_nr_mesh_var - 1]->dims.size() == 4)
            {
                // 4D variable:
                // - time
                // - sediment dimension
                // - xy-space dimension
                // - z-space (nbedlayers)
            }
        }
    }

    free(var_name_c);
    var_name_c = nullptr;

    return (long)status;
}

//------------------------------------------------------------------------------
struct _mesh_variable * SGRID::get_variables()
{
   return this->m_mesh_vars;
}
//------------------------------------------------------------------------------
struct _mesh2d * SGRID::get_mesh_2d()
{
    return m_mesh2d;
}
//------------------------------------------------------------------------------
struct _mesh2d_string ** SGRID::get_mesh2d_string()
{
#ifdef NATIVE_C
    fprintf(stderr, "SGRID::get_mesh2d()\n");
#endif    
    return m_mesh2d_strings;
}
//------------------------------------------------------------------------------
struct _mapping * SGRID::get_grid_mapping()
{
#ifdef NATIVE_C
    fprintf(stderr, "SGRID::get_grid_mapping()\n");
#endif    
    return m_mapping;
}
//------------------------------------------------------------------------------
struct _variable * SGRID::get_var_by_std_name(struct _mesh_variable * vars, std::string mesh, std::string standard_name)
{
    for (int i = 0; i < vars->nr_vars; i++)
    {
        if (standard_name == vars->variable[i]->standard_name && mesh == vars->variable[i]->mesh)
        {
            return vars->variable[i];
        }
    }
    return nullptr;
}
//==============================================================================
// PRIVATE functions
//==============================================================================
int SGRID::read_mesh2d_attributes(struct _mesh2d_string* mesh2d_strings, int i_var, std::string var_name, int topology_dimension)
{
    int status = 1;

    mesh2d_strings->var_name = var_name;

    // Required attributes
    //
    // cf_role == mesh_topology
    // topology_dimension
    // node coordinates
    // face_node_connectivity
    mesh2d_strings->topology_dimension = size_t(topology_dimension);
    status = get_attribute(this->m_ncid, i_var, "node_coordinates", &mesh2d_strings->node_coordinates);
    if (m_ftype == FILE_TYPE::UGRID)
    {
        if (status == NC_NOERR) { status = get_attribute(this->m_ncid, i_var, "face_node_connectivity", &mesh2d_strings->face_node_connectivity); }
        if (status != NC_NOERR)
        {
            fprintf(stderr, "    Mesh \'%s\' does not meet the UGRID standard for 2D meshes. Required attributes are missing.\n", mesh2d_strings->var_name.c_str());
        }
    }

    // optional required attributes
    status = get_attribute(this->m_ncid, i_var, "edge_dimension", &mesh2d_strings->edge_dimension);
    status = get_attribute(this->m_ncid, i_var, "face_dimension", &mesh2d_strings->face_dimension);
    status = get_attribute(this->m_ncid, i_var, "edge_node_connectivity", &mesh2d_strings->edge_node_connectivity);

    // optional attributes
    //status = get_attribute(this->m_ncid, i_var, "coordinate_space", &mesh2d_strings[nr_mesh2d - 1]->coordinate_space);
    status = get_attribute(this->m_ncid, i_var, "edge_coordinates", &mesh2d_strings->edge_coordinates);
    status = get_attribute(this->m_ncid, i_var, "edge_geometry", &mesh2d_strings->edge_geometry);
    status = get_attribute(this->m_ncid, i_var, "edge_length", &mesh2d_strings->edge_length);
    status = get_attribute(this->m_ncid, i_var, "long_name", &mesh2d_strings->long_name);
    //status = get_attribute(this->m_ncid, i_var, "node_dimension", &mesh2d_strings->node_dimension);
    status = get_attribute(this->m_ncid, i_var, "node_edge_exchange", &mesh2d_strings->node_edge_exchange);

    status = get_attribute(this->m_ncid, i_var, "max_face_nodes_dimension", &mesh2d_strings->max_face_nodes_dimension);
    status = get_attribute(this->m_ncid, i_var, "edge_face_connectivity", &mesh2d_strings->edge_face_connectivity);
    status = get_attribute(this->m_ncid, i_var, "face_coordinates", &mesh2d_strings->face_coordinates);
    status = get_attribute(this->m_ncid, i_var, "edge_type", &mesh2d_strings->edge_type);
    if (status != NC_NOERR)
    {
        mesh2d_strings->edge_type = var_name + "_edge_type";
    }

    status = get_attribute(this->m_ncid, i_var, "layer_dimension", &mesh2d_strings->layer_dimension);
    if (status == NC_NOERR)
    {
        m_map_dim_name["zs_dim_layer"] = mesh2d_strings->layer_dimension;
        status = get_attribute(this->m_ncid, i_var, "interface_dimension", &mesh2d_strings->layer_interface_dimension);
        if (status == NC_NOERR)
        {
            m_map_dim_name["zs_dim_interface"] = mesh2d_strings->layer_interface_dimension;
        }
    }
    if (status != NC_NOERR)
    {
        mesh2d_strings->layer_dimension = "";
        mesh2d_strings->layer_interface_dimension = "";
    }
    // split required 'node coordinate' string
    std::vector<std::string> token = tokenize(mesh2d_strings->node_coordinates, ' ');
    if (token.size() == 2)
    {
        mesh2d_strings->x_node_name = token[0];
        mesh2d_strings->y_node_name = token[1];
    }

    // split 'edge coordinate' string
    token = tokenize(mesh2d_strings->edge_coordinates, ' ');
    if (token.size() == 2)
    {
        mesh2d_strings->x_edge_name = token[0];
        mesh2d_strings->y_edge_name = token[1];
        status = get_attribute_by_var_name(this->m_ncid, mesh2d_strings->x_edge_name, "bounds", &mesh2d_strings->x_bound_edge_name);
        status = get_attribute_by_var_name(this->m_ncid, mesh2d_strings->y_edge_name, "bounds", &mesh2d_strings->y_bound_edge_name);
    }

    // split 'face_coordinates' string
    token = tokenize(mesh2d_strings->face_coordinates, ' ');
    if (token.size() == 2)
    {
        mesh2d_strings->x_face_name = token[0];
        mesh2d_strings->y_face_name = token[1];
        status = get_attribute_by_var_name(this->m_ncid, mesh2d_strings->x_face_name, "bounds", &mesh2d_strings->x_bound_face_name);
        status = get_attribute_by_var_name(this->m_ncid, mesh2d_strings->y_face_name, "bounds", &mesh2d_strings->y_bound_face_name);
    }

    return status;
}
//------------------------------------------------------------------------------
std::vector<std::string> SGRID::get_names(int ncid, std::string names, size_t count)
{
    int var_id;
    int status;
    std::vector<std::string> token;

    status = nc_inq_varid(ncid, names.c_str(), &var_id);
    if (status == NC_NOERR)
    {
        int ndims[2];
        char* length_name = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));;
        status = nc_inq_vardimid(ncid, var_id, (int*)ndims);
        status = nc_inq_dimname(ncid, ndims[1], length_name);
        size_t strlen;
        status = get_dimension(ncid, length_name, &strlen);

        char* c = (char*)malloc(sizeof(char) * (count * strlen));
        status = nc_get_var_text(ncid, var_id, c);

        token = tokenize(c, count);
        free(c);
        c = nullptr;
        free(length_name);
        length_name = nullptr;
    }
    return token;
}
//------------------------------------------------------------------------------
std::vector<std::string> SGRID::get_string_var(int ncid, std::string var_name)
{
    int varid;
    int ndims;
    int nr_names = 0;
    nc_type nc_type;
    int status = -1;
    std::vector<std::string> result;

    char* dim_name_c = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));
    dim_name_c[0] = '\0';

    if (var_name.size() != 0)
    {
        status = nc_inq_varid(ncid, var_name.c_str(), &varid);
        status = nc_inq_var(ncid, varid, NULL, &nc_type, &ndims, NULL, NULL);
        long* sn_dims = (long*)malloc(sizeof(long) * ndims);
        status = nc_inq_vardimid(ncid, varid, (int*)sn_dims);

        long mem_length = 1;
        long name_len = -1;
        for (long j = 0; j < ndims; j++)
        {
            size_t length = (size_t)-1;
            status = nc_inq_dim(this->m_ncid, sn_dims[j], dim_name_c, &length);

            if (strstr(dim_name_c, "len") || name_len == -1 && j == 1)  // second dimension is the string length if not already set
            {
                name_len = (long)length;
            }
            else
            {
                nr_names = (long)length;
            }
            mem_length = mem_length * (long)length;
        }
        result.reserve(nr_names);
        // reading 64 strings for each location, length of string??
        if (nc_type == NC_STRING)
        {
            char** location_strings = (char**)malloc(sizeof(char*) * (mem_length)+1);
            status = nc_get_var_string(ncid, varid, location_strings);
            //QString janm = QString("JanM");
            //for (int k = 0; k < nr_names; k++)
            //{
            //    janm = QString("");
            //    for (int k2 = 0; k2 < name_len; k2++)
            //    {
            //        QTextCodec* codec2 = QTextCodec::codecForName("UTF-8");
            //        QString str = codec2->toUnicode(*(location_strings + k * name_len + k2));
            //        janm = janm + str;
            //    }
            //    result.push_back(janm.toStdString());
            //}
        }
        else if (nc_type == NC_CHAR)
        {
            char* location_chars = (char*)malloc(sizeof(char*) * (mem_length)+1);
            location_chars[0] = '\0';
            status = nc_get_var_text(ncid, varid, location_chars);
            char* janm = (char*)malloc(sizeof(char) * (name_len + 1));
            janm[name_len] = '\0';
            for (int k = 0; k < nr_names; k++)
            {
                strncpy(janm, location_chars + k * name_len, name_len);
                result.push_back(std::string(janm));
            }
            free(janm);
        }
        else
        {
            // trying to read unsupported variable
        }
        free(sn_dims);
    }
    free(dim_name_c);
    return result;
}
//------------------------------------------------------------------------------
int SGRID::read_variables_with_cf_role(int i_var, std::string var_name, std::string cf_role, int ndims, int* var_dimids)
{
    int topology_dimension;
    int status = 1;
    int var_id = -1;

    long nr_ntw = 0;
    long nr_geom = 0;
    long nr_mesh1d = 0;
    long nr_mesh2d = 0;

    std::string att_value;
    double fill_value;

    if (cf_role == "grid_topology")
    {
        topology_dimension = 0;
        status = get_attribute(this->m_ncid, i_var, const_cast<char*>("topology_dimension"), &topology_dimension);
        if (status != NC_NOERR && m_ftype == FILE_TYPE::KISS)
        {
            topology_dimension = 2;
        }
        ///////////////////////////////////////////////////////////////////////////////////////////
        if (topology_dimension == 2)  // it is a unstructured mesh
        {
            nr_mesh2d += 1;
            if (nr_mesh2d == 1)
            {
                m_mesh2d_strings = (struct _mesh2d_string**)malloc(sizeof(struct _mesh2d_string*) * nr_mesh2d);
            }
            else
            {
            }
            m_mesh2d_strings[nr_mesh2d - 1] = new _mesh2d_string;

            status = read_mesh2d_attributes(m_mesh2d_strings[nr_mesh2d - 1], i_var, var_name, topology_dimension);
            if (status != NC_NOERR && m_ftype == FILE_TYPE::KISS)
            {
                m_mesh2d_strings[nr_mesh2d - 1]->x_node_name = "x_coordinate";
                m_mesh2d_strings[nr_mesh2d - 1]->y_node_name = "y_coordinate";
            }
            if (nr_mesh2d == 1)
            {
                m_mesh2d = new _mesh2d();
                m_mesh2d->node = (struct _feature**)malloc(sizeof(struct _feature*));
                m_mesh2d->node[nr_mesh2d - 1] = new _feature();

                m_mesh2d->edge = (struct _edge**)malloc(sizeof(struct _edge*));
                m_mesh2d->edge[nr_mesh2d - 1] = new _edge();

                m_mesh2d->face = (struct _feature**)malloc(sizeof(struct _feature*));
                m_mesh2d->face[nr_mesh2d - 1] = new _feature();
            }
            else
            {
            }
            m_mesh2d->nr_mesh2d = nr_mesh2d;

#ifdef NATIVE_C
            fprintf(stderr, "    Variables with \'mesh_topology\' attribute: %s\n", var_name.c_str());
#endif    

            //get edge nodes, optional required
            int* dimids;

            status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_node_name.c_str(), &var_id);
            status = nc_inq_varndims(this->m_ncid, var_id, &ndims);
            dimids = (int*)malloc(sizeof(int) * ndims);
            status = nc_inq_vardimid(this->m_ncid, var_id, dimids);
            int imax_node = m_dimids[dimids[0]];
            int jmax_node = 0;
            if (m_ftype == FILE_TYPE::KISS)
            {
                status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_node_name.c_str(), &var_id);
                status = nc_inq_varndims(this->m_ncid, var_id, &ndims);
                dimids = (int*)malloc(sizeof(int) * ndims);
                status = nc_inq_vardimid(this->m_ncid, var_id, dimids);
                jmax_node = (int)m_dimids[dimids[0]];
            }
            else
            {
                jmax_node = (int)m_dimids[dimids[1]];
            }

            m_mesh2d->edge[nr_mesh2d - 1]->count = imax_node * (jmax_node - 1) + jmax_node * (imax_node - 1);
            mesh2d_edge_nodes = (int*)malloc(sizeof(int) * m_mesh2d->edge[nr_mesh2d - 1]->count * _two);
            m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes = (int**)malloc(sizeof(int*) * m_mesh2d->edge[nr_mesh2d - 1]->count);
            for (int i = 0; i < m_mesh2d->edge[nr_mesh2d - 1]->count; i++)
            {
                m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes[i] = mesh2d_edge_nodes + _two * i;
            }
            // vertical edges
            int k = -1;
            for (int i = 0; i <imax_node; ++i)
            {
                for (int j = 0; j < jmax_node-1; ++j)
                {
                    k += 1;
                    m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes[k][0] = i * jmax_node + j;
                    m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes[k][1] = i * jmax_node + j + 1;
                }
            }
            //horizontal edges
            for (int j = 0; j < jmax_node; ++j)
            {
                for (int i = 0; i < imax_node - 1; ++i)
                {
                    k += 1;
                    m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes[k][0] = i * jmax_node + j;
                    m_mesh2d->edge[nr_mesh2d - 1]->edge_nodes[k][1] = (i + 1) * jmax_node + j;
                }
            }


            /* Read the data (x, y)-coordinate of each node */
            size_t length = 1;
            status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_node_name.c_str(), &var_id);
            status = nc_inq_varndims(this->m_ncid, var_id, &ndims);
            dimids = (int*)malloc(sizeof(int) * ndims);
            status = nc_inq_vardimid(this->m_ncid, var_id, dimids);
            for (int i = 0; i < ndims; i++)
            {
                length *= m_dimids[dimids[i]];
                m_mesh2d->node[nr_mesh2d - 1]->dims.push_back(m_dimids[dimids[i]]);
            }
            m_mesh2d->node[nr_mesh2d - 1]->count = length;
            m_mesh2d->node[nr_mesh2d - 1]->x = std::vector<double>(m_mesh2d->node[nr_mesh2d - 1]->count);
            m_mesh2d->node[nr_mesh2d - 1]->y = std::vector<double>(m_mesh2d->node[nr_mesh2d - 1]->count);
            status = get_attribute(this->m_ncid, var_id, "standard_name", &att_value);
            status = get_attribute(this->m_ncid, var_id, strdup("_FillValue"), &m_mesh2d->node[nr_mesh2d - 1]->fill_value);
            if (att_value == "projection_x_coordinate" || att_value == "longitude")
            {
                status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->node[nr_mesh2d - 1]->x.data());
                status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_node_name.c_str(), &var_id);
                status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->node[nr_mesh2d - 1]->y.data());
            }
            else
            {
                status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->node[nr_mesh2d - 1]->y.data());
                status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_node_name.c_str(), &var_id);
                status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->node[nr_mesh2d - 1]->x.data());
            }



            /* Read the data (x, y)-coordinate of each edge */

            /* Read the data (x, y)-coordinate of each face */
            if (m_mesh2d_strings[nr_mesh2d - 1]->x_face_name != "")
            {
                status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_face_name.c_str(), &var_id);
                status = nc_inq_varndims(this->m_ncid, var_id, &ndims);
                dimids = (int*)malloc(sizeof(int) * ndims);
                status = nc_inq_vardimid(this->m_ncid, var_id, dimids);
                length = 1;
                for (int i = 0; i < ndims; i++)
                {
                    length *= m_dimids[dimids[i]];
                    m_mesh2d->face[nr_mesh2d - 1]->dims.push_back(m_dimids[dimids[i]]);
                }

                m_mesh2d->face[nr_mesh2d - 1]->count = length;
                m_mesh2d->face[nr_mesh2d - 1]->x = std::vector<double>(m_mesh2d->face[nr_mesh2d - 1]->count);
                m_mesh2d->face[nr_mesh2d - 1]->y = std::vector<double>(m_mesh2d->face[nr_mesh2d - 1]->count);

                if (m_mesh2d->face[nr_mesh2d - 1]->count != 0)  // not required attribute
                {
                    m_mesh2d->face[nr_mesh2d - 1]->x = std::vector<double>(m_mesh2d->face[nr_mesh2d - 1]->count);
                    m_mesh2d->face[nr_mesh2d - 1]->y = std::vector<double>(m_mesh2d->face[nr_mesh2d - 1]->count);
                    status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_face_name.c_str(), &var_id);
                    status = get_attribute(this->m_ncid, var_id, "standard_name", &att_value);
                    status = get_attribute(this->m_ncid, var_id, strdup("_FillValue"), &m_mesh2d->face[nr_mesh2d - 1]->fill_value);
                    if (att_value == "projection_x_coordinate" || att_value == "longitude")
                    {
                        status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->face[nr_mesh2d - 1]->x.data());
                        status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_face_name.c_str(), &var_id);
                        status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->face[nr_mesh2d - 1]->y.data());
                    }
                    else
                    {
                        status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->face[nr_mesh2d - 1]->y.data());
                        status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_face_name.c_str(), &var_id);
                        status = nc_get_var_double(this->m_ncid, var_id, m_mesh2d->face[nr_mesh2d - 1]->x.data());
                    }
                    status = get_attribute_by_var_name(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_face_name, "bounds", &m_mesh2d_strings[nr_mesh2d - 1]->x_bound_face_name);
                    status = get_attribute_by_var_name(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->y_face_name, "bounds", &m_mesh2d_strings[nr_mesh2d - 1]->y_bound_face_name);
                }
            }
            else
            {
                // Compute the face coordinates from the node coordinates, mass centre of face. But is it needed?
                // Boundary of the face is determined by the node coordinates (face_nodes)
            }

            /* Read the nodes indices for each face */
            /* Determine, because they are not on the file */
            status = nc_inq_varid(this->m_ncid, m_mesh2d_strings[nr_mesh2d - 1]->x_face_name.c_str(), &var_id);
            status = nc_inq_varndims(this->m_ncid, var_id, &ndims);
            dimids = (int*)malloc(sizeof(int) * ndims);
            status = nc_inq_vardimid(this->m_ncid, var_id, dimids);

            int m_max;
            int n_max;
            if (m_ftype == FILE_TYPE::KISS)
            {
                m_max = imax_node;
                n_max = jmax_node;
            }
            else
            {
                m_max = m_dimids[dimids[0]] + 1;  // HACK: 1 more node then faces, only true for structured grids
                n_max = m_dimids[dimids[1]] + 1;  // HACK: 1 more node then faces, only true for structured grids
            }

            std::vector<int> value;
            for (int m = 0; m < m_max-1; m++)  // faces x-direction
            {
                for (int n = 0; n < n_max-1; n++)  // faces y-direction
                {
                    value.push_back(m * n_max + n);
                    value.push_back((m+1) * n_max + n);
                    value.push_back((m+1) * n_max + n + 1);
                    value.push_back(m * n_max + n + 1);

                    m_mesh2d->face_nodes.push_back(value);
                    value.clear();
                }
            }


            // Do not name the mesh2d nodes, because the 2D mesh can be very large
            length = m_mesh2d->node[nr_mesh2d - 1]->name.size();
            length += m_mesh2d->node[nr_mesh2d - 1]->long_name.size();
            length += m_mesh2d->edge[nr_mesh2d - 1]->name.size();
            length += m_mesh2d->edge[nr_mesh2d - 1]->long_name.size();
            length += m_mesh2d->face[nr_mesh2d - 1]->name.size();
            length += m_mesh2d->face[nr_mesh2d - 1]->long_name.size();
            if (length != 0)
            {
                fprintf(stderr, "Length of the node/edge names of a 2D mesh have to be zero");
                return 1;
            }
        }
    }

    return status;
}

int SGRID::get_attribute_by_var_name(int ncid, std::string var_name, std::string att_name, std::string* att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(ncid, var_name.c_str(), &i_var);
    status = get_attribute(ncid, i_var, att_name, att_value);
    return status;
}
int SGRID::read_grid_mapping(int i_var, std::string var_name, std::string grid_mapping_name)
{
    int status = -1;
    m_mapping->epsg = -1;

    status = get_attribute(this->m_ncid, i_var, "name", &m_mapping->name);
    status = get_attribute(this->m_ncid, i_var, const_cast<char*>("epsg"), &m_mapping->epsg);
    m_mapping->grid_mapping_name = grid_mapping_name; //  == status = get_attribute(this->m_ncid, i_var, "grid_mapping_name", &map->grid_mapping_name);
    status = get_attribute(this->m_ncid, i_var, const_cast<char*>("longitude_of_prime_meridian"), &m_mapping->longitude_of_prime_meridian);
    status = get_attribute(this->m_ncid, i_var, const_cast<char*>("semi_major_axis"), &m_mapping->semi_major_axis);
    status = get_attribute(this->m_ncid, i_var, const_cast<char*>("semi_minor_axis"), &m_mapping->semi_minor_axis);
    status = get_attribute(this->m_ncid, i_var, const_cast<char*>("inverse_flattening"), &m_mapping->inverse_flattening);
    status = get_attribute(this->m_ncid, i_var, "EPSG_code", &m_mapping->epsg_code);
    status = get_attribute(this->m_ncid, i_var, "value", &m_mapping->value);
    status = get_attribute(this->m_ncid, i_var, "projection_name", &m_mapping->projection_name);
    status = get_attribute(this->m_ncid, i_var, "wkt", &m_mapping->wkt);

    if (m_mapping->epsg == -1 && m_mapping->epsg_code.size() != 0)
    {
        std::vector<std::string> token = tokenize(m_mapping->epsg_code, ':');
        m_mapping->epsg = atoi(token[1].c_str());  // second token contains the plain EPSG code
    }

    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, char * att_name, char ** att_value)
{
    size_t length = 0;
    int status = -1;

    status = nc_inq_attlen(ncid, i_var, att_name, &length);
    *att_value = (char *)malloc(sizeof(char) * (length + 1));
    *att_value[0] = '\0';
    if (status == NC_NOERR)
    {
        status = nc_get_att(ncid, i_var, att_name, *att_value);
        att_value[0][length] = '\0';
    }
    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, char * att_name, std::string * att_value)
{
    size_t length = 0;
    int status = -1;

    status = nc_inq_attlen(ncid, i_var, att_name, &length);
    if (status != NC_NOERR)
    {
        *att_value = "";
    }
    else
    {
        char * tmp_value = (char *)malloc(sizeof(char) * (length + 1));
        status = nc_get_att(ncid, i_var, att_name, tmp_value);
        tmp_value[length] = '\0';
        *att_value = std::string(tmp_value, length);
        free(tmp_value);
        tmp_value = nullptr;
    }
    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, std::string att_name, std::string * att_value)
{
    size_t length = 0;
    int status = -1;

    status = nc_inq_attlen(ncid, i_var, att_name.c_str(), &length);
    if (status != NC_NOERR)
    {
        *att_value = "";
    }
    else
    {
        char * tmp_value = (char *)malloc(sizeof(char) * (length + 1));
        status = nc_get_att(ncid, i_var, att_name.c_str(), tmp_value);
        tmp_value[length] = '\0';
        *att_value = std::string(tmp_value, length);
        free(tmp_value);
        tmp_value = nullptr;
    }
    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, char * att_name, double * att_value)
{
    int status = -1;

    status = nc_get_att_double(ncid, i_var, att_name, att_value);

    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, char * att_name, int * att_value)
{
    int status = -1;

    status = nc_get_att_int(ncid, i_var, att_name, att_value);

    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_attribute(int ncid, int i_var, char * att_name, long * att_value)
{
    int status = -1;

    status = nc_get_att_long(ncid, i_var, att_name, att_value);

    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_dimension(int ncid, char * dim_name, size_t * dim_length)
{
    int dimid;
    int status = -1;

    if (dim_name != NULL && strlen(dim_name) != 0)
    {
        status = nc_inq_dimid(ncid, dim_name, &dimid);
        if (status == NC_NOERR)
        {
            status = nc_inq_dimlen(ncid, dimid, dim_length);
        }
    }
    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_dimension(int ncid, std::string dim_name, size_t * dim_length)
{
    int dimid;
    int status = -1;

    if (dim_name.size() != 0)
    {
        status = nc_inq_dimid(ncid, dim_name.c_str(), &dimid);
        if (status == NC_NOERR)
        {
            status = nc_inq_dimlen(ncid, dimid, dim_length);
        }
    }
    return status;
}
//------------------------------------------------------------------------------
int SGRID::get_dimension_var(int ncid, std::string var_name, size_t * dim_length)
{
    // get the total dimension length in bytes of the var_name variable
    int dimid;
    int status = -1;

    if (var_name.size() != 0)
    {
        int janm;
        status = nc_inq_varid(ncid, var_name.c_str(), &dimid);
        if (status == NC_NOERR)
        {
            *dim_length = 0;
            status = nc_inq_vardimid(ncid, dimid, &janm);
            char* tmp_value = (char*)malloc(sizeof(char) * (NC_MAX_NAME + 1));;
            status = nc_inq_dimname(ncid, janm, tmp_value);
            status = get_dimension(ncid, tmp_value, dim_length);
            free(tmp_value);
            tmp_value = nullptr;
        }
    }
    return status;
}
std::vector<std::string> SGRID::get_dimension_names(int ncid, std::string var_name)
{
    int varid;
    int status = -1;
    std::vector<std::string> dim_names;

    if (var_name.size() != 0)
    {
        int nr_dims;
        int * var_dimids = NULL;
        status = nc_inq_varid(ncid, var_name.c_str(), &varid);
        status = nc_inq_varndims(ncid, varid, &nr_dims);
        if (nr_dims > 0) {
            var_dimids = (int *)malloc(sizeof(int) * nr_dims);
        }
        status = nc_inq_var(ncid, varid, NULL, NULL, &nr_dims, var_dimids, NULL);
        for (int i = 0; i < nr_dims; i++)
        {
            //status = nc_inq_vardimid(this->m_ncid, dimid, janm);
            char * tmp_value = (char *)malloc(sizeof(char) * (NC_MAX_NAME + 1));;
            status = nc_inq_dimname(ncid, var_dimids[i], tmp_value);
            dim_names.push_back(std::string(tmp_value));
            free(tmp_value);
            tmp_value = nullptr;
        }
    }
    return dim_names;
}

//------------------------------------------------------------------------------
std::vector<std::string> SGRID::tokenize(const std::string& s, char c) {
    auto end = s.cend();
    auto start = end;

    std::vector<std::string> v;
    for (auto it = s.cbegin(); it != end; ++it) {
        if (*it != c) {
            if (start == end)
                start = it;
            continue;
        }
        if (start != end) {
            v.emplace_back(start, it);
            start = end;
        }
    }
    if (start != end)
        v.emplace_back(start, end);
    return v;
}
//
std::vector<std::string> SGRID::tokenize(const std::string& s, std::size_t count)
{
    size_t minsize = s.size() / count;
    std::vector<std::string> tokens;
    for (size_t i = 0, offset = 0; i < count; ++i)
    {
        size_t size = minsize;
        if ((offset + size) < s.size())
            tokens.push_back(s.substr(offset, size));
        else
            tokens.push_back(s.substr(offset, s.size() - offset));
        offset += size;
    }
    return tokens;
}

char * SGRID::strndup(const char *s, size_t n)
{
    char *result;
    size_t len = strlen(s);

    if (n < len)
        len = n;

    result = (char *)malloc(len + 1);
    if (!result)
        return 0;

    result[len] = '\0';
    return (char *)memcpy(result, s, len);
}
double* SGRID::permute_array(double* in_arr, std::vector<long> order, std::vector<long> dim)
{
    //
    // in_array: array which need to be permuted
    // order: order of the dimensions in the supplied array in_array
    // dim: dimensions in required order
    // return value: permuted array in required order
    //
    if (dim.size() != 4) {
        return nullptr;
    }
    long k_source;
    long k_target;
    long length = 1;
    std::vector<long> cnt(4);
    for (long i = 0; i < dim.size(); ++i)
    {
        length *= dim[i];
    }
    double* out_arr = (double*)malloc(sizeof(double) * length);
    for (long i = 0; i < dim[0]; i++) {
        cnt[0] = i;
        for (long j = 0; j < dim[1]; j++) {
            cnt[1] = j;
            for (long k = 0; k < dim[2]; k++) {
                cnt[2] = k;
                for (long l = 0; l < dim[3]; l++) {
                    cnt[3] = l;
                    k_source = get_index_in_c_array(cnt[order[0]], cnt[order[1]], cnt[order[2]], cnt[order[3]], dim[order[0]], dim[order[1]], dim[order[2]], dim[order[3]]);
                    k_target = get_index_in_c_array(i, j, k, l, dim[0], dim[1], dim[2], dim[3]);
                    out_arr[k_target] = in_arr[k_source];
                }
            }
        }
    }
    free(in_arr);
    return out_arr;
}
long SGRID::get_index_in_c_array(long t, long l, long s, long xy, long time_dim, long layer_dim, long sed_dim, long xy_dim)
{
    return t * layer_dim * sed_dim * xy_dim + l * sed_dim * xy_dim + s * xy_dim + xy;
}
