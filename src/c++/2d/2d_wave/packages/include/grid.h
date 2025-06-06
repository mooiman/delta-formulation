#ifndef __GRID_H__
#define __GRID_H__

#include <map>
#include "data_struct.h"
#include "include/netcdf.h"

////////////////////////////////////////////////////////////////////////////////
enum class FILE_TYPE {
    UGRID,
    SGRID,
    KISS,
    UNKNOWN
};
struct _variable {
public:
    _variable() {};
    double fill_value;
    std::vector<long> dims;
    nc_type nc_type;
    int topology_dimension;
    std::string var_name;
    std::string location;
    std::string mesh;
    std::string coordinates;
    std::string cell_methods;
    std::string standard_name;
    std::string long_name;
    std::string units;
    std::string grid_mapping;
    std::string comment;
    std::string description;
    std::vector<int> flag_values;
    std::vector<std::string> flag_meanings;
    //
    bool draw;
    bool read;
    bool time_series;
    std::vector<std::string> dim_names;
    int nr_hydro_layers = -1;
    int nr_bed_layers = -1;
    int sediment_index = -1;
    int sediment_array = -1;
    std::vector<double> layer_center;
    std::vector<double> layer_interface;
    std::map<std::string, std::string> map_dim_name;
};
struct _feature {
    size_t count;
    double fill_value;
    std::vector<int> dims;
    std::vector<long> branch;  // 1D, on which branch the node is
    std::vector<double> chainage;  // 1D, chainage of the node
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::string> name;
    std::vector<std::string> long_name;
};

////////////////////////////////////////////////////////////////////////////////
struct _mesh_variable {
    long nr_vars;
    struct _variable** variable;
};

struct _mesh2d {
    long nr_mesh2d;
    struct _feature ** node;
    struct _edge ** edge;
    struct _feature ** face;  // face: the location of the mass centre
    std::vector<std::vector<int>> face_nodes;  // face: the location of the mass centre
};

struct _edge {
    size_t count;
    double fill_value;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::vector<double>> x_bounds;  // begin- and endpoint of edge when drawing quantities, not necessarily the begin an end point
    std::vector<std::vector<double>> y_bounds;  // begin- and endpoint of edge when drawing quantities, not necessarily the begin an end point
    int** edge_nodes;
    std::vector<long> edge_branch;
    std::vector<double> edge_length;
    std::vector<std::string> name;
    std::vector<std::string> long_name;
};

struct _mesh2d_string {
    std::string var_name;
    std::string long_name;
    std::string edge_coordinates;  // 2D
    std::string x_edge_name;  // 2D
    std::string y_edge_name;  // 2D
    std::string face_coordinates;  // 2D
    std::string x_face_name;  // 2D
    std::string y_face_name;  // 2D
    std::string node_coordinates;  // 1D, 2D
    std::string x_node_name;  // 2D
    std::string y_node_name;  // 2D
    size_t topology_dimension;
    //
    std::string edge_dimension;
    std::string edge_geometry;
    std::string edge_length;
    std::string edge_face_connectivity;
    std::string edge_node_connectivity;
    std::string edge_type;  // internal_closed, internal, boundary, boundary_closed
    std::string face_dimension;
    std::string face_edge_connectivity;
    std::string face_face_connectivity;
    std::string face_node_connectivity;
    std::string max_face_nodes_dimension;
    std::string node_edge_exchange;
    std::string layer_dimension;
    std::string layer_interface_dimension;
    //
    std::vector<std::string> dim_name;  // in case of 3D-sigma simualtion
    std::string x_bound_edge_name;
    std::string y_bound_edge_name;
    std::string x_bound_face_name;
    std::string y_bound_face_name;
};

class SGRID
{
public:
    SGRID();
    ~SGRID();
    SGRID(std::string filename);
    long open(std::string filename);
    std::string get_filename();

    long read();
    long read_global_attributes();
    long get_global_attribute_value(std::string, std::string*);

    long read_sgrid_mesh();
    long read_sgrid_variables();
    struct _mesh2d_string ** get_mesh2d_string();

    struct _global_attributes * get_global_attributes(void);

    struct _mesh2d * get_mesh_2d();
    struct _mapping * get_grid_mapping();

    long set_grid_mapping_epsg(long, std::string);
    std::vector<std::string> get_names(int, std::string, size_t);

    std::vector<std::string> tokenize(const std::string & , const char );
    std::vector<std::string> tokenize(const std::string & , std::size_t );

    struct _mesh_variable * get_variables();
    struct _variable * get_var_by_std_name(struct _mesh_variable *, std::string, std::string);

    FILE_TYPE get_file_type();

private:
    char * strndup(const char *, size_t);
    int read_mesh2d_attributes(struct _mesh2d_string *, int, std::string, int);

    int get_attribute_by_var_name(int, std::string, std::string, std::string *);
    int get_attribute(int, int, char *, char **);
    int get_attribute(int, int, char *, std::string *);
    int get_attribute(int, int, std::string, std::string *);
    int get_attribute(int, int, char *, double *);
    int get_attribute(int, int, char *, int *);
    int get_attribute(int, int, char *, long *);
    int get_dimension(int, char *, size_t *);
    int get_dimension(int, std::string, size_t *);
    int get_dimension_var(int, std::string, size_t *);
    std::vector<std::string> get_string_var(int, std::string);
    std::vector<std::string> get_dimension_names(int, std::string);

    int read_variables_with_cf_role(int, std::string, std::string, int, int *);
    int read_grid_mapping(int, std::string, std::string);

    double * permute_array(double *, std::vector<long>, std::vector<long>);
    long get_index_in_c_array(long, long, long, long, long, long, long, long);

    struct _mesh2d * m_mesh2d = nullptr;
    struct _mapping * m_mapping = nullptr;
    struct _mesh_variable * m_mesh_vars = nullptr;
    struct _feature * m_nodes = nullptr;
    struct _mesh2d_string** m_mesh2d_strings = nullptr;

    long m_nr_mesh_var;
    size_t _two;
    size_t * m_dimids;
    std::vector<std::string> m_dim_names;
    std::map<std::string, long> m_map_dim;
    std::map<std::string, std::string> m_map_dim_name;

    FILE_TYPE m_ftype;
    int m_ncid;
    int * topo_edge_nodes;
    int * mesh2d_edge_nodes;
    std::string m_fname;
    std::string m_grid_file_name;
    struct _global_attributes * global_attributes;
};

#endif  // __GRID_H__
