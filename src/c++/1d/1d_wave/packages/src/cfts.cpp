//
// programmer: Jan Mooiman
// Email: jan.mooiman@outlook.com
//

#include "cfts.h"
#include "include/netcdf.h"
#include <include/ugrid1d.h>

CFTS::CFTS()
{
    m_ncid = -1;
    m_times = -1;
    m_time_units = "seconds since 2007-01-01 00:00:00";
}
CFTS::~CFTS()
{
}
int CFTS::open(std::string ncfile, std::string model_title)
{
    int status = nc_create(ncfile.c_str(), NC_NETCDF4, &m_ncid);
  
    const std::chrono::zoned_time now{ std::chrono::current_zone(), std::chrono::system_clock::now() };
    //auto date_time = std::format("{:%FT%H:%M:%OSZ%Oz (%Z)}", now);
    auto date_time = std::format("{:%F %H:%M:%OS %Oz}", now);

    // Define global attributes
    status = set_global_attribute("Title", model_title);
    status = set_global_attribute("Model", "Delta-Formulation 1D, C++");
    status = set_global_attribute("Conventions", "CF-1.8");
    status = set_global_attribute("featureType", "timeSeries");
    status = set_global_attribute("file_created", date_time);
    status = set_global_attribute("reference", "https://www.github.com/mooiman");

    int var_id;
    status = nc_def_var(m_ncid, "projected_coordinate_system", NC_INT, 0, nullptr, &var_id);
    int epsg = 28992;
    status = set_attribute(std::string("projected_coordinate_system"), std::string("epsg"), epsg);
    status = set_attribute(std::string("projected_coordinate_system"), std::string("EPSG_CODE"), std::string("EPSG:0"));

    return status;
}
int CFTS::close()
{
    int status = nc_close(m_ncid);
    return status;
}
int CFTS::add_stations(std::vector<std::string> stations, std::vector<double> x, std::vector<double> y)
{
    int dim_id;
    int i_var;
    int status = nc_def_dim(m_ncid, "nr_stations", stations.size(), &dim_id);
    std::vector<int> dimids;
    dimids.push_back(dim_id);

    size_t max_strlen = 0;
    for (int i = 0; i < stations.size(); ++i)
    {
        max_strlen = std::max(max_strlen, stations[i].size());
    }
    status = nc_def_dim(m_ncid, "max_name_len", max_strlen, &dim_id);
    dimids.push_back(dim_id);

    status = nc_def_var(m_ncid, "station_name", NC_CHAR, dimids.size(), dimids.data(), &i_var);
    status = set_attribute(std::string("station_name"), std::string("cf_role"), std::string("timeseries_id"));
    status = set_attribute(std::string("station_name"), std::string("long_name"), std::string("Station name"));
    for (int i = 0; i < stations.size(); ++i)
    {
        size_t start1[] = { i, 0 };
        size_t count1[] = { 1, std::min(max_strlen, stations[i].size()) };

        status = nc_put_vara_text(m_ncid, i_var, start1, count1, stations[i].data());
    }

    status = nc_def_var(m_ncid, "station_x", NC_DOUBLE, 1, &dimids[0], &i_var);
    status = set_attribute(std::string("station_x"), std::string("standard_name"), std::string("projection_x_coordinate"));
    status = set_attribute(std::string("station_x"), std::string("long_name"), std::string("x-coordinate"));
    status = set_attribute(std::string("station_x"), std::string("units"), std::string("m"));
    for (int i = 0; i < x.size(); ++i)
    {
        size_t start1[] = { i };
        size_t count1[] = { 1 };

        status = nc_put_vara_double(m_ncid, i_var, start1, count1, &x[i]);
    }

    status = nc_def_var(m_ncid, "station_y", NC_DOUBLE, 1, &dimids[0], &i_var);
    status = set_attribute(std::string("station_y"), std::string("standard_name"), std::string("projection_y_coordinate"));
    status = set_attribute(std::string("station_y"), std::string("long_name"), std::string("y-coordinate"));
    status = set_attribute(std::string("station_y"), std::string("units"), std::string("m"));
    for (int i = 0; i < y.size(); ++i)
    {
        size_t start1[] = { i };
        size_t count1[] = { 1 };

        status = nc_put_vara_double(m_ncid, i_var, start1, count1, &y[i]);
    }

    return status;
}
int CFTS::add_time_series(void)
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
int CFTS::add_variable(std::string var_name, std::string std_name, std::string long_name, std::string unit)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    status = nc_inq_dimid(m_ncid, "time", &dim_id);
    dimids.push_back(dim_id);
    status = nc_inq_dimid(m_ncid, "nr_stations", &dim_id);
    dimids.push_back(dim_id);

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("coordinates"), "station_x station_y station_name");
    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}
int CFTS::add_variable_without_location(std::string var_name, std::string std_name, std::string long_name, std::string unit)
{
    int dim_id;
    int i_var;
    int status = -1;
    std::vector<int> dimids;
    status = nc_inq_dimid(m_ncid, "time", &dim_id);
    dimids.push_back(dim_id);

    status = nc_def_var(m_ncid, var_name.data(), NC_DOUBLE, dimids.size(), dimids.data(), &i_var);

    status = set_attribute(var_name, std::string("standard_name"), std_name);
    status = set_attribute(var_name, std::string("long_name"), long_name);
    status = set_attribute(var_name, std::string("units"), unit);
    return status;
}
int CFTS::put_variable(std::string var_name, int i_time, std::vector<double> values)
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
            size_t start1[] = { i_time, 0};
            size_t count1[] = { 1, values.size()};

            status = nc_put_vara_double(m_ncid, i_var, start1, count1, (const double *)values.data());
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}
int CFTS::put_time(const int nst, double time)
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
            size_t start1[] = { nst };
            size_t count1[] = { 1 };

            status = nc_put_vara_double(m_ncid, i_var, start1, count1, (const double*) &time);
            break;
        }
    }
    free(tmp_var_name_c);
    //status = nc_put_var1_double(m_ncid, i_var, NC_DOUBLE, &time);
    return status;
}

////////////////////////////////////////////////////////////////////////////////
int CFTS::set_global_attribute(std::string att_name, std::string att_value)
{
    int status = nc_put_att_text(m_ncid, NC_GLOBAL, att_name.data(), att_value.size(), att_value.data());
    return status;
}
int CFTS::set_attribute(std::string var_name, std::string att_name, int att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_int(m_ncid, i_var, att_name.data(), NC_INT, 1, &att_value);
    return status;
}
int CFTS::set_attribute(std::string var_name, std::string att_name, std::string att_value)
{
    int status = -1;
    int i_var;
    status = nc_inq_varid(m_ncid, var_name.data(), &i_var);
    status = nc_put_att_text(m_ncid, i_var, att_name.data(), att_value.size(), att_value.data());
    return status;
}