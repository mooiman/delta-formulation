#ifndef __DEFINITION_MAP_FILE_H__
#define __DEFINITION_MAP_FILE_H__

#include <cstdlib>
#include <vector>
#include "ugrid1d.h"

UGRID1D* create_map_file(std::string nc_mapfile, std::string model_title, std::vector<double>& x, std::vector<std::string>& map_names);

#endif __DEFINITION_MAP_FILE_H__
