#ifndef __DEFINITION_MAP_FILE_H__
#define __DEFINITION_MAP_FILE_H__

#include <cstdlib>
#include <vector>
#include "ugrid1d.h"

UGRID1D * create_map_file(std::string, std::string, std::vector<double> &, std::vector<std::string>&);

#endif __DEFINITION_MAP_FILE_H__
