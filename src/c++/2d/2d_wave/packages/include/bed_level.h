#ifndef __BED_LEVEL_H__
#define __BED_LEVEL_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


class BED_LEVEL
{
public:
    BED_LEVEL();
    ~BED_LEVEL();
    long open(std::string filename);
    long read(int nx, int ny);
    std::vector<double> get_bed_level();

private:
    std::ifstream m_fname;
    std::vector<double> m_bed_given;
};

#endif  // __BED_LEVEL_H__
