#include "bed_level.h"

//------------------------------------------------------------------------------
BED_LEVEL::BED_LEVEL()
{
}
//------------------------------------------------------------------------------
BED_LEVEL::~BED_LEVEL()
{
}
//------------------------------------------------------------------------------
long BED_LEVEL::open(std::string filename)
{
    long status = 1;
    m_fname.open(filename);
    if (m_fname.is_open()) {
        status = 0;
    }
    return status;
}
//------------------------------------------------------------------------------
long BED_LEVEL::read(int nx, int ny)
{
    long status = -1;
    std::istringstream iss;
    double missing_value = -999.0;  // used by quickin

    m_bed_given.resize(nx*ny);

    int k = -1;
    std::string token;
    std::string line;
    std::getline(m_fname, line);
    iss.clear();
    iss.str(line);
    for (;;)
    {
        while (iss >> token)
        {
            if (missing_value != stod(token))
            {
                ++k;
                m_bed_given[k] = stod(token);
            }
        }
        std::getline(m_fname, line);
        iss.clear();
        iss.str(line);
        if (std::fmod(k + 1, nx*ny) == 0) { break; }
    }
    if (k + 1 == nx * ny)
    {
        status = 0;
    }

    return status;
}
//------------------------------------------------------------------------------
std::vector<double> BED_LEVEL::get_bed_level()
{
     return m_bed_given;
}
