//---- LGPL --------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2015.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation version 2.1.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------
//#include <stdio.h>
#include "string.h"
#include "heat_2d_version.h"

#if defined(WIN64)
#   define strdup _strdup
#endif

#if defined(LINUX64)
static char heat_2d_version[] = {heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build" (Linux64)"};
static char heat_2d_version_id[] = {"@(#)Mooiman, "heat_2d_program" Version "heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build" (Linux64), "__DATE__", "__TIME__""};
#elif defined(UCRT64)
static char heat_2d_version[] = { heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build " (Win64)" };
static char heat_2d_version_id[] = {"@(#)Mooiman, " heat_2d_program " Version " heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build " (Win64), " __DATE__ ", " __TIME__ "" };
#elif defined(WIN32) || defined(WIN64)
static char heat_2d_version[] = { heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build " (Win64)" };
static char heat_2d_version_id[] = {"@(#)Mooiman, " heat_2d_program " Version " heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build " (Win64), " __DATE__ ", " __TIME__ ""};
#else
static char heat_2d_version[] = {heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build" (Unknown)"};
static char heat_2d_version_id[] = {"@(#)Mooiman, "heat_2d_program" Version "heat_2d_major "." heat_2d_minor "." heat_2d_revision "." heat_2d_build" (Unknown), "__DATE__", "__TIME__""};
#endif
static char heat_2d_company_name[] = {"Mooiman"};
static char heat_2d_program_name[] = { heat_2d_program };

char * getfullversionstring_heat_2d(void)
{
    return strdup(heat_2d_version_id);
};
char * getversionstring_heat_2d(void)
{
    return strdup(heat_2d_version);
};
char * getcompanystring_heat_2d(void)
{
    return strdup(heat_2d_company_name);
};
char * getprogramstring_heat_2d(void)
{
    return strdup(heat_2d_program_name);
};
char * getsourceurlstring_heat_2d(void)
{
    return strdup(heat_2d_source_url);
};
