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
#include "main_version.h"

#if defined(WIN64)
#   define strdup _strdup
#endif

#if defined(LINUX64)
static char main_version[] = {main_major "." main_minor "." main_revision " (Linux64)"};
static char main_version_id[] = {"@(#)Mooiman, "main_program" Version "main_major "." main_minor "." main_revision "." main_git_build" (Linux64), "__DATE__", "__TIME__""};
#elif defined(UCRT64)
static char main_version[] = { main_major "." main_minor "." main_revision " (UCRT64)" };
static char main_version_id[] = {"@(#)Mooiman, " main_program " Version " main_major "." main_minor "." main_revision "." main_git_build " (Win64), " __DATE__ ", " __TIME__ "" };
#elif defined(WIN32) || defined(WIN64)
static char main_build_string[] = { main_git_build };
static char main_version[] = { main_major "." main_minor "." main_revision " (Win64)" };
static char main_version_id[] = {"@(#)Mooiman, " main_program " Version " main_major "." main_minor "." main_revision "." main_git_build " (Win64), " __DATE__ ", " __TIME__ ""};
#else
static char main_version[] = {main_major "." main_minor "." main_revision " (Unknown)"};
static char main_version_id[] = {"@(#)Mooiman, "main_program" Version "main_major "." main_minor "." main_revision "." main_git_build" (Unknown), "__DATE__", "__TIME__""};
#endif
static char main_company_name[] = {"Mooiman"};
static char main_program_name[] = { main_program };

char * getfullversionstring_main(void)
{
    return strdup(main_version_id);
};
char * getversionstring_main(void)
{
    return strdup(main_version);
};
char * getcompanystring_main(void)
{
    return strdup(main_company_name);
};
char * getprogramstring_main(void)
{
    return strdup(main_program_name);
};
char * getgitbranchstring_main(void)
{
    return strdup(main_git_branch);
};
char * getgiturlstring_main(void)
{
    return strdup(main_git_source_url);
};
char * getgitbuildstring_main(void)
{
    return strdup(main_git_build);
};
char * getgitdatestring_main(void)
{
    return strdup(main_git_date);
};
