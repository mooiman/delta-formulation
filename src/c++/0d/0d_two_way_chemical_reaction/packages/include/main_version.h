#ifndef __main_VERSION_H__
#define __main_VERSION_H__

#include <string>

#define main_major "1"
#define main_minor "00"
#define main_revision "000"
#define main_git_build "M-c096595"
#define main_git_branch "main"
#define main_git_date "2026-10-29 16:15:11"

#define main_company "Mooiman"
#define main_company_url "https://github.com/mooiman"
#define main_git_source_url "https://github.com/mooiman/delta-formulation.git"
#define main_program "Two_way_chemical_reaction"

/*=================================================== DO NOT MAKE CHANGES BELOW THIS LINE ===================================================================== */

extern "C" {
    extern char * getcompanystring_main(void);
    extern char * getfullversionstring_main(void);
    extern char * getprogramstring_main(void);
    extern char * getversionstring_main(void);
    extern char * getgitbranchstring_main(void);
    extern char * getgiturlstring_main(void);
    extern char * getgitbuildstring_main(void);
    extern char * getgitdatestring_main(void);
}

#endif // __main_VERSION_H__
