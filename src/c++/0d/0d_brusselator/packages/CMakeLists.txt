message("*** Enter packages directory") 

set(exe_name main)

# include_directories(include)
set(headers_files
    include/cfts.h
    include/perf_timer.h
)    

set(no_group_source_files
    "include/version_number.ini"
)

set(source_files
    src/main.cpp
    src/cfts.cpp
    src/perf_timer.cpp
)

add_executable(${exe_name} ${source_files} ${headers_files})

#add_custom_target(
#    myHeaderTarget
#    COMMAND cmd /c ${CMAKE_CURRENT_SOURCE_DIR}\\..\\..\\..\\scripts\\git_insert_hash.cmd ..\\..\\packages\\include main_version.h 
#    )
#add_dependencies(${exe_name} myHeaderTarget)
#add_custom_command(
#    TARGET ${exe_name}
#    PRE_BUILD
#    COMMAND cmd /c ${CMAKE_CURRENT_SOURCE_DIR}\\..\\..\\..\\scripts\\git_insert_hash.cmd ..\\..\\packages\\include main_version.rc 
#    )

if (MSVC)
    message("*** MSVC ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "c:/Program Files/netCDF 4.9.2/"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/eigen"
        "c:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Professional\\VC\\Tools\\MSVC\\14.29.30133\\include\\"
        "c:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.22621.0\\ucrt\\"
    )
else()
    message("*** UCRT64 ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "c:/msys64/ucrt64/include/Qt6/"
        "c:/msys64/ucrt64/include/boost/"
    )
endif()

if (MSVC)
    target_link_libraries(${exe_name} PRIVATE "c:/Program Files/netCDF 4.9.2/lib/netcdf.lib")
else()
    target_link_libraries(${exe_name} PRIVATE netcdf)
endif()

set_property(TARGET ${exe_name} PROPERTY "/SUBSYSTEM:WINDOWS")

set(install_dir ${CMAKE_BINARY_DIR})
set(build_dir ${CMAKE_BINARY_DIR})
install(TARGETS ${exe_name} RUNTIME DESTINATION bin)

message("*** Exiting packages directory") 
