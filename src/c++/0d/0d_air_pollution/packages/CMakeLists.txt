message("********************************") 
message("*** Enter packages directory ***") 
message("********************************") 

set(exe_name main)

# include_directories(include)
set(headers_files
    include/cfts.h
    include/compile_date_and_time.h
    include/data_input_struct.h
    include/main_version.h
    include/perf_timer.h
    include/read_input_toml_file.h
)    

set(no_group_source_files
    "include/version_number.ini"
)

set(source_files
    src/main.cpp
    src/cfts.cpp
    src/main_version.cpp
    src/perf_timer.cpp
    src/read_input_toml_file.cpp
    src/write_input_toml_to_logfile.cpp
)

add_executable(${exe_name} ${source_files} ${headers_files})

# Combine the two custom target definitions into one
add_custom_target(
    myHeaderTarget
    COMMAND powershell -ExecutionPolicy Bypass -File "${CMAKE_CURRENT_SOURCE_DIR}\\..\\scripts\\git_insert_hash.ps1" ..\\..\\packages\\include main_version.h
)

# Ensure myHeaderTarget runs before the executable is built
add_dependencies(${exe_name} myHeaderTarget)

if (MSVC)
    message("*** MSVC ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/eigen"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/toml++"
        "$ENV{BOOST_ROOT_DIR}/."
        "$ENV{NETCDF_DIR}/."
    )
    message("--- MSVC ")
else()
    message("*** UCRT64 ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "c:/msys64/ucrt64/lib"
        "c:/msys64/ucrt64/include/boost/"
    )
    message("--- UCRT64 ")
endif()

if (MSVC)
    target_link_libraries(${exe_name} PRIVATE "$ENV{NETCDF_DIR}/lib/netcdf.lib")
else()
    target_link_libraries(${exe_name} PRIVATE netcdf)
endif()

set_property(TARGET ${exe_name} PROPERTY true)

set(install_dir ${CMAKE_BINARY_DIR})
set(build_dir ${CMAKE_BINARY_DIR})
install(TARGETS ${exe_name} DESTINATION .)

set(netcdf_dlls
    $ENV{NETCDF_DIR}/bin/netcdf.dll
    $ENV{NETCDF_DIR}/bin/hdf5_hl.dll
    $ENV{NETCDF_DIR}/bin/libcurl.dll
    $ENV{NETCDF_DIR}/bin/hdf5.dll
    $ENV{NETCDF_DIR}/bin/zlib1.dll
)
