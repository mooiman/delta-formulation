message("********************************") 
message("*** Enter packages directory ***") 
message("********************************") 

set(exe_name burgers_eq)

# include_directories(include)
set(headers_files
    include/adv_diff_init_concentration.h
    include/adv_diff_init_velocity.h
    include/adv_diff_linear_operator.h
    include/boundary_condition.h
    include/cfts.h
    include/compile_date_and_time.h
    include/definition_map_file.h
    include/data_input_struct.h
    include/definition_map_file.h
    include/main_version.h
    include/observation_stations.h
    include/perf_timer.h
    include/print_matrix.h
    include/read_input_toml_file.h
    include/regularization.h
    include/ugrid1d.h
)    

set(no_group_source_files
    "include/version_number.ini"
)
set(source_files
    src/main.cpp
    src/adv_diff_boundary_condition.cpp
    src/adv_diff_init_concentration.cpp
    src/adv_diff_init_velocity.cpp
    src/cfts.cpp
    src/definition_map_file.cpp
    src/main_version.cpp
    src/observation_stations.cpp
    src/perf_timer.cpp
    src/print_matrix.cpp
    src/read_input_toml_file.cpp
    src/regularization.cpp
    src/ugrid1d.cpp
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
        "$ENV{KDTREE_DIR}/."
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
    target_link_libraries(${exe_name} PRIVATE "$ENV{KDTREE_DIR}/lib/kdtree.lib")
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
