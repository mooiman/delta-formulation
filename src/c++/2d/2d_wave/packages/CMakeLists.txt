set(exe_name 2d_wave)

# include_directories(include)
set(headers_files
    include/bed_level.h
    include/bed_shear_stress.h
    include/boundary_condition.h
    include/build_matrix_pattern.h
    include/build_matrix_pattern_regularization.h
    include/cfts.h
    include/compile_date_and_time.h
    include/convection.h
    include/data_input_struct.h
    include/data_struct.h
    include/diffusion2d.h
    include/grid.h
    include/initial_conditions.h
    include/interpolations.h
    include/jacobians.h
    include/janm
    include/main_version.h
    include/main_version.h.vcs
    include/matrix_assembly_boundaries.h
    include/matrix_assembly_corners.h
    include/matrix_assembly_interior.h
    include/matrix_assembly_utilde_boundaries.h
    include/matrix_assembly_utilde_corners.h
    include/matrix_assembly_utilde_interior.h
    include/observation_stations.h
    include/perf_timer.h
    include/read_input_toml_file.h
    include/regularization.h
    include/ugrid2d.h
    include/version_number.ini
    include/viscosity.h
)    

set(source_files
    src/bed_level.cpp
    src/bed_shear_stress.cpp
    src/boundary_condition.cpp
    src/build_matrix_pattern.cpp
    src/build_matrix_pattern_regularization.cpp
    src/cfts.cpp
    src/convection.cpp
    src/diffusion2d.cpp
    src/grid.cpp
    src/initial_conditions.cpp
    src/interpolations.cpp
    src/jacobians.cpp
    src/main_2d.cpp
    src/main_version.cpp
    src/matrix_assembly_boundaries.cpp
    src/matrix_assembly_corners.cpp
    src/matrix_assembly_interior.cpp
    src/matrix_assembly_utilde_boundaries.cpp
    src/matrix_assembly_utilde_corners.cpp
    src/matrix_assembly_utilde_interior.cpp
    src/observation_stations.cpp
    src/perf_timer.cpp
    src/read_input_toml_file.cpp
    src/regularization.cpp
    src/ugrid2d.cpp
    src/viscosity.cpp
    src/write_input_toml_to_logfile.cpp
)

add_executable(${exe_name} ${source_files} ${headers_files})

if (MSVC)
    message("*** MSVC ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/eigen"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/toml++"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/amgcl"
        "$ENV{BOOST_ROOT_DIR}/."
        "$ENV{NETCDF_DIR}/."
        "$ENV{KDTREE_DIR}/."
        "c:/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools/VC/Tools/MSVC/14.40.33807/include/"
        "c:/Program Files (x86)/Windows Kits/10/Include/10.0.22621.0/ucrt/"
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
    target_link_libraries(${exe_name} PRIVATE "$ENV{NETCDF_DIR}/lib/netcdf.lib")
    target_link_libraries(${exe_name} PRIVATE "$ENV{KDTREE_DIR}/lib/kdtree.lib")
else()
    target_link_libraries(${exe_name} PRIVATE netcdf)
endif()

set_property(TARGET ${exe_name} PROPERTY "/SUBSYSTEM:WINDOWS")

set(install_dir ${CMAKE_BINARY_DIR})
set(build_dir ${CMAKE_BINARY_DIR})
install(TARGETS ${exe_name} RUNTIME DESTINATION bin)

set(netcdf_dlls
    $ENV{NETCDF_DIR}/bin/netcdf.dll
    $ENV{NETCDF_DIR}/bin/hdf5_hl.dll
    $ENV{NETCDF_DIR}/bin/libcurl.dll
    $ENV{NETCDF_DIR}/bin/hdf5.dll
    $ENV{NETCDF_DIR}/bin/zlib1.dll
)
add_custom_command(
        TARGET ${exe_name}  
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${netcdf_dlls} 
                ${EXECUTABLE_OUTPUT_PATH}
)
