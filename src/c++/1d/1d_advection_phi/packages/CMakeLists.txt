set(exe_name 1d_adv_phi)

# include_directories(include)
set(headers_files
    include/adv_diff_boundary_condition.h
    include/adv_diff_init_concentration.h
    include/adv_diff_init_velocity.h
    include/adv_diff_linear_operator.h
    include/adv_diff_regularization.h
    include/adv_diff_source.h
    include/cfts.h
    include/definition_map_file.h
    include/perf_timer.h
    include/regularization.h
    include/ugrid1d.h
)    

set(source_files
    src/main.cpp
    src/adv_diff_boundary_condition.cpp
    src/adv_diff_init_concentration.cpp
    src/adv_diff_init_velocity.cpp
    src/adv_diff_regularization.cpp
    src/adv_diff_source.cpp
    src/cfts.cpp
    src/definition_map_file.cpp
    src/perf_timer.cpp
    src/regularization.cpp
    src/ugrid1d.cpp
)

add_executable(${exe_name} ${source_files} ${headers_files})

if (MSVC)
    message("*** MSVC ")
    target_include_directories(${exe_name} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/."
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/eigen"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/toml++"
        "$ENV{NETCDF_DIR}/."
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
