set(exe_name heat_2d)

# include_directories(include)
set(headers_files
    include/bed_level.h
    include/bed_shear_stress.h
    include/boundary_condition.h
    include/cfts.h
    include/convection.h
    include/data_struct.h
    include/grid.h
    include/initial_conditions.h
    include/perf_timer.h
    include/regularization.h
    include/ugrid2d.h
)    

set(source_files
    src/bed_level.cpp
    src/bed_shear_stress.cpp
    src/boundary_condition.cpp
    src/cfts.cpp
    src/convection.cpp
    src/grid.cpp
    src/initial_conditions.cpp
    src/main_2d.cpp
    src/perf_timer.cpp
    src/regularization.cpp
    src/ugrid2d.cpp
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
