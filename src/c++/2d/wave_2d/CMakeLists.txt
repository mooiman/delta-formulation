cmake_minimum_required(VERSION 3.30.0)

project(Main VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
    message("*** MSVC ")
    set(CMAKE_CXX_COMPILER cl)
    set(CMAKE_BUILD_TYPE Release)
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin/x64/${CMAKE_BUILD_TYPE})
else()
    message("*** UCRT64 ")
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin/ucrt64)
    set(CMAKE_CXX_COMPILER g++)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(netCDF 4.9.2 REQUIRED)
if(netCDF_FOUND)
    message("*** netCDF found")
else()
    message("*** netCDF not found")
endif()

add_subdirectory(packages)

