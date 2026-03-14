echo off

set netCDF_DIR=c:\Program Files\netCDF 4.9.2\
set PATH=c:\Qt\Qt6.10.2\6.10.2\msvc2022_64\;%PATH%
set PATH=c:\boost\Boost-1.90.0\;%PATH%

echo on

del _build
cmake -B _build > cmake_configure.log 2>&1
cmake --build _build --verbose > cmake_build.log 2>&1

