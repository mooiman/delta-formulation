echo off
set NETCDF_DIR=c:\Program Files\netCDF 4.9.2\
set PATH=c:\Qt\Qt6.8.2\6.8.2\msvc2022_64\;%PATH%
set PATH=c:\boost\Boost-1.90.0\;%PATH%

echo on
del _build
cmake -B _build > cmake_configure.log 2>&1
cmake --build _build --verbose --config Release > cmake_build.log 2>&1

