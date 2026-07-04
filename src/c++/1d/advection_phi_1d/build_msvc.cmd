echo off

set netCDF_DIR=c:\Program Files\netCDF 4.9.3\
set PATH=c:\Qt\Qt6.7.2\6.7.2\msvc2019_64\;%PATH%
set PATH=c:\boost\Boost-1.87.0\;%PATH%

del %exec%
del _build

echo on
cmake -B _build > cmake_configure.log 2>&1
cmake --build _build --verbose > cmake_build.log 2>&1

