echo off

set exec=".\bin\x64\burgers_eq.exe"
set netCDF_DIR=c:\Program Files\netCDF 4.9.2\
set PATH=c:\Qt\Qt6.10.3\6.10.3\msvc2022_64\;%PATH%
set PATH=c:\boost\Boost-1.90.0\;%PATH%

del %exec%
del _build

copy packages\include\main_version.h.vcs packages\include\main_version.h 

echo on
cmake -B _build > cmake_configure.log 2>&1
cmake --build _build --verbose --config Release > cmake_build.log 2>&1
cmake --install _build > cmake_install.log 2>&1

%exec%
