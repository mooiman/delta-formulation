echo off

set exec=".\bin\x64\main.exe"
set netCDF_DIR=c:\Program Files\netCDF 4.9.2\
set PATH=c:\Qt\Qt6.7.2\6.7.2\msvc2019_64\;%PATH%
set PATH=c:\boost\Boost-1.85.0\;%PATH%

del %exec%
del _build

copy packages\include\main_version.h.vcs packages\include\main_version.h 
copy packages\include\main_version.rc.vcs packages\include\main_version.rc 

echo on
cmake -B _build > cmake_configure.log 2>&1
cmake --build _build --verbose > cmake_build.log 2>&1

%exec% --toml ./input_cpp.toml
