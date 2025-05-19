set exe_dir=..\..\..\..\0d_air_pollution\bin\x64\Release
del /q .\output\*.*
%exe_dir%\air_pollution.exe --toml .\input_cpp.toml
