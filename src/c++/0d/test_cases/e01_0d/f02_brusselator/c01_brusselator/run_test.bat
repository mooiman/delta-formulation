set exe_dir=..\..\..\..\0d_brusselator\bin\x64\Release
del /q .\output\*.*
%exe_dir%\brusselator.exe --toml .\input_cpp.toml
