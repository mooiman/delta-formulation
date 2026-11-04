set exe_dir=..\..\analytical_functions_1d\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_analytic.exe --toml .\input_boundary_layer.toml

del /q .\analytic_boundary_layer\*.*
mkdir .\analytic_boundary_layer
xcopy /s /y .\output  .\analytic_boundary_layer

del /q .\output
