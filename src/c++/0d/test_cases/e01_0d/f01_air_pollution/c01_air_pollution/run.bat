set exe_dir=..\..\..\..\0d_air_pollution\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\air_pollution.exe --toml .\input_0d.toml
