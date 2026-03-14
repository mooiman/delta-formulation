set exe_dir=..\..\..\..\wave_1d\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\wave_1d.exe --toml .\input_cpp.toml

