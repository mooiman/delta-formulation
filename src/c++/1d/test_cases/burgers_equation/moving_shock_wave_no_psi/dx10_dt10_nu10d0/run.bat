set exe_dir=..\..\..\..\burgers_equation\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\burgers_eq.exe --toml .\input_1d.toml
