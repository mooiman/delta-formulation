set exe_dir=..\..\..\..\1d_wave\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_1d.toml
