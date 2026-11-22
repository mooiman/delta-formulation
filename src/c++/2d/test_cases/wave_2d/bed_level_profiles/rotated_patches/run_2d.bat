set exe_dir=..\..\..\..\wave_2d\bin\x64\Release

del /q .\output\*.*
%exe_dir%\wave_2d.exe --toml .\input_2d.toml
