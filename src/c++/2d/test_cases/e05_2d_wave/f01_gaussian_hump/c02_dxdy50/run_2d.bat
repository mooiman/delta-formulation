set exe_dir=..\..\..\2d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\2d_wave.exe --toml .\input_2d.toml
