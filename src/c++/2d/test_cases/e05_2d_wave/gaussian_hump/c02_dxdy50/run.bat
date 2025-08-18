set exe_dir=..\..\..\..\..\engines_to_compare\x64\Release

del /q .\output\*.*
%exe_dir%\2d_wave.exe --toml .\input_2d.toml
