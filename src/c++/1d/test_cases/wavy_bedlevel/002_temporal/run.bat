set exe_dir=..\..\..\1d_wave\bin\x64\Release

rem del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp.toml
