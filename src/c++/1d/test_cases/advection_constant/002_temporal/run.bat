set exe_dir=..\..\..\1d_advection\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_adv.exe --toml .\input_cpp.toml
