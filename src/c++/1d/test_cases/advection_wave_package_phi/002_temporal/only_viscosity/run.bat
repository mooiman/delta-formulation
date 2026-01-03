set exe_dir=..\..\..\1d_advection_phi\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_adv_phi.exe --toml .\input_cpp.toml
