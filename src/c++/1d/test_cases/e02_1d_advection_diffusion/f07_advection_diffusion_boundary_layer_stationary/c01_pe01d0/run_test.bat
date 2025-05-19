set exe_dir=..\..\..\..\1d_advection_diffusion\bin\x64\Release
del /q .\output\*.*
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp.toml
