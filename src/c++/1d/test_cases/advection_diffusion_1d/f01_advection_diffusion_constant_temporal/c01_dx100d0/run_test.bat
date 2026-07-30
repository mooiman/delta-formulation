set exe_dir=..\..\..\..\advection_diffusion_1d\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\adv_diff_1d.exe --toml .\input_cpp.toml
