set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dirichlet_dirichlet.toml
del /q .\output_dirichlet_dirichlet\*.*
xcopy /s /y output output_dirichlet_dirichlet
