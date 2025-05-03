set exe_dir=..\..\..\1d_advection_diffusion\bin\x64\Release

del /q .\output\*.*
del /q ".\output_dx100d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_dx100d0.toml
mkdir ".\output_dx100d0"
xcopy /s /y output output_dx100d0

del /q .\output\*.*
del /q ".\output_dx050d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_dx050d0.toml
mkdir ".\output_dx050d0"
xcopy /s /y output output_dx050d0

del /q .\output\*.*
del /q ".\output_dx025d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_dx025d0.toml
mkdir ".\output_dx025d0"
xcopy /s /y output output_dx025d0

del /q .\output\*.*
del /q ".\output_dx010d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_dx010d0.toml
mkdir ".\output_dx010d0"
xcopy /s /y output output_dx010d0

del /q .\output\*.*
del /q ".\output_dx005d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_dx005d0.toml
mkdir ".\output_dx005d0"
xcopy /s /y output output_dx005d0
