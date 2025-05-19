set exe_dir=..\..\..\1d_advection_diffusion\bin\x64\Release

del /q .\output\*.*
del /q ".\output_pe05d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe05d0.toml
mkdir ".\output_pe05d0"
xcopy /s /y output output_pe05d0

del /q .\output\*.*
del /q ".\output_pe02d6"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe02d6.toml
mkdir ".\output_pe02d6"
xcopy /s /y output output_pe02d6

del /q .\output\*.*
del /q ".\output_pe02d5"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe02d5.toml
mkdir ".\output_pe02d5"
xcopy /s /y output output_pe02d5

del /q .\output\*.*
del /q ".\output_pe02d4"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe02d4.toml
mkdir ".\output_pe02d4"
xcopy /s /y output output_pe02d4

del /q .\output\*.*
del /q ".\output_pe02d1"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe02d1.toml
mkdir ".\output_pe02d1"
xcopy /s /y output output_pe02d1

del /q .\output\*.*
del /q ".\output_pe01d9"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe01d9.toml
mkdir ".\output_pe01d9"
xcopy /s /y output output_pe01d9

del /q .\output\*.*
del /q ".\output_pe01d0"
%exe_dir%\1d_adv_diff.exe --toml .\input_cpp_pe01d0.toml
mkdir ".\output_pe01d0"
xcopy /s /y output output_pe01d0
