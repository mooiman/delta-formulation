set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
del /q ".\output_dx20d0_dt4d0"
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx20d0_dt4d0.toml
mkdir ".\output_dx20d0_dt4d0"
xcopy /s /y output output_dx20d0_dt4d0

del /q .\output\*.*
del /q ".\output_dx10d0_dt2d0"
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx10d0_dt2d0.toml
mkdir ".\output_dx10d0_dt2d0"
xcopy /s /y output output_dx10d0_dt2d0

del /q .\output\*.*
del /q ".\output_dx05d0_dt1d0"
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx05d0_dt1d0.toml
mkdir ".\output_dx05d0_dt1d0"
xcopy /s /y output output_dx05d0_dt1d0

del /q .\output\*.*
del /q ".\output_dx02d5_dt0d5"
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx02d5_dt0d5.toml
mkdir ".\output_dx02d5_dt0d5"
xcopy /s /y output output_dx02d5_dt0d5

del /q .\output\*.*
del /q ".\output_dx01d25_dt0d25"
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx01d25_dt0d25.toml
mkdir ".\output_dx01d25_dt0d25"
xcopy /s /y output output_dx01d25_dt0d25
