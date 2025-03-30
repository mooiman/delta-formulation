set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx00d625.toml
del /q ".\output_dx00d625"
mkdir ".\output_dx00d625"
xcopy /s /y output output_dx00d625

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx01d25.toml
del /q ".\output_dx01d25"
mkdir ".\output_dx01d25"
xcopy /s /y output output_dx01d25

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx02d50.toml
del /q ".\output_dx02d50"
mkdir ".\output_dx02d50"
xcopy /s /y output output_dx02d50

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx05d00.toml
del /q ".\output_dx05d00"
mkdir ".\output_dx05d00"
xcopy /s /y output output_dx05d00

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx10d00.toml
del /q ".\output_dx10d00"
mkdir ".\output_dx10d00"
xcopy /s /y output output_dx10d00
