set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx=01d25.toml
del /q ".\output_dx=01d25"
mkdir ".\output_dx=01d25"
xcopy /s /y output output_dx=01d25

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx=02d50.toml
del /q ".\output_dx=02d50"
mkdir ".\output_dx=02d50"
xcopy /s /y output output_dx=02d50

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx=05d00.toml
del /q ".\output_dx=05d00"
mkdir ".\output_dx=05d00"
xcopy /s /y output output_dx=05d00

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx=10d00.toml
del /q ".\output_dx=10d00"
mkdir ".\output_dx=10d00"
xcopy /s /y output output_dx=10d00
