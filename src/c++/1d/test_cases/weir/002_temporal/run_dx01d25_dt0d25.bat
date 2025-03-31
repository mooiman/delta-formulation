set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_dx01d25_dt0d25.toml
del /q ".\output_dx01d25_dt0d25"
mkdir ".\output_dx01d25_dt0d25"
xcopy /s /y output output_dx01d25_dt0d25
