set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_mm.toml
del /q .\output_mm\*.*
xcopy /s /y output output_mm

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_mm_qz.toml
del /q .\output_mm_qz\*.*
xcopy /s /y output output_mm_qz

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_mm_zq.toml
del /q .\output_mm_zq\*.*
xcopy /s /y output output_mm_zq
