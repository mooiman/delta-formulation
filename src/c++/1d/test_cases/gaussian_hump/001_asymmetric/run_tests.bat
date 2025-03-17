set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_zz.toml
del /q .\output_bb_zz\*.*
mkdir .\output_bb_zz
xcopy /s /y output output_bb_zz

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_qz.toml
del /q .\output_bb_qz\*.*
mkdir .\output_bb_qz
xcopy /s /y output output_bb_qz

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_zq.toml
del /q .\output_bb_zq\*.*
mkdir .\output_bb_zq
xcopy /s /y output output_bb_zq

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_qq.toml
del /q .\output_bb_qq\*.*
mkdir .\output_bb_qq
xcopy /s /y output output_bb_qq
