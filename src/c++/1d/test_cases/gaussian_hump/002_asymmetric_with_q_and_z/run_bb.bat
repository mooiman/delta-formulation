set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb.toml
del /q .\output_bb\*.*
xcopy /s /y output output_bb

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_qz.toml
del /q .\output_bb_qz\*.*
xcopy /s /y output output_bb_qz

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bb_zq.toml
del /q .\output_bb_zq\*.*
xcopy /s /y output output_bb_zq
