set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_bm.toml
del /q .\output_bm\*.*
xcopy /s /y output output_bm
