set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp.toml
rem del /q .\output_mm\*.*
rem xcopy /s /y output output_mm
