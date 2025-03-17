set exe_dir=..\..\..\1d_wave\bin\x64\Release

del /q .\output\*.*
%exe_dir%\1d_wave.exe --toml .\input_cpp_mb.toml
del /q .\output_mb\*.*
mkdir .\output_mb
xcopy /s /y output output_mb
