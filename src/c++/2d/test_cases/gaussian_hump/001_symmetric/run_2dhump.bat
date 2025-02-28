set exe_dir=..\..\..\2d_wave\bin\x64\Release\

del /q ".\output"
%exe_dir%\2d_wave.exe --toml .\input_cpp_2dhump.toml
del /q ".\output_2dhump"
xcopy /s /y output output_2dhump