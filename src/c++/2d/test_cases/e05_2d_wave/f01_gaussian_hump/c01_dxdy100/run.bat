set exe_dir=..\..\..\2d_wave\bin\x64\Release\

del /q ".\output"
del /q ".\output_2dhump"
%exe_dir%\2d_wave.exe --toml .\input_cpp.toml
copy ..\..\..\2d_wave\packages\src\main_2d.cpp output\main_2d.cpp
mkdir ".\output_2dhump"
xcopy /s /y output output_2dhump
S