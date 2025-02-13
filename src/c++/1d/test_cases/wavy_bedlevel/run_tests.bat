set exe_dir=..\..\1d_wave\bin\x64\Release

set test_name=001_stationary
del /q .\%test_name%\output\*.*
%exe_dir%\1d_wave.exe --toml .\%test_name%\input_cpp_dx=00d00.toml
del /q .\output_dx=00d00\*.*
xcopy /s /y output output_dx=00d00


set test_name=002_temporal
del /q .\%test_name%\output\*.*
%exe_dir%\1d_wave.exe --toml .\%test_name%\input_cpp_dx=01d25.toml
del /q .\output_dx=01d25\*.*
xcopy /s /y output output_dx=01d25

set test_name=002_temporal
del /q .\%test_name%\output\*.*
%exe_dir%\1d_wave.exe --toml .\%test_name%\input_cpp_dx=02d50.toml
del /q .\output_dx=02d50\*.*
xcopy /s /y output output_dx=02d50

set test_name=002_temporal
del /q .\%test_name%\output\*.*
%exe_dir%\1d_wave.exe --toml .\%test_name%\input_cpp_dx=05d00.toml
del /q .\output_dx=05d00\*.*
xcopy /s /y output output_dx=05d00

set test_name=002_temporal
del /q .\%test_name%\output\*.*
%exe_dir%\1d_wave.exe --toml .\%test_name%\input_cpp_dx=10d00.toml
del /q .\output_dx=10d00\*.*
xcopy /s /y output output_dx=10d00
