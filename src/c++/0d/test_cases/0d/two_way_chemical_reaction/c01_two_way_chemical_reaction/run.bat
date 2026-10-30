set exe_dir=..\..\..\..\0d_two_way_chemical_reaction\bin\x64\Release

del /q .\output_01\*.*
mkdir .\output_01
xcopy /s /y .\output  .\output_01

del /q .\output\*.*
%exe_dir%\two_way_chemical_reaction.exe --toml .\input_0d.toml
