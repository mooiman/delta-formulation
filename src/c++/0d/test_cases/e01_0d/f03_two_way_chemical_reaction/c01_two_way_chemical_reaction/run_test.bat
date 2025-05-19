set exe_dir=..\..\..\..\0d_two_way_chemical_reaction\bin\x64\Release
del /q .\output\*.*
%exe_dir%\two_way_chemical_reaction.exe --toml .\input_cpp.toml
