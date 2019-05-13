@echo off
cls

REM loading visual studio 2015 (v14) environment (to get nmake)
call "%vs140comntools%vsvars32.bat"

set G4_bat_file_dir=%~dp0

cd %G4_bat_file_dir%\build

if exist CMakeCache.txt del CMakeCache.txt

%G4_bat_file_dir%\cmake_win32\bin\cmake.exe ^
-DCMAKE_CONFIGURATION_TYPES=%build_type% ^
-DCMAKE_PREFIX_PATH=%G4_bat_file_dir%\install\lib\Geant4-10.4.3 ^
%G4_bat_file_dir%

%G4_bat_file_dir%\cmake_win32\bin\cmake.exe --build . --config %build_type%

start "" %G4_bat_file_dir%\ExampleB1\build\%build_type%\TGF_propa.exe

PAUSE
