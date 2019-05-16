set bat_file_dir=%~dp0
cd %bat_file_dir%

REM aguments passed to the executable are (in order)
REM - number of particles to record before stoping
REM - TGF altitude (km)
REM - TGF latitude (deg)
REM - TGF longitude (deg)
REM - TGF time distribution (sigma of gaussian)
REM - TGF opening angle (deg) : sigma of gaussian if "Gaussian", half cone angle if "Uniform"
REM - TGF beaming tilt angle (deg)
REM - TGF beaming type (Uniform or Gaussian)
REM - record altitude (km)

start "" "%bat_file_dir%TGF_Propa.exe" 1000 14 11.01 -95.4 0 35 0 Uniform 400