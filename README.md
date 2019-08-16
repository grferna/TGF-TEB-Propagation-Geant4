TGF-TEB-Propagation-Geant4
=======
Geant4 model for Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons propagation in Earth atmosphere and environment (magnetic field).
=======

contact : <david.sarria@uib.no>

## Generalities
- Propagation of photons, electrons and positron in Earth's environment (atmosphere, ionosphere, magnetosphere), in the context of Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons beams.
- This code is probably not perfect. Feel free to suggest improvements.
- Uses mostly Geant4 features. See [documentation](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/index.html "Geant4 documentation").
- Integrates the [NRL-MSISE-00 model](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/) for the atmosphere and the [IGRF-12](http://wdc.kugi.kyoto-u.ac.jp/igrf/index.html) model for the magnetic field (now, IGRF-12 is for Linux only, as it requires a Fortran compiler). Alternatively, the [WMM](https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml) can be used (in both Linux and Windows, only available option for this last).

## Compilation, installation
- The source code of `TGF-TEB-Propagation-Geant4` is located in `src/` and the build should be done in the folder `build/`.
- Requires [Geant4](https://geant4.web.cern.ch/) compiled and [installed](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html) properly. Recommended to use geant4-10-04-patch-02 (25-May-2018). Minor changes in the source code may be required for other versions. Geant4 source code is included in folder `geant4/geant4_source/`.
### Linux
- Geant4 needs to be installed and set-up properly. Easy installation scripts of Geant4 for Linux are provided [here](https://github.com/DavidSarria89/GEANT4-easy-install-scripts)
- Once Geant4 is installed and set-up properly, open a terminal in the `build/` folder and run `cmake ../src/` (to generate the `makefile` using CMake) and then `make` to compile. It will make the executable file `build/tgf_propa`. 

**Alternatively**, for Ubuntu users, an installation script is provided:
- Go to main/root directory (that must be called `TGF-TEB-Propagation-Geant4`) and run in terminal `bash compile_install.bash`.
- it will compile and install Geant4 [from source](https://geant4.web.cern.ch/node/1604), set up the environement variables and compile the `TGF-TEB-Propagation-Geant4` code in the `build/` folder.
- Will require super user priviledges (`sudo`) to download missing dependencies. 
- alternatively, run before-hand the command `sudo apt-get install build-essential qt4-default qtcreator cmake-qt-gui gcc g++ gfortran zlib1g-dev libxerces-c-dev libx11-dev libexpat1-dev libxmu-dev libmotif-dev libboost-filesystem-dev libeigen3-dev qt4-qmake` to install dependencies, and `sudo` should not be required.
- It was successfully tested on Ubuntu 16.04 and 18.04, but is probably not free of bugs.
- If it works properly, it will launch the program that should open a Geant4 GUI window displaying a 3D sphere.

For non-Ubuntu users look at the code in the file `compile_install.bash` to check the compile and set-up steps, and read the [Geant4 installation instructions](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html), or look to the [easy install scripts](https://github.com/DavidSarria89/GEANT4-easy-install-scripts)

### Windows
- You can try the an experimental installer (`TGF_TEB_propa_1.2.1.msi`), provided in the release part of this repository. Once installed it makes two shortcuts in the Desktop : `run_TGF_TEB_propa_noGui` to launch the program with default parameters, and `TGF_TEB_Propa_output_dir` to open the output folder where output files (particle lists) are produced. You can edit `[INSTALL DIR]\dsarria\TGF_TEB_propa\Bin\run_windows_no_gui.bat` to change the simulation parameters (record altitude, TGF altitude, latitude, longitude, beaming, etc...).
- It is not recommended to install it inside `C:\Program Files (x86)\`  because in this case, admin rights will be required to modify the installed files. You may install inside `C:\Dev\` for example.

## Simulation Settings:
Most of settings can be adjusted in `include/Settings.hh`. In particular:
- `settings->RECORD_ELEC_POSI_ONLY` = boolean to record only electrons and positrons (not recording photons)
- `settings->RECORD_PHOT_ONLY` = boolean to record only photons (not recording electrons and positrons)
- If both are set to `false`, all particles are recorded.
- `settings->SOURCE_ALT` = TGF source altitude in km
- `settings->SOURCE_LAT` = TGF source latitude in deg
- `settings->SOURCE_LONG` = TGF source longitude in deg
- `settings->OPENING_ANGLE` = half-cone TGF opening angle in degrees. If "Gaussian" is selected for `settings->BEAMING_TYPE`, it is the sigma of the gaussian distribution.
- `settings->TILT_ANGLE` = TGF tilt angle in degrees
- `settings->BEAMING_TYPE` = TGF beaming type, that is a string that values "Uniform" or "Gaussian" for isotropic or gaussian distribution
- `settings->SOURCE_SIGMA_TIME` = TGF sigma time. Assumes the TGF has an intrinsic duration, that has Gaussian (=normal) distribution. The parameter is the sigma of this distribution, in microseconds
### Other settings:
- Two modes are possible: `visualization` and `run`. `visualization` will show the 3D geometry (simplified Earth) and particle track. `run` will not show any 3D visualization, to run the code as quickly as possible. By default, the mode is set to `visu` if no input argument for the executable is specified and `run` otherwise. This can be changed by editing the `G4String` variable `Mode` in the main function located in the source file `src/tgf_propa.cc`, that can be set to `"visu"` or `"run"`.
- Primary Generator is a point source, with adjustable altitude and geometry. See `src/src/PrimaryGeneratorAction.hh` and `src/src/PrimaryGeneratorAction.cc`
- The simulation stops when the number of recorded particles has reached `nb_to_get_per_run` (located in `tgf_propa.cc`), that can be changed in the source file or specified as input parameter of the executable.
- Atmosphere density is not constant with altitude, it evolves ~exponentially. However, Geant4 can only handle volumes with constant density, therefore the atmosphere is simulated by 256 exponentially spaced layers, each with constant density, covering altitude from 1 km to 150 km (negligible above). This can be changed in the source code, with the `src/src/DetectorConstruction.hh` and `src/src/DetectorConstruction.cc` files.
- If required, Magnetic Field can be turned ON with the Setting : `settings->MAG_FIELD_ON` set to `true`. Magnetic field is always turned OFF below 45 km altitude (where it is negligible), for performance.
- `settings->MAGNETIC_FIELD_MODEL` = magnetic field model to use. `"IGRF"` or `"WMM"`. On Windows, only WMM works.

## The program executable can accept input parameters in this order: 
- `nb_to_get_per_run` = number of particles to record before stopping run
- `settings->SOURCE_ALT` = TGF source altitude in km
- `settings->SOURCE_LAT` = TGF source latitude in deg
- `settings->SOURCE_LONG` = TGF source longitude in deg
- `settings->SOURCE_SIGMA_TIME` = TGF sigma time. Assumes the TGF has an intrinsic duration, that has Gaussian (=normal) distribution. The parameter is the sigma of this distribution, in microseconds
- `settings->OPENING_ANGLE` = half-cone TGF opening angle in degrees. If "Gaussian" is selected for `settings->BEAMING_TYPE`, it is the sigma of the gaussian distribution.
- `settings->TILT_ANGLE` = TGF tilt angle in degrees
- `settings->BEAMING_TYPE` = TGF beaming type, that is a string that values "Uniform" or "Gaussian" for isotropic or gaussian distribution
- `settings->record_altitude` = record altitude (in km) of the TGF, default is 400 km.

## Additional information:
- Recorded particles are outputed as a list (one by one) in files located in the `/output/` dir (`build/output/` on Linux). See `README_output.txt` to find what quantity is in what column.
- By default, the code uses the `G4EmStandardPhysics_option1` physics list, which is fast and accurate enough for this problem. This can be changed inside the source file `src/src/PhysicsList.cc`.
- The python script `build/run_on_multiple_cpu.py` makes it possible to run the code on multiple threads (CPU cores) by running several times the executable (possibly with different settings). Implementation is straightforward since every initial particle is independent. See comments inside the file. It requires `mpi4py`, `numpy`, and possibly other python libraries. Communication between python script and executable is done with the help of the parameters `int argv` and `char** argc`  of the main function in `src/tgf_propa.cc`.
- The code is made so that each run will have a different random seed (that is a `long` integer storing the current time given by the `std::chrono::high_resolution_clock` function, in nanoseconds). It assumes that the program cannot be launched twice during the exact same nanosecond.

## ToDo / 
- Improve 3D view to have more relevant automatic viewing angle
