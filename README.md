TGF-Propagation-Geant4
=======
Geant4 model for Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons propagation in Earth atmosphere and environment (magnetic field).
=======

contact : <david.sarria@uib.no>

## Generalities
- Propagation of photons, electrons and positron in Earth's environment (atmosphere, ionosphere, magnetosphere), in the context of Terrestrial Gamma-ray Flashes (TGF) and associated electrons and positrons beams.
- This code is a first test and is far from being perfect. Feel free to clone and suggest improvements for the code.
- Uses mostly Geant4 features. See [documentation](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/index.html "Geant4 documentation").
- Integrates the [NRL-MSISE-00 model](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/) for the atmosphere and the [IGRF-12](http://wdc.kugi.kyoto-u.ac.jp/igrf/index.html) model for the magnetic field.

## Compilation, installation
- The source code of `TGF-Propagation-Geant4` is localted in `src/` and the build should be done in the folder `build/`.
- Requires [Geant4](https://geant4.web.cern.ch/) compiled and [installed](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html) properly. Recommended to use geant4-10-04-patch-02 (25-May-2018). Minor changes in the source code may be required for other versions. Geant4 source code is included in folder `geant4/geant4_source/`.

For Ubuntu users, an installation script is provided:
- it will compile and install Geant4 [from source](https://geant4.web.cern.ch/node/1604), set up the environement variables and compile the `TGF-Propagation-Geant4` code in the `build/` folder.
- Go to `TGF-Propagation-Geant4` directory and run in terminal `bash compile_install.bash`.
- Will require super user priviledges (`sudo`) to download missing dependencies. 
- alternatively, run command `sudo apt-get install build-essential qt4-default qtcreator cmake-qt-gui gcc g++ gfortran zlib1g-dev libxerces-c-dev libx11-dev libexpat1-dev libxmu-dev libmotif-dev libboost-filesystem-dev libeigen3-dev qt4-qmake` to install dependencies before-hand, and `sudo` should not be required.
- It was successfully tested on Ubuntu 16.04 and 18.04, but is probably not free of bugs.

For non-Ubuntu users look at the code in the file `compile_install.bash` to check the compile and set-up steps, and read the [Geant4 installation instructions](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html)

## Simulation Settings:
- Most of settings can be adjusted in `src/Settings.cc`
- Two mode are possible "visualization" and "run". "visualization" will show the 3D geometry (simplified Earth) and particle track. "run" will not show any 3D visualization, to run the code as quickly as possible. This can be changed by editing the `G4String` variable `Mode` in the main function located in the source file `src/tgf_propa.cc`, that can be set to `"visu"` or `"run"`.
- Primary Generator is a point source, with adjustable altitude and geometry. See `src/src/PrimaryGeneratorAction.hh` and `src/src/PrimaryGeneratorAction.cc`
- Record is made in a layer at chosen altitude(s). Record altitudes have to be set at the beginning of the main function inside `src/tgf_propa.cc`, e.g. `Settings::record_altitudes.push_back(400.);` (input is altitude in km)
- The simulation stops when the number of recorded particles has reached `nb_to_get_per_run`, that can be changed.
- Atmosphere density is not constant with altitude, it evolves ~exponentially. However, Geant4 can only handle volumes with constant density, therefore the atmosphere is simulated by 256 exponentially spaced layers, each with constant density, covering altitude from 1 km to 150 km (negligible above). This can be changed in the source code, with the `src/src/DetectorConstruction.hh` and `src/src/DetectorConstruction.cc` files.
- Recorded particles are outputed as a list (one by one) in files located in `build/output/`. See `src/src/Analysis.cc` to find which quantity is in which column.
- If required, Magnetic Field can be turned ON with the Setting : `Settings::MAG_FIELD_ON` set to `true`. Magnetic field is always turned OFF below 45 km altitude (where it is negligible), for performance. The model of magnetic field is [IGRF-12](http://wdc.kugi.kyoto-u.ac.jp/igrf/index.html).

The code is built so that the executable can accept input parameters: 
- `Settings::SOURCE_ALT` = TGF source altitude in km
- `Settings::OPENING_ANGLE` = TGF beaming angle on degrees. If "Gaussian" is selected for `Settings::BEAMING_TYPE`, it is the sigma of the gaussian
- `Settings::TILT_ANGLE` = TGF tilt angle in degrees
- `Settings::BEAMING_TYPE` = TGF beaming type, that is a string that values "Uniform" or "Gaussian" for isotropic or gaussian distribution
- `Settings::SOURCE_SIGMA_TIME` = TGF sigma time. Assumes the TGF has an intrinsic duration, that has Gaussian (=normal) distribution. The parameter is the sigma of this distribution, in microseconds

- The python script `build/run_on_multiple_cpu.py` makes it possible to run the code on multiple threads (CPU cores) by running several times the executable (possibly with different settings). Implementation is straightforward since every initial particle is independent. See comments inside the file. It requires `mpi4py`, `numpy`, and possibly other python libraries. Communication between python script and executable is done with the help of the parameters `int argv` and `char** argc`  of the main function in `src/tgf_propa.cc`.
- The code is made so that each run will have a different random seed (that is a `long` integer storing the current time given by the `std::chrono::high_resolution_clock` function, in nanoseconds). It assumes that the program cannot be launched twice during the exact same nanosecond.

- By default, it uses the `G4EmStandardPhysics_option1` physics list, which is fast and accurate enough for this problem. This can be changed inside the source file `src/src/PhysicsList.cc`.

## ToDo / 

- Magnetic field track of particles (using IGRF) is not working properly yet, by comparison to other reference model (MCPEP).
- Improve 3D view to have more relevant automatic viewing angle
