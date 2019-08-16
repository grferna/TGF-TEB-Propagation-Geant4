#!/bin/bash
set -e

echo "checking if OS is Ubuntu..."
# Checking if OS is Ubuntu
if [ -f /etc/os-release ]; then
  . /etc/os-release
  OS=$NAME
  VER=$VERSION_ID
fi

if [ ! "$OS" = "Ubuntu" ]; then
  echo "Error: OS is not Ubuntu. Script works only for Ubuntu. Aborting."
  exit 1
else
  echo "... OS is Ubuntu"
fi

echo "checking if current folder name is "TGF-TEB-Propagation-Geant4" ..."
# checking if current folder name is " TGF-TEB-Propagation-Geant4 "
result=${PWD##*/}
if [ ! [$result = "TGF-TEB-Propagation-Geant4" || $result = "TGF-TEB-Propagation-Geant4-master" ] ]; then
  echo "Error: Current folder should be TGF-TEB-Propagation-Geant4. The bash script compile_install.sh is made to run from it. Aborting."
  exit 1
else
  echo "...current folder name is TGF-TEB-Propagation-Geant4"
fi

############# VARIABLES
current_dir=$PWD
srcdir=$current_dir/geant4/geant4_source/
builddir=$current_dir/geant4/geant4_build/
installdir_g4=$current_dir/geant4/geant4_install/
datadir=$current_dir/geant4/data/

datadir2=$current_dir/geant4/data/

build_dir_tgfprop=$current_dir/build

geant4_lib_dir=$current_dir/geant4/geant4_install/lib/Geant4-10.4.2/

ubuntu_dependences_list=("build-essential"
  "qt4-default"
  "qtcreator"
  "cmake-qt-gui"
  "gcc"
  "g++"
  "gfortran"
  "zlib1g-dev"
  "libxerces-c-dev"
  "libx11-dev"
  "libexpat1-dev"
  "libxmu-dev"
  "libmotif-dev"
  "libboost-filesystem-dev"
  "libeigen3-dev"
  "qt4-qmake"
  "libtool"
  "m4"
  "automake"
)

entered_one_time=true

############# FUNCTIONS
run_install() {
  echo "Some missing dependencies were detected."
  ## Prompt the user
  if [ entered_one_time=true ]; then
    entered_one_time=false
    read -p "Do you have (root) sudo access ? [Y/n]. It is required to install missing dependencies: " answer
    ## Set the default value if no answer was given
    answer=${answer:N}
    if [[ $answer =~ [Nn] ]]; then
      echo "root access is required to install missing dependencies. Aborting."
      exit 1
    fi
  fi
  ## Prompt the user
  read -p "Do you want to install missing dependencies? [Y/n]: " answer
  ## Set the default value if no answer was given
  answer=${answer:Y}
  ## If the answer matches y or Y, install
  if [[ $answer =~ [Yy] ]]; then
    sudo apt-get install ${ubuntu_dependences_list[@]}
  else
    echo "Missing dependencies are required for proper compilation and installation. Aborting."
    exit 0
  fi
}

check_dependencies() {
  echo "checking dependencies..."

  dpkg -s "${ubuntu_dependences_list[@]}" >/dev/null 2>&1 || run_install

}

check_dependencies

echo "dependencies are OK"

##############

build_geant4() {

  cd ${builddir}

  echo "build_geant4: Attempt to execute CMake"

  env -i \
    QT_SELECT=4 \
    PATH=/usr/bin \
    ../../cmake/bin/cmake \
    -DCMAKE_INSTALL_PREFIX=$installdir_g4 \
    -DCMAKE_BUILD_TYPE=Release \
    -DGEANT4_BUILD_MULTITHREADED=OFF \
    -DGEANT4_BUILD_CXXSTD=c++11 \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_INSTALL_DATADIR=$datadir \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_USE_G3TOG4=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_FORCE_QT4=ON \
    -DGEANT4_USE_XM=ON \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_INVENTOR=OFF \
    -DGEANT4_USE_RAYTRACER_X11=ON \
    -DGEANT4_USE_SYSTEM_CLHEP=OFF \
    -DGEANT4_USE_SYSTEM_EXPAT=ON \
    -DGEANT4_USE_SYSTEM_ZLIB=OFF \
    -DCMAKE_INSTALL_LIBDIR=lib \
    ../geant4_source/

  echo "Attempt to compile Geant4"
  G4VERBOSE=1 make -j4
  make install

}

build_geant4

set_environement() {

  echo "Attempt to setup up environement variables..."

  # clean environement that was previously set by this script
  first_line=$(grep -n "## --> Added by G4-TGF-TEB-propa installation script" ~/.bashrc | awk -F ":" '{print $1}')
  echo $first_line
  last_line=$(grep -n "## <-- Added by G4-TGF-TEB-propa installation script" ~/.bashrc | awk -F ":" '{print $1}')
  echo $last_line

  re='^[0-9]+$'
  if [[ $first_line =~ $re ]]; then # if $first_line is a number (i.e. it was found)
    if [[ $first_line =~ $re ]]; then # if $last_line is a number (i.e. it was found)
      sed -i.bak "${first_line},${last_line}d" ~/.bashrc # delete text in .bashrc from first-line to last-line
    fi
  fi

  # set environement
  echo "## --> Added by G4-TGF-TEB-propa installation script" >>~/.bashrc
  echo " " >>~/.bashrc

  cd $current_dir
  cp ./geant4/geant4.sh ./geant4/geant4_install/bin/geant4.sh

  if grep -Fxq "source $current_dir/geant4/geant4_install/bin/geant4.sh" ~/.bashrc; then
    echo "< source $current_dir/geant4/geant4_install/bin/geant4.sh > already set up in ~/.bashrc"
  else
    echo "source $current_dir/geant4/geant4_install/bin/geant4.sh" >>~/.bashrc
  fi

  echo " " >>~/.bashrc
  echo "## <-- Added by G4-TGF-TEB-propa installation script" >>~/.bashrc

  echo -e "${RED}Please excecute command < ${GREEN}source ./geant4/geant4_install/bin/geant4.sh${RED} > for the system to be able to find the databases and libraries.${NC}"

}

set_environement

build_tgf_propa() {

  source ~/.bashrc

  cd ${build_dir_tgfprop}

  echo "build_tgf_propa : Attempt to execute CMake"

  ../cmake/bin/cmake \
    -DCMAKE_BUILD_TYPE=DEBUG \
    -DGeant4_DIR=$geant4_lib_dir \
    ../src/

  RED='\033[0;31m'
  GREEN='\033[0;32m'
  NC='\033[0m'

  echo "Attempt to compile TGF propa project..."
  G4VERBOSE=1 make -j4
  echo -e "...Done. Executable is ./build/TGF_Propa"

  echo "  "
  echo -e "${RED}Please excecute command < ${GREEN}source ./geant4/geant4_install/bin/geant4.sh${RED} > for the system to be able to find the databases and libraries. Alternatively, open a new terminal.${NC}"
  echo -e "${RED}The executable is ./build/TGF_propa .${NC}"

}

build_tgf_propa

## executing G4-TGF-TEB-propagation program
cd $current_dir
source ./geant4/geant4_install/bin/geant4.sh
cd build
./TGF_Propa
