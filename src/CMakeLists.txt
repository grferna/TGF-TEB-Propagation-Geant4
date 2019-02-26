#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TGF_propa C CXX Fortran)

execute_process(COMMAND "source ~/.bashrc")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
endif ()

#set(CMAKE_BUILD_TYPE Release)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
        ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB sourcesS ${PROJECT_SOURCE_DIR}/src/steppers/*.cc)
file(GLOB sourcesSH ${PROJECT_SOURCE_DIR}/include/steppers/*.hh)
file(GLOB sourcesBL ${PROJECT_SOURCE_DIR}/src/Bline_tracer/*.cc)
file(GLOB sourcesBLH ${PROJECT_SOURCE_DIR}/include/Bline_tracer/*.hh)

file(GLOB sourcesc ${PROJECT_SOURCE_DIR}/src/C/*.c)
file(GLOB sourcesh ${PROJECT_SOURCE_DIR}/include/C/*.h)
file(GLOB sourcesf ${PROJECT_SOURCE_DIR}/src/fortran/*.f)
file(GLOB sourcesf2 ${PROJECT_SOURCE_DIR}/src/msis/*.for)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(TGF_Propa tgf_propa.cc ${sourcesf} ${sourcesf2} ${sourcesc} ${sourcesh} ${sources} ${headers} ${sourcesS} ${sourcesSH} ${sourcesBL} ${sourcesBLH})
target_link_libraries(TGF_Propa ${Geant4_LIBRARIES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    target_compile_options(TGF_Propa PRIVATE -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    target_compile_options(TGF_Propa PRIVATE -fp-model=extended)
endif ()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(TGF_propa DEPENDS TGF_Propa)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TGF_Propa DESTINATION bin)

