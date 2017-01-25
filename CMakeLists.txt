cmake_minimum_required(VERSION 2.6)
project(moldycat)
# The version number.
set(moldycat_VERSION_MAJOR 0)
set(moldycat_VERSION_MINOR 8)


# CMake scripts MacroEnsureVersion and OptimizeForArchitecture were reproduced with permission from Freddie Witherden's teatree project 
# https://gitweb.freddie.witherden.org/?p=teatree.git;a=tree


set(CMAKE_MODULE_PATH "${moldycat_SOURCE_DIR}/cmake")
INCLUDE(MacroEnsureVersion)
INCLUDE(OptimizeForArchitecture)
include_directories(${Armadillo_INCLUDE_DIRS})
INCLUDE(FindArmadillo)
OptimizeForArchitecture()



IF(CMAKE_BUILD_TYPE MATCHES "Debug")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -DDBG")
set(MOLDY_EXE_NAME moldycat_dbg)
ELSE()
set(MOLDY_EXE_NAME moldycat)
ENDIF()

  
IF(CMAKE_BUILD_TYPE MATCHES "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}  -O3 -DARMA_NO_DEBUG -ffast-math -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Wall  -fopenmp")
ENDIF()

FIND_PACKAGE(OpenMP)

if(OPENMP_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

FIND_PACKAGE(Armadillo REQUIRED)

#SET(Armadillo_VERSION_OK TRUE)
#MACRO_ENSURE_VERSION("2.4.3" Armadillo_VERSION_STRING Armadillo_VERSION_OK)
if(Armdadilo_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Armadillo_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Armadillo_EXE_LINKER_FLAGS}")
INCLUDE_DIRECTORIES($(Armadillo_INCLUDE_DIR))
endif()


add_executable(${MOLDY_EXE_NAME} 
    src/moldycat.cpp 
    src/modules/moduleA.hpp
    src/modules/moduleB.hpp
    src/modules/moduleC.hpp
    src/modules/moduleDnew.hpp
    src/modules/moduleE.hpp
    src/modules/moduleG.hpp
    src/atom/atom.cpp
    src/atom/vect3.cpp
    src/sim/unitcell.cpp
    src/sim/simulation.hpp
    src/sim/forceeval.hpp
    src/stepper/stepperbase.hpp
    src/stepper/verlet.hpp
    src/stepper/sprk.hpp
    src/stepper/factory.hpp
    src/thermo/rand.cpp
    src/thermo/thermobase.cpp
    src/thermo/velocityscaler.cpp
    src/thermo/andersen.cpp
    src/thermo/nothermo.cpp
    src/opt/optbase.hpp
    src/opt/dmdcm.hpp
    src/opt/dmdum.hpp
    src/opt/cg.hpp
    src/opt/sd.hpp
    src/opt/factory.hpp
    src/anal/radial.hpp
    src/anal/elasticConstants.hpp
    src/anal/graphiteGenerator.hpp
    src/anal/evStretcher.hpp)


target_link_libraries(${MOLDY_EXE_NAME} armadillo)


