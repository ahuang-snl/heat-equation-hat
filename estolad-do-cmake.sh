#!/bin/bash


# Note that in this MixedPoisson example (built against Trilinos), the include "Phalanx_ConfigDefs.hpp" is replaced with "Phalanx_config.hpp"

rm -rf CMakeCache.txt CMakeFiles src test

module purge

#module load  gcc/4.8.4/openmpi/1.8.7
module load sems-gcc/4.8.4 sems-openmpi/1.8.7
module load sems-cmake/2.8.12

module load sems-gcc/4.8.4 sems-openmpi/1.8.7 sems-cmake/2.8.12 boost/1.59.0/gcc/4.8.4/base netcdf/4.3.2/gcc/4.8.4/openmpi/1.8.7 hdf5/1.8.12/gcc/4.8.4/openmpi/1.8.7 sems-numpy/1.9.1/base

#env | grep SEMS >& modules-loaded.txt

#-D CMAKE_CXX_COMPILER:FILEPATH="/projects/sems/install/rhel6-x86_64/sems/compiler/gcc/4.8.4/openmpi/1.8.7/bin/mpic++" \
#\
#-D CMAKE_CXX_FLAGS="-g -std=c++11 -Wall" \
#-D CMAKE_C_FLAGS="-g -Wall" \
#-D CMAKE_Fortan_FLAGS="-g"\
#-D CMAKE_C_COMPILER:FILEPATH="/projects/sems/install/rhel6-x86_64/sems/compiler/gcc/4.8.4/openmpi/1.8.7/bin/mpicc" \
#-D CMAKE_Fortran_COMPILER:FILEPATH="/projects/sems/install/rhel6-x86_64/sems/compiler/gcc/4.8.4/openmpi/1.8.7/bin/mpif77" \
#\


TRILINOS_INSTALL_DIR=/home/ahuang/playground/install-trilinos/

cmake \
\
  -D Trilinos_INSTALL_DIR=${TRILINOS_INSTALL_DIR} \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
\
  ..