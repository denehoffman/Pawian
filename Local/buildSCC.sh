#!/bin/bash

# Set some variables
export ROOTSYS=/opt/root
export CMAKE_SOURCE_DIR=/home/rberlich/Pawian
export LD_LIBRARY_PATH=${ROOTSYS}/lib
export PATH=${PATH}:${ROOTSYS}/bin
export BOOST_BUILD_PATH=${CMAKE_SOURCE_DIR}/build/
export Minuit2_ROOT=/opt/minuit2
export GENEVA_ROOT=/opt/geneva
export BUILD_TYPE=Debug
export VERBOSE=1
export PAWIAN_BOOSTROOT=/opt/boost146

cmake -DCMAKE_VERBOSE_MAKEFILE=${VERBOSE} -DPAWIAN_BOOSTROOT=${PAWIAN_BOOSTROOT} -DPAWIAN_BUILD_TYPE=${BUILD_TYPE} ..
