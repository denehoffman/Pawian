#!/bin/bash

# Set some variables
export ROOTSYS=/opt/root
export CMAKE_SOURCE_DIR=/home/rberlich/Pawian
export LD_LIBRARY_PATH=${ROOTSYS}/lib
export PATH=${PATH}:${ROOTSYS}/bin
export BOOST_BUILD_PATH=${CMAKE_SOURCE_DIR}/build/
export BOOST_ROOT=/opt/boost146
export BOOST_INCLUDEDIR=/opt/boost146/include
export BOOST_LIBRARYDIR=/opt/boost146/lib
export Minuit2_ROOT=/opt/minuit2
export GENEVA_ROOT=/opt/geneva

cmake ..
