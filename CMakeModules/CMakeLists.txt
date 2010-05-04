# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

SET ( CMAKEFILES
        FindROOT.cmake
)

INSTALL ( FILES ${CMAKEFILES} DESTINATION CMakeModules )
