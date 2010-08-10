# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

# Set name of our project to "PAWIAN" (= PArtial Wave Interactive ANalysis)                  
project(PAWIAN)


set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

# Load some basic macros which are needed later on
include(WriteConfigFile)

# searches for needed packages
# REQUIRED means that cmake will stop if this packages are not found
# For example the framework can run without GEANT4, but ROOT is
# mandatory
find_package(ROOT REQUIRED)

SET (Boost_USE_MULTITHREAD ON)

SET (Boost_ADDITIONAL_VERSIONS "1.42" "1.43.0")

FIND_PACKAGE(Boost 1.42 REQUIRED COMPONENTS
        date_time
        filesystem
        program_options
        regex
        serialization
        system
        thread
)


# The test component is not necessary for the library, but needed for the tests
FIND_PACKAGE(Boost 1.42 COMPONENTS
        test_exec_monitor
)

MESSAGE("Boost_INCLUDE_DIRS: "  ${Boost_INCLUDE_DIRS})
MESSAGE("Boost_LIBRARIES: "  ${Boost_LIBRARIES})
MESSAGE("Boost_LIBRARY_DIRS: "  ${Boost_LIBRARY_DIRS})


FIND_PACKAGE(Minuit2 REQUIRED)
MESSAGE("Minuit2_INCLUDE_DIRS: " ${Minuit2_INCLUDE_DIRS})
MESSAGE("Minuit2_LIBRARY_DIRS: " ${Minuit2_LIBRARY_DIRS})
MESSAGE("Minuit2_LIBRARY: " ${Minuit2_LIBRARY})
MESSAGE("Minuit2_LIBRARIES: " ${Minuit2_LIBRARIES})

FIND_PACKAGE(log4cpp REQUIRED)
MESSAGE("log4cpp_INCLUDE_DIRS: " ${log4cpp_INCLUDE_DIRS})
MESSAGE("log4cpp_LIBRARY_DIRS: " ${log4cpp_LIBRARY_DIRS})
MESSAGE("log4cpp_LIBRARY: " ${log4cpp_LIBRARY})
MESSAGE("log4cpp_LIBRARIES: " ${log4cpp_LIBRARIES})

FIND_PACKAGE(Geneva REQUIRED)
MESSAGE("GENEVA_INCLUDE_DIRS: " ${GENEVA_INCLUDE_DIRS})
MESSAGE("GENEVA_LIBRARY_DIRS: " ${GENEVA_LIBRARY_DIRS})
MESSAGE("GENEVA_LIBRARIES: " ${GENEVA_LIBRARIES})
MESSAGE("GENEVA_LIBRARY: " ${GENEVA_LIBRARY})
# Set the library type to static by default
IF( NOT PAWIAN_STATIC )
        SET( PAWIAN_STATIC TRUE )
ENDIF()

SET (CMAKE_BUILD_TYPE Debug)

IF( PAWIAN_STATIC )
#        SET (CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -lpthread")
#        SET (LINK_FLAGS "-static-libgcc -lpthread")
	 SET (CMAKE_EXE_LINKER_FLAGS "-lpthread -Wall")
	 SET (LINK_FLAGS "-lpthread -Wall")
        SET (BUILD_SHARED_LIBS OFF)
        SET (Boost_USE_STATIC_LIBS ON)
        LINK_LIBRARIES (dl)
#        LINK_LIBRARIES (z)      
        LINK_LIBRARIES (pthread)
ELSE()
        SET (BUILD_SHARED_LIBS ON)
        SET (Boost_USE_STATIC_LIBS OFF)
        LINK_LIBRARIES (pthread)
ENDIF()

ADD_SUBDIRECTORY (CMakeModules)
add_subdirectory (qft++)
add_subdirectory (Utils)
add_subdirectory (ErrLogger)
add_subdirectory (Particle)
add_subdirectory (DecayTree)
add_subdirectory (Event)
add_subdirectory (Setup)
add_subdirectory (PwaUtils)
add_subdirectory (Examples/qft++)
add_subdirectory (Examples/pbarp)
#add_subdirectory (Examples/EtacToapi0Fit)
add_subdirectory (Examples/geneva/GMultiPopulation)
add_subdirectory (Examples/geneva/GStartProject)
#add_subdirectory (Examples/pbarpToOmegaPi)
add_subdirectory (Examples/MATpbarpToOmegaPi)


EXECUTE_PROCESS(COMMAND ps OUTPUT_VARIABLE PSOUTPUT)

IF(${PSOUTPUT} MATCHES tcsh )
 SET(CONFIG_TCSH_FILE CONFIG_TCSH_FILE-NOTFOUND)
 FIND_FILE(CONFIG_TCSH_FILE config.csh PATHS ${CMAKE_BINARY_DIR}) 
 IF(CONFIG_TCSH_FILE MATCHES config.csh) 
  EXECUTE_PROCESS(COMMAND rm config.csh)
 ENDIF() 
 WRITE_CONFIG_FILE(config.csh)
ELSEIF(${PSOUTPUT} MATCHES bash)
 SET(CONFIG_SH_FILE CONFIG_SH_FILE-NOTFOUND)
 FIND_FILE(CONFIG_SH_FILE config.sh PATHS ${CMAKE_BINARY_DIR}) 
 IF(CONFIG_SH_FILE MATCHES config.sh) 
  EXECUTE_PROCESS(COMMAND rm config.sh)
 ENDIF() 
  WRITE_CONFIG_FILE(config.sh)
ELSE()
 MESSAGE("No config file: tcsh or bash not found !!!")
ENDIF()

 


