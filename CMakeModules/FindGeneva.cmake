# If Minuit2_ROOT was defined in the environment, use it.
if(NOT GENEVA_ROOT AND NOT $ENV{GENEVA_ROOT} STREQUAL "")
    set(GENEVA_ROOT $ENV{GENEVA_ROOT})
endif()

if(GENEVA_ROOT)
    set(_GENEVA_INCLUDE_SEARCH_DIRS
        ${GENEVA_INCLUDEDIR}
        ${GENEVA_ROOT}/include
        ${GENEVA_ROOT}/inc
        ${GENEVA_ROOT})
    set(_GENEVA_LIBRARIES_SEARCH_DIRS
        ${GENEVA_LIBRARYDIR}
        ${GENEVA_ROOT}/lib
        ${GENEVA_ROOT}/src/.libs
        ${GENEVA_ROOT})
endif()

find_path(GENEVA_INCLUDE_DIRS geneva/GIndividual.hpp
    HINTS ${_GENEVA_INCLUDE_SEARCH_DIRS})

find_library(GENEVA_LIBRARY gemfony-geneva
    HINTS ${_GENEVA_LIBRARIES_SEARCH_DIRS})
find_library(GENEVA_IND_LIBRARY gemfony-geneva-individuals
    HINTS ${_GENEVA_LIBRARIES_SEARCH_DIRS})
find_library(GENEVA_COM_LIBRARY gemfony-common
    HINTS ${_GENEVA_LIBRARIES_SEARCH_DIRS})
find_library(GENEVA_HAP_LIBRARY gemfony-hap
    HINTS ${_GENEVA_LIBRARIES_SEARCH_DIRS})
find_library(GENEVA_COU_LIBRARY gemfony-courtier
    HINTS ${_GENEVA_LIBRARIES_SEARCH_DIRS})

set(GENEVA_LIBRARIES 
  ${GENEVA_LIBRARY}
  ${GENEVA_IND_LIBRARY}
  ${GENEVA_COM_LIBRARY} 
  ${GENEVA_HAP_LIBRARY} 
  ${GENEVA_COU_LIBRARY} 
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GENEVA DEFAULT_MSG
    GENEVA_INCLUDE_DIRS GENEVA_LIBRARY)

# If Minuit2 is built with GCC and OpenMP enabled, it contains
# undefined references to some omp_* and GOMP_* symbols. 3rd-party
# programs therefore need to link to gomp on their own.
#if(UNIX)
#     execute_process(COMMAND ldd -r ${Minuit2_LIBRARY}
#         ERROR_VARIABLE Minuit2_LDD_ERROR OUTPUT_QUIET)
#     if (Minuit2_LDD_ERROR MATCHES "undefined.*GOMP")
#         set(Minuit2_LIBRARIES ${Minuit2_LIBRARIES} gomp)
#     endif()
#endif()

get_filename_component(GENEVA_LIBRARY_DIRS ${GENEVA_LIBRARY} PATH)

mark_as_advanced(
    GENEVA_INCLUDE_DIRS
    GENEVA_LIBRARY_DIRS
    GENEVA_LIBRARIES)

