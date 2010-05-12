# If Minuit2_ROOT was defined in the environment, use it.
if(NOT Minuit2_ROOT AND NOT $ENV{Minuit2_ROOT} STREQUAL "")
    set(Minuit2_ROOT $ENV{Minuit2_ROOT})
endif()

if(Minuit2_ROOT)
    set(_Minuit2_INCLUDE_SEARCH_DIRS
        ${Minuit2_INCLUDEDIR}
        ${Minuit2_ROOT}/include
        ${Minuit2_ROOT}/inc
        ${Minuit2_ROOT})
    set(_Minuit2_LIBRARIES_SEARCH_DIRS
        ${Minuit2_LIBRARYDIR}
        ${Minuit2_ROOT}/lib
        ${Minuit2_ROOT}/src/.libs
        ${Minuit2_ROOT})
endif()

find_path(Minuit2_INCLUDE_DIRS Minuit2/FCNBase.h
    HINTS ${_Minuit2_INCLUDE_SEARCH_DIRS})

set(Minuit2_LIBRARIES Minuit2)
find_library(Minuit2_LIBRARY ${Minuit2_LIBRARIES}
    HINTS ${_Minuit2_LIBRARIES_SEARCH_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Minuit2 DEFAULT_MSG
    Minuit2_INCLUDE_DIRS Minuit2_LIBRARY)

# If Minuit2 is built with GCC and OpenMP enabled, it contains
# undefined references to some omp_* and GOMP_* symbols. 3rd-party
# programs therefore need to link to gomp on their own.
# if(UNIX)
#     execute_process(COMMAND ldd -r ${Minuit2_LIBRARY}
#         ERROR_VARIABLE Minuit2_LDD_ERROR OUTPUT_QUIET)
#     if (Minuit2_LDD_ERROR MATCHES "undefined.*GOMP")
#         set(Minuit2_LIBRARIES ${Minuit2_LIBRARIES} gomp)
#     endif()
# endif()

get_filename_component(Minuit2_LIBRARY_DIRS ${Minuit2_LIBRARY} PATH)

mark_as_advanced(
    Minuit2_INCLUDE_DIRS
    Minuit2_LIBRARY_DIRS
    Minuit2_LIBRARIES)

