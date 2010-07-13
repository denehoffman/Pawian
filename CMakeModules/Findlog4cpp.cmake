# If log4cpp_ROOT was defined in the environment, use it.
if(NOT log4cpp_ROOT AND NOT $ENV{log4cpp_ROOT} STREQUAL "")
    set(log4cpp_ROOT $ENV{log4cpp_ROOT})
endif()

if(log4cpp_ROOT)
    set(_log4cpp_INCLUDE_SEARCH_DIRS
        ${log4cpp_INCLUDEDIR}
        ${log4cpp_ROOT}/include
        ${log4cpp_ROOT}/inc
        ${log4cpp_ROOT})
    set(_log4cpp_LIBRARIES_SEARCH_DIRS
        ${log4cpp_LIBRARYDIR}
        ${log4cpp_ROOT}/lib
        ${log4cpp_ROOT}/src/.libs
        ${log4cpp_ROOT})
endif()

find_path(log4cpp_INCLUDE_DIRS log4cpp/Category.hh
    HINTS ${_log4cpp_INCLUDE_SEARCH_DIRS})

set(log4cpp_LIBR log4cpp)
find_library(log4cpp_LIBRARY log4cpp
    HINTS ${_log4cpp_LIBRARIES_SEARCH_DIRS})

set(log4cpp_LIBRARIES ${log4cpp_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(log4cpp DEFAULT_MSG
    log4cpp_INCLUDE_DIRS log4cpp_LIBRARY)

get_filename_component(log4cpp_LIBRARY_DIRS ${log4cpp_LIBRARY} PATH)

mark_as_advanced(
    log4cpp_INCLUDE_DIRS
    log4cpp_LIBRARY_DIRS
    log4cpp_LIBRARIES)
