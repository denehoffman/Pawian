cmake_minimum_required(VERSION 3.27)
project(Pawian)

set(CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

find_package(Boost)
find_package(Boost COMPONENTS program_options unit_test_framework serialization filesystem timer REQUIRED)
find_package(ROOT)
find_package(GSL REQUIRED)

#include_directories(
# ${Boost_INCLUDE_DIRS}
# ${ROOT_INCLUDE_DIRS}  
#)

#link_libraries(
# ${Boost_LIBRARIES}
# ${ROOT_LIBRARIES}
#)

add_subdirectory(ErrLogger)
add_subdirectory(qft++)
add_subdirectory(qft++Extension)
add_subdirectory(ConfigParser)
add_subdirectory(Particle)
add_subdirectory(Utils)
add_subdirectory(Event)
add_subdirectory(FitParams)
add_subdirectory(PwaDynamics)
add_subdirectory(PwaUtils)
add_subdirectory(MinFunctions)
add_subdirectory(AppUtils)
add_subdirectory(pbarpUtils)
add_subdirectory(epemUtils)
add_subdirectory(resUtils)
add_subdirectory(ggUtils)
add_subdirectory(pipiScatteringUtils)
add_subdirectory(gslUtils)
add_subdirectory(qaErrorExtract)
add_subdirectory(PspGen)
add_subdirectory(KMatrixExtract)
add_subdirectory(PwaApps)
