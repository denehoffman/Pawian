add_library(MinFunctions
  AbsFcn.cc
  AbsPawianMinimizer.cc
  EvoMinimizer.cc
  MinuitMinimizer.cc
  PwaFcnBase.cc
  PwaFcnServer.cc
)

target_include_directories(MinFunctions
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(MinFunctions
  Utils
  qft++
  ErrLogger  
  ConfigParser
  FitParams
  PwaDynamics
  PspGen
  PwaUtils
  ROOT::Minuit2
  ${Boost_LIBRARIES}
)

install(TARGETS MinFunctions
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)