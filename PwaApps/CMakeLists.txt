add_library(PwaApps  
  Dummy.cc  
)

target_include_directories(PwaApps
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    #${ROOT_INCLUDE_DIRS}
    #${GSL_INCLUDE_DIRS}  
)

#target_link_libraries(PwaApps)

install(TARGETS PwaApps
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(coupledChannelApp
  coupledChannelApp.cc
)

target_include_directories(coupledChannelApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(coupledChannelApp
  PwaApps
  pbarpUtils
  epemUtils
  ConfigParser
  FitParams
  PwaUtils
  AppUtils
  Particle
  Event
  qft++  
  ErrLogger
  ${Boost_LIBRARIES}  
)

add_executable(singleChannelApp
singleChannelApp.cc
)

target_include_directories(singleChannelApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/
    ${Boost_INCLUDE_DIRS}  
)

target_link_libraries(singleChannelApp
  PwaApps
  pbarpUtils
  epemUtils
  ConfigParser
  FitParams
  PwaUtils
  AppUtils
  Utils
  Particle
  Event
  qft++  
  ErrLogger
  ${Boost_LIBRARIES} 
)

install(TARGETS coupledChannelApp singleChannelApp
  RUNTIME DESTINATION bin
)