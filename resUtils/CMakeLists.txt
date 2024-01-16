add_library(resUtils
  ResChannelEnv.cc
  res1DMassShapeLh.cc
  resBaseLh.cc
  resReaction.cc 
)

target_include_directories(resUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(resUtils  
  Utils
  qft++
  ErrLogger
  ConfigParser
  PwaUtils
  Particle
)

install(TARGETS resUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)