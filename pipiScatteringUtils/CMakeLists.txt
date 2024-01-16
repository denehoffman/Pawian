add_library(pipiScatteringUtils
  PiPiScatteringChannelEnv.cc
  pipiScatteringBaseLh.cc
)

target_include_directories(pipiScatteringUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(pipiScatteringUtils  
  Utils
  qft++
  ErrLogger
  ConfigParser
  Event
  PwaUtils
  Particle
)

install(TARGETS pipiScatteringUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)