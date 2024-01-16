add_library(qaErrorExtract  
  AbsCovMatErrorExtract.cc
  ScatteringContribution.cc
  WaveContribution.cc  
)

target_include_directories(qaErrorExtract
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(qaErrorExtract  
  ErrLogger
  ConfigParser
  FitParams
  PwaUtils
  pipiScatteringUtils
)

install(TARGETS qaErrorExtract
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)