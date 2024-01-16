add_library(pbarpUtils
  PbarpChannelEnv.cc
  pbarpBaseLh.cc
  pbarpCanoLh.cc
  pbarpHeliLh.cc
  pbarpReaction.cc
  pbarpStatesLS.cc
  pbarpTensorLh.cc
  spinDensityHist.cc
)

target_include_directories(pbarpUtils
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(pbarpUtils  
  Utils
  qft++
  ErrLogger
  ConfigParser
  FitParams
  PwaUtils
)

install(TARGETS pbarpUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)