add_library(epemUtils
  EpemChannelEnv.cc
  epemBaseLh.cc
  epemCanoLh.cc
  epemHeliLh.cc
  epemReaction.cc
  epemTensorLh.cc
)

target_include_directories(epemUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(epemUtils  
  Utils
  qft++
  ErrLogger
  ConfigParser
  PwaUtils
)

install(TARGETS epemUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)