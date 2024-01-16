add_library(ConfigParser  
  epemParser.cc
  gammapParser.cc
  ggParser.cc
  globalParser.cc
  KMatrixParser.cc
  ParserBase.cc
  pbarpParser.cc
  pipiScatteringParser.cc
  resParser.cc  
)

target_include_directories(ConfigParser
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(ConfigParser  
  ErrLogger
  ${Boost_LIBRARIES}    
)

install(TARGETS ConfigParser
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
