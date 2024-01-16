file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(qaErrorExtract  
  ${sourcefiles} 
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