file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(epemUtils
  ${sourcefiles}
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