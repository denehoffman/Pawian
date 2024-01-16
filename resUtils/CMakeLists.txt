file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(resUtils
  ${sourcefiles} 
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