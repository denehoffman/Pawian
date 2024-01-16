file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(pipiScatteringUtils
  ${sourcefiles}
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