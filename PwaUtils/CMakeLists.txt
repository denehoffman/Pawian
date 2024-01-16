file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(PwaUtils
  ${sourcefiles}
)

target_include_directories(PwaUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(PwaUtils
  Utils
  qft++
  qft++Extension
  ErrLogger
  Event
  FitParams
  ConfigParser
  PwaDynamics
  PspGen
)

install(TARGETS PwaUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
