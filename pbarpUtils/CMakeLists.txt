file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(pbarpUtils
  ${sourcefiles}
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