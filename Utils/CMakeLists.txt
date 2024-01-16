file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(Utils  
  ${sourcefiles}
)

target_include_directories(Utils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(Utils  
  ErrLogger  
)

add_library(Utils_ErrValue
  ErrValue.cc
)

target_include_directories(Utils_ErrValue
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

install(TARGETS Utils Utils_ErrValue
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
