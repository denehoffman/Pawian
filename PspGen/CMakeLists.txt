file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(PspGen
  ${sourcefiles} 
)

target_include_directories(PspGen
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(PspGen  
  ErrLogger
)

install(TARGETS PspGen
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(PspGenTestApp
  PspGenTestApp.cc
)

target_include_directories(PspGenTestApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}    
)

target_link_libraries(PspGenTestApp
  PspGen  
  ErrLogger
  ${ROOT_LIBRARIES}
)

install(TARGETS PspGenTestApp
  RUNTIME DESTINATION bin
)