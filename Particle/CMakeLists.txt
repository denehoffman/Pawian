file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(Particle  
  ${sourcefiles}
)

target_include_directories(Particle
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(Particle  
  ErrLogger
  Utils_ErrValue
)

install(TARGETS Particle
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(PdtTestApp
  PdtTestApp.cc
)

target_include_directories(PdtTestApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(PdtTestApp
  Particle
  Utils_ErrValue
  ErrLogger  
)

install(TARGETS PdtTestApp
  RUNTIME DESTINATION bin  
)