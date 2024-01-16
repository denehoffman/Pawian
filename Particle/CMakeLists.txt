add_library(Particle  
  Parity.cc
  Particle.cc
  ParticleTable.cc
  PdtParser.cc  
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