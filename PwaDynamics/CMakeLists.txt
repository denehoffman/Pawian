file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(PwaDynamics  
  ${sourcefiles}
)

target_include_directories(PwaDynamics
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(PwaDynamics  
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger
  ROOT::Minuit2
)

install(TARGETS PwaDynamics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(BBUnstableParFitApp
  BBUnstableParFitApp.cc
)

target_include_directories(BBUnstableParFitApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(BBUnstableParFitApp
  PRIVATE
    PwaDynamics
    Utils
    Particle
    qft++
    qft++Extension
    ErrLogger
)

add_executable(LUTFileApp
  LUTFileApp.cc
)

target_include_directories(LUTFileApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(LUTFileApp
  PRIVATE
    PwaDynamics
    Utils
    Particle
    qft++
    qft++Extension
    ErrLogger
)

install(TARGETS BBUnstableParFitApp LUTFileApp
  RUNTIME DESTINATION bin  
)