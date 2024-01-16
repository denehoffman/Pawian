file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(FitParams
  ${sourcefiles}
)

target_include_directories(FitParams
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(FitParams  
  ROOT::Minuit2
)

install(TARGETS FitParams
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(StreamerTestApp
  StreamerTestApp.cc
)

target_include_directories(StreamerTestApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(StreamerTestApp
  FitParams
  qft++
  ConfigParser
  qft++Extension
  ErrLogger
  ROOT::Minuit2
)

add_executable(PwaCovMatrixApp
  PwaCovMatrixApp.cc
)

target_include_directories(PwaCovMatrixApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(PwaCovMatrixApp
  FitParams
  qft++
  ConfigParser
  qft++Extension
  ErrLogger 
)

add_executable(ChangeCovMatrixApp
  ChangeCovMatrixApp.cc
)

target_include_directories(ChangeCovMatrixApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}      
)

target_link_libraries(ChangeCovMatrixApp
  FitParams
  qft++
  ConfigParser
  qft++Extension
  ErrLogger  
)

add_executable(CreateCovMatrixApp
  CreateCovMatrixApp.cc
)

target_include_directories(CreateCovMatrixApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(CreateCovMatrixApp
  FitParams
  qft++
  ConfigParser
  qft++Extension
  ErrLogger  
)

install(TARGETS StreamerTestApp PwaCovMatrixApp ChangeCovMatrixApp CreateCovMatrixApp
  RUNTIME DESTINATION bin
)
