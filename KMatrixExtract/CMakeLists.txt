file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(KMatrixExtract  
  ${sourcefiles}
)

target_include_directories(KMatrixExtract
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}   
)

target_link_libraries(KMatrixExtract    
  ${ROOT_LIBRARIES}/libMinuit2.so 
)

install(TARGETS KMatrixExtract
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(TMatrixErrorExtrApp
  TMatrixErrorExtrApp.cc
)

target_include_directories(TMatrixErrorExtrApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(TMatrixErrorExtrApp
  KMatrixExtract
  pipiScatteringUtils
  PwaDynamics
  ConfigParser
  PwaUtils
  Event
  FitParams
  Particle
  qft++Extension
  qft++
  ErrLogger
)

add_executable(TMatrixResidueExtrApp
  TMatrixResidueExtrApp.cc
)

target_include_directories(TMatrixResidueExtrApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(TMatrixResidueExtrApp
  KMatrixExtract
  pipiScatteringUtils
  PwaDynamics
  ConfigParser
  PwaUtils
  Event
  FitParams
  Particle
  qft++Extension
  qft++
  ErrLogger
)

add_executable(ArgandModuloConverterApp
  ArgandModuloConverterApp.cc
)

target_include_directories(ArgandModuloConverterApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(ArgandModuloConverterApp
  KMatrixExtract
  pipiScatteringUtils
  PwaDynamics
  ConfigParser
  PwaUtils
  Event
  FitParams
  Particle
  qft++Extension
  qft++
  ErrLogger
)

add_executable(TMatrixGeneralApp
  TMatrixGeneralApp.cc
)

target_include_directories(TMatrixGeneralApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(TMatrixGeneralApp
  KMatrixExtract
  pipiScatteringUtils
  PwaDynamics
  ConfigParser
  PwaUtils
  Event
  FitParams
  Particle
  qft++Extension
  qft++
  ErrLogger
  Utils
)

add_executable(TMatrixPiPiSWaveSimple4piPhpApp
  TMatrixPiPiSWaveSimple4piPhpApp.cc
)

target_include_directories(TMatrixPiPiSWaveSimple4piPhpApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(TMatrixPiPiSWaveSimple4piPhpApp
  KMatrixExtract
  pipiScatteringUtils
  PwaDynamics
  ConfigParser
  PwaUtils
  Event
  FitParams
  Particle
  qft++Extension
  qft++
  ErrLogger
  Utils
)

install(TARGETS TMatrixErrorExtrApp TMatrixResidueExtrApp ArgandModuloConverterApp TMatrixGeneralApp TMatrixPiPiSWaveSimple4piPhpApp
  RUNTIME DESTINATION bin  
)