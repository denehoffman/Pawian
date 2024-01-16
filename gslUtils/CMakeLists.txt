add_library(gslUtils  
  CMIntegration.cc
  CMIntegration2D.cc  
  CMIntegration2DReid.cc  
  CMIntegrationDudek.cc
  CMIntegrationReid.cc  
)

target_include_directories(gslUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/
    ${ROOT_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}  
)

target_link_libraries(gslUtils  
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger
  ROOT::Minuit2
  GSL::gsl
  GSL::gslcblas
)

install(TARGETS gslUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(CMIntegrationApp
  CMIntegrationApp.cc
)

target_include_directories(CMIntegrationApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(CMIntegrationApp
  gslUtils
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger  
)

add_executable(CMIntegration2DApp
  CMIntegration2DApp.cc
)

target_include_directories(CMIntegration2DApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(CMIntegration2DApp
  gslUtils
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger  
)

add_executable(convertToBinApp
  convertToBinApp.cc
)

target_include_directories(convertToBinApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(convertToBinApp
  gslUtils
  Utils
  Particle
  qft++
  qft++Extension
  ErrLogger  
)

install(TARGETS CMIntegrationApp CMIntegration2DApp convertToBinApp
  RUNTIME DESTINATION bin 
)