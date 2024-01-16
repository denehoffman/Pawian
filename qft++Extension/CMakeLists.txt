add_library(qft++Extension
  AbsComplexFunction.cc
  InputOmnesMatrixPiPiS5Channel.cc
  InputSelfEnergyMatrixPiPiS5Channel.cc
  OmnesMatrixFactory.cc
  OmnesMatrixGeneral.cc
  OmnesMatrixPiPiS3Channelwrhorho.cc
  OmnesMatrixPiPiS5Channel.cc
  PawianUtils.cc
  SphericalHarmonic.cc
)

target_include_directories(qft++Extension
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS} 
)

target_link_libraries(qft++Extension
  qft++
  ErrLogger
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES}
)

add_executable(OmnesMatrixGeneralApp
  OmnesMatrixGeneralApp.cc
)

target_include_directories(OmnesMatrixGeneralApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}  
)

target_link_libraries(OmnesMatrixGeneralApp
  ErrLogger
  qft++
  qft++Extension
  ConfigParser
  Particle
  ${ROOT_LIBRARIES}
)

install(TARGETS qft++Extension
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(TARGETS OmnesMatrixGeneralApp
  RUNTIME DESTINATION lib  
)