add_library(PspGen
  EvtComplex.cc
  EvtConst.cc
  EvtGenKine.cc
  EvtKine.cc
  EvtMTRandomEngine.cc
  EvtRandom.cc
  EvtReport.cc
  EvtSimpleRandomEngine.cc
  EvtTensor4C.cc
  EvtVector3C.cc
  EvtVector3R.cc
  EvtVector4C.cc
  EvtVector4R.cc
  EvtdFunction.cc
  EvtdFunctionSingle.cc  
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