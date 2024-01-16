add_library(AppUtils
  AppBase.cc
  LhFactory.cc  
)

target_include_directories(AppUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}    
)

target_link_libraries(AppUtils
  Utils
  qft++
  ErrLogger
  ConfigParser
  FitParams
  Event
  PwaDynamics
  PspGen
  PwaUtils
  MinFunctions
  pbarpUtils
  epemUtils
  resUtils
  ggUtils
  pipiScatteringUtils
  qaErrorExtract
)

install(TARGETS AppUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
