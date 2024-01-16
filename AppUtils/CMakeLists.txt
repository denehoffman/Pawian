file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(AppUtils
  ${sourcefiles}  
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
