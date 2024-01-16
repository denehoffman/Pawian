file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(MinFunctions
  ${sourcefiles}
)

target_include_directories(MinFunctions
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(MinFunctions
  Utils
  qft++
  ErrLogger  
  ConfigParser
  FitParams
  PwaDynamics
  PspGen
  PwaUtils
  ${ROOT_LIBRARIES}/libMinuit2.so
  ${Boost_LIBRARIES}
)

install(TARGETS MinFunctions
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)