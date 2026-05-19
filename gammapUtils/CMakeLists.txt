file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(gammapUtils
  ${sourcefiles}
)

target_include_directories(gammapUtils
  PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(gammapUtils
  Utils
  qft++
  ErrLogger
  ConfigParser
  FitParams
  PwaUtils
  ROOT::Hist
)

install(TARGETS gammapUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(gammapStatesApp
  gammapStatesApp.cc
)

target_include_directories(gammapStatesApp
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(gammapStatesApp
  gammapUtils
)

install(TARGETS gammapStatesApp
  RUNTIME DESTINATION bin
)
