file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(ErrLogger
  ${sourcefiles}
)

install(TARGETS ErrLogger
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(ErrLoggerTestApp
  ErrLoggerTestApp.cc
)

target_link_libraries(ErrLoggerTestApp
  ErrLogger
)

add_executable(testLogger
  TestLogger.cpp
)

target_link_libraries(testLogger
  ErrLogger
)

install(TARGETS ErrLoggerTestApp testLogger
  RUNTIME DESTINATION bin
)
