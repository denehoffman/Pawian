add_library(ErrLogger
  ErrLogger.cc
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
