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

add_executable(test_ReflectivityBasis
  test/TestReflectivityBasis.cc
)

target_include_directories(test_ReflectivityBasis
  PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
    ${Boost_UNIT_TEST_FRAMEWORK_INCLUDE_DIRS}
)

target_link_libraries(test_ReflectivityBasis
  PRIVATE
    gammapUtils
    ${Boost_LIBRARIES}
    ${ROOT_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(NAME test_ReflectivityBasis COMMAND test_ReflectivityBasis)
