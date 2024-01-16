file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(Event  
  ${sourcefiles}
)

target_include_directories(Event
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(Event  
  ErrLogger
  Particle
  ${Boost_LIBRARIES}    
  ${ROOT_LIBRARIES}      
)

add_executable(test_Event
  test/TestEventReaderDefault.cc
)

target_include_directories(test_Event
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
    ${Boost_UNIT_TEST_FRAMEWORK_INCLUDE_DIRS}
)

target_link_libraries(test_Event
  PRIVATE 
    ErrLogger
    Particle
    Event
    ${Boost_LIBRARIES}
    ${ROOT_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_test(NAME test_Event COMMAND test_Event)

install(TARGETS Event
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)


