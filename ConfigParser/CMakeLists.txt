file(GLOB sourcefiles *.cc)
foreach(item IN LISTS sourcefiles)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefiles ${item})
  endif()
endforeach()

add_library(ConfigParser  
  ${sourcefiles}  
)

target_include_directories(ConfigParser
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(ConfigParser  
  ErrLogger
  ${Boost_LIBRARIES}    
)

install(TARGETS ConfigParser
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
