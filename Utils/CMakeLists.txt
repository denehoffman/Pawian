add_library(Utils  
  Faddeeva.cc
  FunctionUtils.cc
  IdStringMapRegistry.cc
  PawianIOUtils.cc
)

target_include_directories(Utils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

target_link_libraries(Utils  
  ErrLogger  
)

add_library(Utils_ErrValue
  ErrValue.cc
)

target_include_directories(Utils_ErrValue
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/  
)

install(TARGETS Utils Utils_ErrValue
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
