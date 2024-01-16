file(GLOB sourcefilesRQM relativistic-quantum-mechanics/*.cc)
file(GLOB sourcefilesTensor tensor/*.cc)
file(GLOB sourcefilesMatrix matrix/*.cc)

foreach(item IN LISTS sourcefilesRQM)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefilesRQM ${item})
  endif()
endforeach()

foreach(item IN LISTS sourcefilesTensor)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefilesTensor ${item})
  endif()
endforeach()

foreach(item IN LISTS sourcefilesMatrix)  
  if(${item} MATCHES "App.cc")
    list(REMOVE_ITEM sourcefilesMatrix ${item})
  endif()
endforeach()

add_library(qft++
  ${sourcefilesRQM}
  ${sourcefilesTensor} 
  ${sourcefilesMatrix} 
)

target_include_directories(qft++
  PRIVATE 
    ${CMAKE_SOURCE_DIR}/
    topincludes/  
)

install(TARGETS qft++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
