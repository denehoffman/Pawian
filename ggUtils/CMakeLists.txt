add_library(ggUtils
  GGChannelEnv.cc
  centralProdBaseLh.cc
  gamgamBaseLh.cc
  ggReaction.cc
  ggStates.cc  
  pbarpProductionLh.cc 
)

target_include_directories(ggUtils
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}    
)

target_link_libraries(ggUtils
  ROOT::Minuit2
)

install(TARGETS ggUtils
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable(ggStatesApp
  ggStatesApp.cc
)

target_include_directories(ggStatesApp
  PRIVATE 
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(ggStatesApp
  ggUtils
  ConfigParser
  FitParams
  PwaUtils
  Utils
  Event
  pbarpUtils
  qft++
  ErrLogger
)

install(TARGETS ggStatesApp
  RUNTIME DESTINATION bin
)