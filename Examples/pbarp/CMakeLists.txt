CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

MESSAGE("cmake installation for package: Examples/pbarp")

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET (INCLUDE_DIRECTORIES
     ${ROOT_INCLUDE_DIR}
     ${Minuit2_INCLUDE_DIRS}
     ${CMAKE_SOURCE_DIR} )

INCLUDE_DIRECTORIES (AFTER
 ${ROOT_INCLUDE_DIR}
 ${Minuit2_INCLUDE_DIRS}
 ${CMAKE_SOURCE_DIR}
)


######## Linking pbarpStatesApp ###########################
SET (EXECUTABLENAME pbarpStatesApp)


SET (${EXECUTABLENAME}_SRCS 
	pbarpStatesApp.cc
)

ADD_EXECUTABLE(${EXECUTABLENAME}
	${${EXECUTABLENAME}_SRCS}
)

TARGET_LINK_LIBRARIES (${EXECUTABLENAME} PwaUtils)
TARGET_LINK_LIBRARIES (${EXECUTABLENAME} qft++)
TARGET_LINK_LIBRARIES (${EXECUTABLENAME} ErrLogger)


INSTALL ( TARGETS ${EXECUTABLENAME} DESTINATION bin/${EXECUTABLENAME})

######## End of Linking pbarpStatesApp ###########################

