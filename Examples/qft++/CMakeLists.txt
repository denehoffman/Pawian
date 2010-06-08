CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3 FATAL_ERROR)

MESSAGE("cmake installation for package: Examples/qft++")

SET (INCLUDE_DIRECTORIES
     ${ROOT_INCLUDE_DIR}
     ${Minuit2_INCLUDE_DIRS}
     ${CMAKE_SOURCE_DIR} )

INCLUDE_DIRECTORIES (AFTER
 ${ROOT_INCLUDE_DIR}
 ${Minuit2_INCLUDE_DIRS}
 ${CMAKE_SOURCE_DIR}
)

LINK_DIRECTORIES(
   ${Minuit2_LIBRARY_DIRS}
   ${ROOT_LIBRARY_DIR})

######## Linking "compton_scatteringApp ###########################

SET (EXECUTABLENAME compton_scatteringApp)

SET (${EXECUTABLENAME}_SRCS 
	compton_scattering.cc
)

ADD_EXECUTABLE(${EXECUTABLENAME}
	${${EXECUTABLENAME}_SRCS}
)

TARGET_LINK_LIBRARIES (${EXECUTABLENAME}
        qft++
)

INSTALL ( TARGETS ${EXECUTABLENAME} DESTINATION bin/${EXECUTABLENAME})

######## End: Linking "compton_scatteringApp ###########################

######## Linking "SeqDecayHeliAmpsApp ###########################
SET (EXECUTABLENAME SeqDecayHeliAmpsApp)


SET (${EXECUTABLENAME}_SRCS 
	SeqDecHeliAmps.cc
	SeqDecayHeliAmpsMain.cc
)

ADD_EXECUTABLE(${EXECUTABLENAME}
	${${EXECUTABLENAME}_SRCS}
)

TARGET_LINK_LIBRARIES (${EXECUTABLENAME}
        qft++
	${ROOT_LIBRARIES}
)

INSTALL ( TARGETS ${EXECUTABLENAME1} DESTINATION bin/${EXECUTABLENAME})

######## End of Linking "SeqDecayHeliAmpsApp ###########################


