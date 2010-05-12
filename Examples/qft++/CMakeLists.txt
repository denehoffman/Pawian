CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

MESSAGE("cmake installation for package: Examples/qft++")

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

######## Linking "EtacToa1320pi0fit" ###########################

SET (EXECUTABLENAME EtacToa1320pi0fitApp)

SET (${EXECUTABLENAME}_SRCS
        EtacToa1320pi0Fcn.cc 
	EtacToa1320pi0fit.cc
	EtacToa1320pi0fitMain.cc
)

ADD_EXECUTABLE(${EXECUTABLENAME}
	${${EXECUTABLENAME}_SRCS}
)
MESSAGE("ROOT_LIBRARIES: " ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES (${EXECUTABLENAME}
        ${Minuit2_LIBRARIES}
        ${ROOT_LIBRARIES}
        qft++
)

INSTALL ( TARGETS ${EXECUTABLENAME} DESTINATION bin/${EXECUTABLENAME})

######## End of Linking "EtacToa1320pi0fit" ###########################
