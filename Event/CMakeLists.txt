CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

set(PACKAGE_NAME Event)

MESSAGE("cmake installation for package: " ${PACKAGE_NAME})

include_directories (
                    ${CMAKE_SOURCE_DIR}
                    ${log4cpp_INCLUDE_DIRS}
                    ${Boost_INCLUDE_DIRS}
)



AUX_SOURCE_DIRECTORY("." SRC)
Set(CCSRC "")

FOREACH (THEFILE ${SRC})
 IF (NOT ${THEFILE} MATCHES "(.*)App(.*)")
  SET (CCSRC ${CCSRC} ${THEFILE})
 ENDIF()
ENDFOREACH(THEFILE)

# BUILD_SHARED_LIBS controls the behavior of ADD_LIBRARY if STATIC/SHARED omitted
IF( PAWIAN_STATIC )
        ADD_LIBRARY ( ${PACKAGE_NAME} STATIC ${CCSRC} )
ELSE()
        ADD_LIBRARY ( ${PACKAGE_NAME} SHARED ${CCSRC} )
ENDIF( PAWIAN_STATIC )


INSTALL (TARGETS ${PACKAGE_NAME} DESTINATION ${CMAKE_BINARY_DIR}/lib)

if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET (EXECUTABLENAME EvtReaderTestApp)
ADD_EXECUTABLE( ${EXECUTABLENAME}
   EvtReaderTestApp.cc
)

target_link_libraries(${EXECUTABLENAME} Event)
target_link_libraries(${EXECUTABLENAME} Particle)
target_link_libraries(${EXECUTABLENAME} ErrLogger)
target_link_libraries(${EXECUTABLENAME} Utils)
target_link_libraries(${EXECUTABLENAME} qft++)
target_link_libraries(${EXECUTABLENAME} ${log4cpp_LIBRARIES})

INSTALL ( TARGETS ${EXECUTABLENAME} DESTINATION bin/${EXECUTABLENAME})

