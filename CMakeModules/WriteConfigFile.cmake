MACRO (WRITE_CONFIG_FILE filename)

  MESSAGE(STATUS "filename ${filename}")

  IF(${filename} MATCHES "[.]csh")
    SET(CSH_OUT TRUE)
  ELSE(${filename} MATCHES "[.]csh")
    SET(CSH_OUT FALSE)
  ENDIF(${filename} MATCHES "[.]csh")
  
 SET(PAWIANLIBDIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
  
  IF(CMAKE_SYSTEM_NAME MATCHES Linux)
    SET(LD_LIBRARY_PATH ${PAWIANLIBDIR} ${LD_LIBRARY_PATH} /usr/lib /usr/X11R6/lib ) 
  ELSE(CMAKE_SYSTEM_NAME MATCHES Linux)
    IF(CMAKE_SYSTEM_NAME MATCHES Darwin)
      SET(LD_LIBRARY_PATH ${FAIRLIBDIR} ${LD_LIBRARY_PATH})
    ENDIF(CMAKE_SYSTEM_NAME MATCHES Darwin)
  ENDIF(CMAKE_SYSTEM_NAME MATCHES Linux)
  

# configure_file(${PROJECT_SOURCE_DIR}/CMakeModules/scripts/check_system.sh.in
#                   ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh
#                  )

#  IF(CMAKE_SYSTEM_NAME MATCHES Linux)
#    FILE(READ /etc/issue _linux_flavour)
#    STRING(REGEX REPLACE "[\\]" " " _result1 "${_linux_flavour}")
#    STRING(REGEX REPLACE "\n" ";" _result "${_result1}")
#    SET(_counter 0)
#    FOREACH(_line ${_result})
#      if (_counter EQUAL 0)
#        SET(_counter 1)
#        set(_linux_flavour ${_line})
#      endif (_counter EQUAL 0)
#    ENDFOREACH(_line ${_result})
#    EXECUTE_PROCESS(COMMAND uname -m 
#                    OUTPUT_VARIABLE _system 
#                    OUTPUT_STRIP_TRAILING_WHITESPACE
#                   )
#   
#    WRITE_TO_FILE(${filename} Linux_Flavour_ ${_linux_flavour} "")
#    WRITE_TO_FILE(${filename} System_ ${_system} APPEND)
#    WRITE_FILE(${filename} ". ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh" APPEND)
#    WRITE_FILE(${filename} " if [ \"$same_system\" == \"1\" ]; then" APPEND)
#  ENDIF(CMAKE_SYSTEM_NAME MATCHES Linux)



  MESSAGE("config filename: " ${filename})

  WRITE_TO_FILE(${filename} CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR} APPEND)

  WRITE_TO_FILE(${filename} ROOTSYS ${ROOTSYS} APPEND) 
 
  CONVERT_LIST_TO_STRING(${LD_LIBRARY_PATH})

  IF(CMAKE_SYSTEM_NAME MATCHES Linux)
    WRITE_TO_FILE(${filename} LD_LIBRARY_PATH ${output} APPEND)
  ELSE(CMAKE_SYSTEM_NAME MATCHES Linux)
    IF(CMAKE_SYSTEM_NAME MATCHES Darwin)
      WRITE_TO_FILE(${filename} DYLD_LIBRARY_PATH ${output} APPEND)
    ENDIF(CMAKE_SYSTEM_NAME MATCHES Darwin)
  ENDIF(CMAKE_SYSTEM_NAME MATCHES Linux)



  IF($ENV{PATH} MATCHES ^${ROOT_BINARY_DIR}) 
    SET (PATH "$ENV{PATH}")
  ELSE()
   SET (PATH "${ROOT_BINARY_DIR}:$ENV{PATH}")      
  ENDIF($ENV{PATH} MATCHES ^${ROOT_BINARY_DIR})

  CONVERT_LIST_TO_STRING(${PATH})
  WRITE_TO_FILE(${filename} PATH ${output} APPEND)
ENDMACRO (WRITE_CONFIG_FILE)


MACRO (CONVERT_LIST_TO_STRING)

  set (tmp "")
  foreach (_current ${ARGN})

    set(tmp1 ${tmp})
    set(tmp "")
    set(tmp ${tmp1}:${_current})

  endforeach (_current ${ARGN})
  STRING(REGEX REPLACE "^:(.*)" "\\1" output ${tmp}) 

ENDMACRO (CONVERT_LIST_TO_STRING LIST)

MACRO (WRITE_TO_FILE FILENAME ENVVARIABLE VALUE OPTION)

  IF(CSH_OUT)
    WRITE_FILE(${FILENAME} "setenv ${ENVVARIABLE} \"${VALUE}\"" ${OPTION})
  ELSE(CSH_OUT)
    WRITE_FILE(${FILENAME} "export ${ENVVARIABLE}=\"${VALUE}\"" ${OPTION})
  ENDIF(CSH_OUT)

ENDMACRO (WRITE_TO_FILE)

MACRO(UNIQUE var_name list)

  #######################################################################
  # Make the given list have only one instance of each unique element and
  # store it in var_name.
  #######################################################################

  SET(unique_tmp "")
  FOREACH(l ${list})
    STRING(REGEX REPLACE "[+]" "\\\\+" l1 ${l})
    IF(NOT "${unique_tmp}" MATCHES "(^|;)${l1}(;|$)")
      SET(unique_tmp ${unique_tmp} ${l})
    ENDIF(NOT "${unique_tmp}" MATCHES "(^|;)${l1}(;|$)")
  ENDFOREACH(l)
  SET(${var_name} ${unique_tmp})
ENDMACRO(UNIQUE)







