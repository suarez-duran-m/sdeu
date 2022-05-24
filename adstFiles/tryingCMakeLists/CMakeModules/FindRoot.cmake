################################################################################
# Module to find Root                                                          #
#                                                                              #
# This sets the following variables:                                           #
#   - Root_FOUND                                                               #
#   - ROOT_LIBRARIES                                                           #
#   - ROOT_GUI_LIBRARIES                                                       #
#   - ROOT_INCLUDE_DIR                                                         #
#   - ROOT_BIN_DIR                                                             #
#   - ROOTSYS                                                                  #
#   - ROOTCINT                                                                 #
#   - ROOT_CPPFLAGS                                                            #
#   - ROOT_LDFLAGS                                                             #
#   - ROOT_VERSION                                                             #
#                                                                              #
# And these might be needed as well... we'll see:                              #
#   - ROOT_LIBS                                                                #
#   - ROOT_GLIBS                                                               #
#   - ROOT_LIBDIR                                                              #
################################################################################

# If root was found already, don't find it again.
# This works because Root_FOUND is not cached
IF (NOT Root_FOUND)
  INCLUDE (CMakeCompareVersionStrings)

  # Minimum ROOT requirements.
  SET (ROOT_MIN_REQUIRED_VERSION "5.16.00")

  # the following disables all default paths (either from cmake, from environment)
  FIND_PATH (ROOT_BIN_DIR root-config
    PATHS $ENV{ROOTSYS}/bin
    $ENV{AUGER_BASE}/External/ROOT/pro/bin
    NO_DEFAULT_PATH
  )

  # now defaults are allowed but if nothing is found it is not overwritten
  # (because it is cached)
  FIND_PATH (ROOT_BIN_DIR root-config)

  GET_FILENAME_COMPONENT (ROOTSYS ${ROOT_BIN_DIR} PATH)

  IF (NOT ENV{ROOTSYS})
    SET (ENV{ROOTSYS} ${ROOTSYS})
  ENDIF ()

  IF (ROOTSYS)
    # ----------------------------------------------------------------------------
    # Get ROOT version, re-express in form XX.YY.ZZ, compare to requirements.
    # ----------------------------------------------------------------------------
    EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --version
      OUTPUT_VARIABLE ROOT_VERSION
    )
    STRING (REGEX REPLACE "[ \t\r\n]+" "" ROOT_VERSION "${ROOT_VERSION}")
    MESSAGE (STATUS "root-config --version ${ROOT_VERSION}")
    STRING (REGEX REPLACE "/" "." ROOT_VERSION "${ROOT_VERSION}")

    # Note: letters sometimes appended to ROOT patches (e.g., 5.14/00e) will be
    # ignored in the version comparison.
    COMPARE_VERSION_STRINGS (
      "${ROOT_VERSION}"
      ${ROOT_MIN_REQUIRED_VERSION}
      _have_MIN_ROOT_VER
    )

    IF (${_have_MIN_ROOT_VER} GREATER 0 OR ${_have_MIN_ROOT_VER} EQUAL 0)
      # --------------------------------------------------------------------------
      # Set ROOT compilation flags.
      # --------------------------------------------------------------------------
      SET (Root_FOUND TRUE)
      SET (HAVE_ROOT 1)

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --noldflags --libs
        OUTPUT_VARIABLE ROOT_LIBS
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" " " ROOT_LIBS "${ROOT_LIBS}")

      # Check if ROOT has Minuit2
      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --has-minuit2
        OUTPUT_VARIABLE HAVE_MINUIT2)

      IF (${HAVE_MINUIT2} STREQUAL "yes\n")
        SET (ROOT_LIBS "${ROOT_LIBS} -lMinuit2")
        SET (ROOT_MINUIT2_FOUND TRUE)
      ENDIF ()

      # Check if ROOT has MathMore
      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --has-mathmore
        OUTPUT_VARIABLE HAVE_MATHMORE)

      IF (${HAVE_MATHMORE} STREQUAL "yes\n")
        SET (ROOT_LIBS "${ROOT_LIBS} -lMathMore")
        SET (ROOT_MATHMORE_FOUND TRUE)
      ENDIF ()

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --libs
        OUTPUT_VARIABLE ROOT_LDFLAGS
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" " " ROOT_LDFLAGS "${ROOT_LDFLAGS}")

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --glibs
        OUTPUT_VARIABLE ROOT_GLIBS
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" "" ROOT_GLIBS "${ROOT_GLIBS}")

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --libdir
        OUTPUT_VARIABLE _ROOT_LIBDIR
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" "" ROOT_LIBDIR "${_ROOT_LIBDIR}")

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --incdir
        OUTPUT_VARIABLE ROOT_INCLUDE_DIR
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" "" ROOT_INCLUDE_DIR "${ROOT_INCLUDE_DIR}")

      EXECUTE_PROCESS (COMMAND ${ROOT_BIN_DIR}/root-config --cflags
        OUTPUT_VARIABLE ROOT_CPPFLAGS
      )
      STRING (REGEX REPLACE "[ \t\r\n]+" " " ROOT_CPPFLAGS "${ROOT_CPPFLAGS}")

      SET (ROOTCINT ${ROOT_BIN_DIR}/rootcint)

      # --------------------------------------------------------------------------
      # Find ROOT libraries.
      # --------------------------------------------------------------------------
      STRING (REGEX REPLACE "-l" " " CRAP "${ROOT_LIBS}")
      STRING (REGEX REPLACE "-" " " CRAP "${CRAP}")
      SEPARATE_ARGUMENTS (CRAP)
      FOREACH (_file ${CRAP})
        FIND_LIBRARY (_LIBY${_file} ${_file} ${ROOT_LIBDIR} NO_DEFAULT_PATH)
        IF (_LIBY${_file})
          LIST (APPEND ROOT_LIBRARIES ${_LIBY${_file}})
        ENDIF ()
        SET (_LIBY${_file} ${_LIBY${_file}} CACHE INTERNAL "" FORCE)
      ENDFOREACH ()

      STRING (REGEX REPLACE "-l" " " CRAP "${ROOT_GLIBS}")
      STRING (REGEX REPLACE "-" " " CRAP "${CRAP}")
      SEPARATE_ARGUMENTS (CRAP)
      FOREACH (_file ${CRAP})
        FIND_LIBRARY (_LIBY${_file} ${_file} ${ROOT_LIBDIR} NO_DEFAULT_PATH)
        IF (_LIBY${_file})
          LIST (APPEND ROOT_GUI_LIBRARIES ${_LIBY${_file}})
        ENDIF ()
        SET (_LIBY${_file} ${_LIBY${_file}} CACHE INTERNAL "" FORCE)
      ENDFOREACH ()

      FIND_LIBRARY (_ROOT_GQT_LIBRARY GQt ${ROOT_LIBDIR} NO_DEFAULT_PATH)
      FIND_LIBRARY (_ROOT_QTROOT_LIBRARY QtRoot ${ROOT_LIBDIR} NO_DEFAULT_PATH)
      IF (_ROOT_GQT_LIBRARY AND _ROOT_QTROOT_LIBRARY)
        SET (ROOT_QTROOT_LIBRARIES ${_ROOT_GQT_LIBRARY} ${_ROOT_QTROOT_LIBRARY} ${ROOT_GUI_LIBRARIES})
      ENDIF ()

    ELSE ()

      SET (Root_FOUND FALSE)
      SET (HAVE_ROOT 0)

      IF (NOT Root_FIND_QUIETLY)
        MESSAGE (STATUS
          "Compatible version of ROOT not found. "
          "Minimum required version: ${ROOT_MIN_REQUIRED_VERSION}")
      ENDIF ()

    ENDIF ()

  ENDIF ()

  SET (_DIRECTORIES)

  IF (NOT DL_LIBRARY)
    FIND_LIBRARY (DL_LIBRARY dl)
    MARK_AS_ADVANCED (DL_LIBRARY)
  ENDIF ()
  SET (ROOT_LIBRARIES ${ROOT_LIBRARIES} ${DL_LIBRARY})

  IF (Root_FOUND)
    #add minuit
    FIND_LIBRARY (ROOT_MINUIT_LIB Minuit ${ROOT_LIBDIR} NO_DEFAULT_PATH)

    IF (Root_FIND_REQUIRED_Minuit2 AND NOT ROOT_MINUIT2_FOUND)
      MESSAGE(FATAL_ERROR "
  ------------------------------------------------
    Could not find required Minuit2 library
  ------------------------------------------------\n")
    ENDIF ()

    IF (ROOT_MINUIT_LIB)
      SET (ROOT_LDFLAGS "${ROOT_LDFLAGS} -lMinuit")
      SET (ROOT_MINUIT_LIB ${ROOT_MINUIT_LIB} CACHE INTERNAL "" FORCE)
    ELSE ()
      MESSAGE(FATAL_ERROR "Could not find TMinuit library in ${ROOT_LIBDIR}")
    ENDIF ()
    LIST (APPEND ROOT_LIBRARIES ${ROOT_MINUIT_LIB})
    IF (NOT Root_FIND_QUIETLY)
      MESSAGE (STATUS "Root: ${ROOTSYS}")
    ENDIF ()
  ELSE ()
    IF (Root_FIND_REQUIRED)
      MESSAGE (FATAL_ERROR "Could not find Root")
    else ()
      if (NOT Root_FIND_QUIETLY)
        MESSAGE (STATUS "Could not find Root")
      endif ()
    ENDIF ()
  ENDIF ()


  ###########################################
  #
  #       Macros for building ROOT dictionary
  #
  ###########################################

  # The following macro creates rules for a dictionary
  # If the LinkDef file is called SomeLinkDef.h
  # it will create the files SomeDictionary.cc and SomeDictionary.h
  # it defines the variables $Some_DICT_DOURCE and $Some_DICT_HEADER
  # $Some_INFILES $Some_INCLUDE_DIRS $Some_VERBOSE

  # local variables: _LINKDEF _INCLUDE_DIRS _NAME _current_FILE

  MACRO (ROOT_GENERATE_DICTIONARY _NAME)
    PARSE_ARGUMENTS (${_NAME}
      "LINKDEF;INFILES;INCLUDE_DIRS;DEFINE_FLAGS"
      "VERBOSE"
      ${ARGN}
    )
    SET (_INCLUDE_DIRS)
    SET (${_NAME}_DICT_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${_NAME}Dictionary.cc)
    SET (${_NAME}_DICT_HEADER ${CMAKE_CURRENT_BINARY_DIR}/${_NAME}Dictionary.h)

    # Set up the list of includes
    IF (ROOT_VERSION VERSION_LESS "6")
      SET (_INCLUDE_DIR_PROP DIRECTORY "")
    ELSE()
      GET_PROPERTY (_INCLUDE_DIR_PROP DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
    ENDIF ()
    FOREACH (_current_FILE ${${_NAME}_INCLUDE_DIRS} ${_DIR_INCLUDES} ${_INCLUDE_DIR_PROP})
      SET (_INCLUDE_DIRS ${_INCLUDE_DIRS} -I${_current_FILE})
    ENDFOREACH ()

    # Set up the list of compiler definitions
    FOREACH (_current_DEF ${${_NAME}_DEFINE_FLAGS} ${_DEF_FLAGS})
      SET (_DEFINE_FLAGS ${_DEFINE_FLAGS} -D${_current_DEF})
    ENDFOREACH ()

    # Set up the call to rootcint
    IF (${_NAME}_VERBOSE)
      MESSAGE (STATUS "Root dictionary:\n   ${_NAME}_INFILES: ${${_NAME}_INFILES}")
      MESSAGE (STATUS "   OutFILEs: ${${_NAME}_DICT_SOURCE} ${${_NAME}_DICT_HEADER}")
      MESSAGE (STATUS "   LINKDEF: ${${_NAME}_LINKDEF}")
      MESSAGE (STATUS "   _DEFINE_FLAG: ${_DEFINE_FLAGS}")
    ENDIF ()

    ADD_CUSTOM_COMMAND (OUTPUT ${${_NAME}_DICT_SOURCE} ${${_NAME}_DICT_HEADER}
      COMMAND ${ROOTCINT} -f ${${_NAME}_DICT_SOURCE} -c ${_INCLUDE_DIRS} ${_DEFINE_FLAGS} ${${_NAME}_INFILES} ${${_NAME}_LINKDEF}
      DEPENDS ${${_NAME}_INFILES}
    )
  ENDMACRO ()

ENDIF ()
