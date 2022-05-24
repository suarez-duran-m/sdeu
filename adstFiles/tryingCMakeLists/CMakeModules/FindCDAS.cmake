################################################################################
# Module to find CDAS packages                                                 #
# (This file looks for the new (post v5r0) CMAKE-ified cdas installation       #
#                                                                              #
# This module sets:                                                            #
#      CDAS_FOUND                                                              #
#      CDAS_LIBRARIES                                                          #
#      CDAS_INCLUDE_DIR                                                        #
#      CDAS_DEFINITIONS                                                        #
#      CDAS_V5R0                                                               #
#      CDAS_LDFLAGS                                                            #
#      CDAS_LIB_DIR                                                            #
#      IOAUGER_CPPFLAGS                                                        #
################################################################################
INCLUDE (AugerCheckCxxSourceRuns)
INCLUDE (CheckCXXSourceCompiles)

IF (NOT CDAS_FIND_QUIETLY)
  MESSAGE (STATUS "Searching for CDAS libraries equal to or newer than v5r0")
ENDIF ()

SET (DIRECTORIES
  $ENV{CDASHOME}
)

FOREACH (DIRECTORY ${DIRECTORIES})
  LIST (APPEND INC_DIRECTORIES ${DIRECTORY}/include)
  LIST (APPEND LIB_DIRECTORIES ${DIRECTORY}/lib)
ENDFOREACH ()

SET (LIBS Er;Ec;Es;IoAuger;IoSd;STCoordinates)
SET (CDAS_LIBRARIES)
FOREACH (LIB ${LIBS})
  FIND_LIBRARY (FOUND${LIB} ${LIB} ${LIB_DIRECTORIES})
  LIST (APPEND CDAS_LIBRARIES ${FOUND${LIB}})
ENDFOREACH ()

FIND_LIBRARY (FOUNDIoMd IoMd ${LIB_DIRECTORIES})
IF (FOUNDIoMd)
  MESSAGE (STATUS "Adding IoMd CDAS library")
  LIST (APPEND CDAS_LIBRARIES ${FOUNDIoMd})
ENDIF ()

FIND_LIBRARY (FOUNDUniv Univ ${LIB_DIRECTORIES})
IF (FOUNDUniv)
  MESSAGE (STATUS "Adding Universality CDAS library")
  LIST (APPEND CDAS_LIBRARIES ${FOUNDUniv})
ENDIF ()

IF (INC_DIRECTORIES AND CDAS_LIBRARIES)

  SET (CDAS_INCLUDE_DIR ${INC_DIRECTORIES})
  SET (CDAS_FOUND TRUE)
  SET (HAVE_IOSD TRUE)
##################################################
#### Check if CDAS was compiled with -DIOMD=ON or with -DIORD=ON ###

  LIST (APPEND AERAROOTIO_INCLUDE_DIR $ENV{AERAROOTIO}/inc)
  FIND_LIBRARY (AERAROOTIO_LIBRARIES aerarootio $ENV{AERAROOTIO}/lib NO_DEFAULT_PATH)
  FIND_LIBRARY (AERAROOTIO_LIBRARIES aerarootio $ENV{AERAROOTIO}/lib)

  # do not force to have AERA libraries
  IF (AERAROOTIO_LIBRARIES)
    MESSAGE (STATUS "AERA libs found")
    SET (CMAKE_REQUIRED_INCLUDES ${CDAS_INCLUDE_DIR} ${FDEVENTLIB_INCLUDE_DIR} ${ROOT_INCLUDE_DIR} ${AERAROOTIO_INCLUDE_DIR})
    SET (CMAKE_REQUIRED_LIBRARIES ${CDAS_LIBRARIES} ${FDEVENTLIB_LIBRARIES} ${ROOT_LIBRARIES} ${AERAROOTIO_LIBRARIES})
  ELSE ()
    MESSAGE (STATUS "AERA libs not found")
    SET (CMAKE_REQUIRED_INCLUDES ${CDAS_INCLUDE_DIR} ${FDEVENTLIB_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})
    SET (CMAKE_REQUIRED_LIBRARIES ${CDAS_LIBRARIES} ${FDEVENTLIB_LIBRARIES} ${ROOT_LIBRARIES})
  ENDIF ()

  SET (CMAKE_REQUIRED_DEFINITIONS "-DHAVE_IOMD -DHAVE_IORD") # do not remove!
  CHECK_CXX_SOURCE_COMPILES (
    "
    #include <IoAuger.h>
    #include <MdEvent.h>

    int
    main()
    {
      AugerEvent ev;
      md::Event& md = ev.Md();
      return !&md;
    }
    "
    HAVE_IOMD
  )

  IF (HAVE_IOMD)
    SET (CMAKE_REQUIRED_DEFINITIONS "-DHAVE_IOMD -DHAVE_IORD") # do not remove!
  
    CHECK_CXX_SOURCE_COMPILES (
      "
      #include <MdEvent.h>

      int
      main()
      {
        md::Module mdMod;
        return !&mdMod.GetIntegrator();
        
      }
      "
      IOMD_V1R7
    )
     
    CHECK_CXX_SOURCE_COMPILES (
      "
      #include <MdEvent.h>

      int
      main()
      {
      
        md::Module mdMod;
        return !&mdMod.GetIntegratorA();
      }
      "
      IOMD_V1R8
    )
  ENDIF ()

  SET (CMAKE_REQUIRED_DEFINITIONS "-DHAVE_IOMD -DHAVE_IORD") # do not remove!
  CHECK_CXX_SOURCE_COMPILES (
    "
    #include <IoAuger.h>
    #include <AERAevent.h>

    int
    main()
    {
      AugerEvent ev;
      AERAevent& ae = ev.GiveREvent();
      return !&ae;
    }
    "
    HAVE_IORD
  )

##################################################
  SET (CDAS_V5R0 TRUE)

  IF (HAVE_IOMD)
    LIST (APPEND CDAS_DEFINITIONS -DHAVE_IOMD)
    IF (IOMD_V1R8)
    	LIST (APPEND CDAS_DEFINITIONS -DIOMD_V1R8)
  	ELSEIF (IOMD_V1R7)
    	LIST (APPEND CDAS_DEFINITIONS -DIOMD_V1R7)
  	ENDIF ()
  ENDIF ()
  IF (HAVE_IORD)
    LIST (APPEND CDAS_DEFINITIONS -DHAVE_IORD)
  ENDIF ()
  STRING (REPLACE ";" " " IOAUGER_CPPFLAGS "${CDAS_DEFINITIONS}")
  SET (IOAUGER_CPPFLAGS "-I${CDAS_INCLUDE_DIR} ${IOAUGER_CPPFLAGS}")

  IF (NOT CDAS_FIND_QUIETLY)
    MESSAGE (STATUS "CDAS definitions : ${CDAS_DEFINITIONS}")
  ENDIF ()

ENDIF ()

IF (CDAS_LIBRARIES)
  MAKE_LDFLAGS (CDAS_LDFLAGS "${CDAS_LIBRARIES}")
  SET (CDAS_LIB_DIR ${LIB_DIRECTORIES})
ENDIF ()
