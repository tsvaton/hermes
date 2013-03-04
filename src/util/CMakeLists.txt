project(util)
    
  # Check for stacktrace and math library.
  include(CheckFunctionExists)
  CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
  CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)

	# Additional building steps required by MSVC.
  include(BuildAndInstallScripts)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/solvers)
  
  #
  # Source files for the Hermes_common library.
  #
  set(SRC
    mixins.cpp
    callstack.cpp
    api.cpp
    qsort.cpp
    c99_functions.cpp
    exceptions.cpp
  )

  set(HEADERS
    ${CMAKE_HOME_DIRECTORY}/include/util/mixins.h
    ${CMAKE_HOME_DIRECTORY}/include/util/callstack.h
    ${CMAKE_HOME_DIRECTORY}/include/util/common.h
    ${CMAKE_HOME_DIRECTORY}/include/util/api.h
    ${CMAKE_HOME_DIRECTORY}/include/util/array.h
    ${CMAKE_HOME_DIRECTORY}/include/util/qsort.h
    ${CMAKE_HOME_DIRECTORY}/include/util/c99_functions.h
    ${CMAKE_HOME_DIRECTORY}/include/util/exceptions.h
    ${CMAKE_HOME_DIRECTORY}/include/util/vector.h
  )
  
  #
  # Configure building and installation of all versions of hermes_common needed by
  # the selected Hermes configurations.
  #
  if(MSVC)
    BUILD_LIB(util ${DEBUG_FLAGS})
    set_target_properties(util PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${FLAGS}")
    INSTALL_LIB(util)
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${DEBUG_FLAGS}")
    endif(DEBUG)
    BUILD_COMMON_LIB(util ${BUILD_FLAGS})
    INSTALL_LIB(util)
  endif(MSVC)

  # Install header files.

  file(GLOB INC_COMMON    "${CMAKE_HOME_DIRECTORY}/include/util/*.h")
  install(FILES ${INC_COMMON}    DESTINATION ${TARGET_ROOT}/include/util)