project(base)
    
  # Check for stacktrace and math library.
  
	# Additional building steps required by MSVC.
  include(BuildAndInstallScripts)

  #
  # Source files for the Hermes_common library.
  #
  set(SRC
    ord.cpp
    hermes_function.cpp
  )

  set(HEADERS
    ${CMAKE_HOME_DIRECTORY}/include/base/ord.h
    ${CMAKE_HOME_DIRECTORY}/include/base/hermes_function.h
  )
  
  #
  # Configure building and installation of all versions of hermes_common needed by
  # the selected Hermes configurations.
  #
  if(MSVC)
    BUILD_LIB(base ${DEBUG_FLAGS})
    target_link_libraries(base util)
    set_target_properties(base PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${FLAGS}")
    INSTALL_LIB(base)
  else(MSVC)
    if(HERMES_COMMON_DEBUG)
      set(BUILD_FLAGS "${DEBUG_FLAGS}")
    endif(HERMES_COMMON_DEBUG)
    BUILD_COMMON_LIB(base ${BUILD_FLAGS})
    INSTALL_LIB(base)
  endif(MSVC)

  # Install header files.

  file(GLOB INC_COMMON    "${CMAKE_HOME_DIRECTORY}/include/base/*.h")
  install(FILES ${INC_COMMON}    DESTINATION ${TARGET_ROOT}/include/base)