set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_HOME_DIRECTORY}/cmake/cxx_flag_overrides.cmake)
project(hermes)
  cmake_minimum_required(VERSION 2.8)

  #------------------------------------------------------------------------------
  # Main configuration options.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ### BEGIN DEFAULT CONFIGURATION ###

  # The values below are sane defaults to produce a build on all platforms,
  # if you want to override them, set them in CMake.vars

  # Global
    # Optional stacktrace
    set(WITH_STACKTRACE NO)
  
  # Debug or release.
    set(DEBUG     YES)

  # Hermes2D:
  set(WITH_H2D                  YES)
    # Optional parts of the library.
    set(H2D_WITH_GLUT           YES)
    set(H2D_WITH_TEST_EXAMPLES  YES)
	
	# Advanced settings.
	# Number of solution / filter components.
	# Minimum is 1, maximum 999.
	set(H2D_MAX_COMPONENTS 10)
	# Number of quadratures, recommended is 4.
	# Minimum is 2, maximum 999.
	set(H2D_MAX_QUADRATURES 4)
	#  This number has to be smaller than 21.
	set(H2D_MAX_HANGING_NODES_LEVEL 15)
	#  Use second derivatives (that are turned off by default).
	set(H2D_USE_SECOND_DERIVATIVES YES)

  # Installation path:
    set(TARGET_ROOT             /usr/local) #default installation path for the library will be /usr/local/lib,
                                            #for include files /usr/local/include
  #### Solvers ###
  # Standard UMFPACK
    # UMFPACK is chosen by default.
    set(WITH_UMFPACK YES)
    
  ################################## CURRENTLY NOT SUPPORTED ###################################
  #  \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/   \/
  
  # SuperLU
    # First option enables sequential SuperLU, which should reside in SUPERLU_ROOT/lib
    # and SUPERLU_ROOT/include. If SUPERLU_ROOT is not provided, the library will be searched
    # for in system directories. If SUPERLU_MT is set to YES, the multithreaded version of
    # SuperLU will be used. In this case, there should be either superlu_mt_PTHREAD.a and/or
    # superlu_mt_OPENMP.a libraries present in the search path (the former will be chosen if
    # WITH_OPENMP is set to NO, the latter otherwise). There should also be a subdirectory
    # SLU_MT in the include directory since SuperLU_MT uses different versions of the same-named
    # header files as the serial library.
    set(WITH_SUPERLU            NO)
    set(SUPERLU_ROOT            /opt/SuperLU_4.0/)
    set(SUPERLU_MT              NO)

  # PETSc
    # If PETSc is enabled, PETSC_ROOT and PETSC_ARCH must be set as well;
    # if PETSc is disabled, PETSC_ROOT and PETSC_ARCH will be ignored.
    set(WITH_PETSC              NO)
    set(PETSC_ROOT              /opt/petsc/petsc-3.1-p7)
    set(PETSC_ARCH              petsc-arch-without-scalar-type)

  # MUMPS
    # Sequential version of MUMPS may be enabled by setting YES on the following line
    # and pointing to the package installation directory with MUMPS_ROOT. To enable the
    # parallel version, setting WITH_MPI to YES in addition to that should suffice. If link
    # problems arise, it may be that the automatically found SCALAPACK libraries (reported
    # during the CMake configuration phase) do not match those used when compiling MUMPS.
    # In this case, set SCALAPACK_LIBRARIES manually in CMake.vars and reconfigure Hermes.
    set(WITH_MUMPS              NO)
    # set(MUMPS_ROOT            root/dir/of/mumps/)
    # set(SCALAPACK_LIBRARIES   -lscalapack -lblacs)

  # Trilinos
    # Enable support for Trilinos solvers.
    set(WITH_TRILINOS           NO)
    
    # Experimental
    set(WITH_ZOLTAN             NO)

  ### Mesh format ###
    set(WITH_EXODUSII           NO)
    set(WITH_HDF5               NO)

  ### Others ###
  # Parallel execution.
    # (tells the linker to use parallel versions of the selected solvers, if available):
    #
    # Enable OpenMP for shared memory processing.
    #
    set(WITH_OPENMP             YES)

    # If MPI is enabled, the MPI library installed on the system should be found by
    # CMake automatically. If the found library doesn't match the one used to compile the
    # particular MPI-dependent package, the other two options should be used to specify it.
    #
    set(WITH_MPI                NO)
    # set(MPI_LIBRARIES         -lmpi)
    # set(MPI_INCLUDE_PATH      /usr/include/openmpi

  #  /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\   /\
	################################## CURRENTLY NOT SUPPORTED ###################################
  
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)
  
	  ### END DEFAULT CONFIGURATION ###
		
	  set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
	  set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR})
	  
	  if(${H2D_MAX_HANGING_NODES_LEVEL} MATCHES "[1-20]")
	  else(${H2D_MAX_HANGING_NODES_LEVEL} MATCHES "[1-20]")
		message(FATAL_ERROR "H2D_MAX_HANGING_NODES_LEVEL has to be a number smaller than 21.")
	  endif(${H2D_MAX_HANGING_NODES_LEVEL} MATCHES "[1-20]")
	  
	  if(${H2D_MAX_COMPONENTS} MATCHES "[1-999]")
	  else(${H2D_MAX_COMPONENTS} MATCHES "[1-999]")
		message(FATAL_ERROR "H2D_MAX_COMPONENTS has to be a number between 1 and 999.")
	  endif(${H2D_MAX_COMPONENTS} MATCHES "[1-999]")
	  
	  if(${H2D_MAX_QUADRATURES} MATCHES "[1-999]")
	  else(${H2D_MAX_QUADRATURES} MATCHES "[1-999]")
		message(FATAL_ERROR "H2D_MAX_QUADRATURES has to be a number between 1 and 999.")
	  endif(${H2D_MAX_QUADRATURES} MATCHES "[1-999]")
	  
	  #ADDED BY COLMAN
	  IF (APPLE)
    IF (NOT CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
      IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        ## It's Mac OS X.

        #This is causing problems when building the debug version.
        #SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Release is default." FORCE)

        #SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING
        #  "compiling for OSX only in i386 format. Library with the
        #   necessary symbol '_glutSetOption' was only available
        #   in /usr/lib/libglut.dylib, which happened to be a pure
        #   i386 library. Maybe fixed later." FORCE)
        
		IF (CMAKE_VERBOSE_MAKEFILE)
          MESSAGE(STATUS)
          MESSAGE(STATUS "*******************************************************************************")
          MESSAGE(STATUS "*        SUMMARY OF USED VARIABLES -> Details in AllVariables.txt             *")
          MESSAGE(STATUS "*******************************************************************************")
          MESSAGE(STATUS "COMPILING FOR OSX, platforms: ${CMAKE_OSX_ARCHITECTURES}")
          MESSAGE(STATUS "CMAKE_ROOT ........................... = ${CMAKE_ROOT}")
          MESSAGE(STATUS "CMAKE_HOME_DIRECTORY ................. = ${CMAKE_HOME_DIRECTORY}")
          MESSAGE(STATUS "CMAKE_VERSION ........................ = ${CMAKE_VERSION}")
          MESSAGE(STATUS "CMAKE_GENERATOR ...................... = ${CMAKE_GENERATOR}")
          MESSAGE(STATUS "CMAKE_MODULE_PATH .................... = ${CMAKE_MODULE_PATH}")
          MESSAGE(STATUS "CMAKE_HOST_SYSTEM .................... = ${CMAKE_HOST_SYSTEM}")
          MESSAGE(STATUS "CMAKE_HOST_SYSTEM_VERSION ............ = ${CMAKE_HOST_SYSTEM_VERSION}")
          MESSAGE(STATUS "CMAKE_HOST_SYSTEM_NAME ............... = ${CMAKE_HOST_SYSTEM_NAME}")
          MESSAGE(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR .......... = ${CMAKE_HOST_SYSTEM_PROCESSOR}")
          MESSAGE(STATUS "CMAKE_HOST_UNIX ...................... = ${CMAKE_HOST_UNIX}")
          MESSAGE(STATUS "CMAKE_BUILD_TYPE ..................... = ${CMAKE_BUILD_TYPE}")
          MESSAGE(STATUS "CMAKE_C_COMPILER ..................... = ${CMAKE_C_COMPILER}")
          MESSAGE(STATUS "CMAKE_CXX_COMPILER ................... = ${CMAKE_CXX_COMPILER}")
        ENDIF()
      ENDIF ()
    ENDIF (NOT CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
    # Use modified search procedures for some libraries on MacOS.
    include(${CMAKE_HOME_DIRECTORY}/OSX.cmake)
  ENDIF (APPLE)
  # end ADDED BY COLMAN

  # Determine the number of cores.
  if(WIN32)
    set(NUM_THREADS "$ENV{NUMBER_OF_PROCESSORS}")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    execute_process(COMMAND sysctl -a COMMAND grep "machdep.cpu.thread_count" COMMAND grep -o "[0-9]" OUTPUT_VARIABLE NUM_THREADS)
  else()
    set(cpuinfo_file "/proc/cpuinfo")
    if(EXISTS "${cpuinfo_file}")
      file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
      list(LENGTH procs NUM_THREADS)
    else()
      set(NUM_THREADS 1)
    endif()
  endif()
  
  if(NOT(WITH_OPENMP))
	set(NUM_THREADS 1)
  endif(NOT(WITH_OPENMP))
  
  
  MESSAGE(STATUS "Hermes will by default be using ${NUM_THREADS} threads for calculation.")

  # Enable debugging symbols for all files all the time:
  # set(CMAKE_CXX_FLAGS "-g ${CMAKE_CXX_FLAGS}")

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Disable all warnings and turn on only important ones:
    set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

    set(RELEASE_FLAGS "-DNDEBUG -O3")
    set(FLAGS ${RELEASE_FLAGS})
    set(DEBUG_FLAGS   "-g")
    if(DEBUG)
      set(FLAGS ${DEBUG_FLAGS})
    endif(DEBUG)
  endif()
  
  # Enabling multiprocessor build on MSVC
  if(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /openmp")
    SET(FLAGS "${FLAGS} /openmp")
  endif(MSVC)

  # This overrides CXX flags for MSVC
  if(MSVC)
    set(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dpopen=_popen /Dpclose=_pclose /D__value=_value /Dfinite=_finite /Dhypot=_hypot /Disatty=_isatty /Dfileno=_fileno /D_CRT_SECURE_NO_WARNINGS /DYY_NO_UNISTD_H /D_USE_MATH_DEFINES /DIMPLEMENT_C99 /wd4275 /wd4251")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT          "/D_DEBUG /Od /Ob2 /MDd /Zi /RTC1 ${MSVC_DEFINES}")
    set(DEBUG_FLAGS ${DEBUG_FLAGS} "${CMAKE_CXX_FLAGS_DEBUG_INIT}")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "/DNDEBUG /O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(RELEASE_FLAGS ${RELEASE_FLAGS} "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/DNDEBUG /O2 /Ob2 /MD /Zi ${MSVC_DEFINES}")
  else(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG_INIT          "${CMAKE_CXX_FLAGS_DEBUG_INIT} -D_DEBUG")
  endif(MSVC)


  if(H2D_WITH_TESTS OR HERMES_COMMON_WITH_TESTS)
    enable_testing()
    add_custom_target(test COMMAND ctest)
  endif(H2D_WITH_TESTS OR HERMES_COMMON_WITH_TESTS)

  #------------------------------------------------------------------------------
  # Find necessary packages according to the settings above.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    find_package(PTHREAD REQUIRED)
    find_package(XSD REQUIRED)
    find_package(XERCES REQUIRED)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      if(WITH_STACKTRACE)
        find_package(EXECINFO)
      endif(WITH_STACKTRACE)
    endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    
    # Solvers.
    if(WITH_UMFPACK)
      find_package(UMFPACK REQUIRED)
      include_directories(${UMFPACK_INCLUDE_DIRS})
    endif(WITH_UMFPACK)

    if(MSVC)
      if(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
          # The search for CLAPACK and associated libraries will be performed
          # every time the project is configured in MSVC (TODO: Find out if
          # CLAPACK_LIBRARY, BLAS_LIBRARY, F2C_LIBRARY and CLAPACK_INCLUDE_DIRS
          # could be cached to improve subsequent configuration performance).
          find_package(CLAPACK REQUIRED)
          include_directories(${CLAPACK_INCLUDE_DIRS})
      endif(WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    else(MSVC)
      # All third party solvers require BLAS and LAPACK (both found in the package LAPACK).
      # FindLAPACK is Fortran-only so Fortran must be enabled as well.
      if(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
        if (NOT LAPACK_FOUND)
            enable_language(Fortran)
            find_package(LAPACK REQUIRED)
                
            # If no error occured, LAPACK library has been found. Save the path to
            # it to cache, so that it will not be searched for during next 'cmake .'
            set(LAPACK_LIBRARIES  ${LAPACK_LIBRARIES}
                    CACHE STRING  "Path to LAPACK/BLAS libraries.")
            set(LAPACK_FOUND      YES
                    CACHE STRING  "Have LAPACK/BLAS libraries been found?")
        endif (NOT LAPACK_FOUND)
        
        add_definitions(-DWITH_BLAS)  
      endif(WITH_UMFPACK OR WITH_PETSC OR WITH_MUMPS OR WITH_TRILINOS OR WITH_SUPERLU)
    endif(MSVC)

    
    if(WITH_PETSC)
      find_package(PETSC REQUIRED)
      # For some header files, PETSc uses the same name for both real and complex
      # versions (they only reside in different directories), so we can't globally
      # set the architecture-specific include directory here, but rather on a
      # per-target basis.
      include_directories(${PETSC_INCLUDE_DIRS})
      include_directories(${PETSC_CPLX_INCLUDE_DIRS})
      set(PETSC_LIBRARIES ${PETSC_CPLX_LIBRARIES})
    endif(WITH_PETSC)

    if(WITH_MUMPS)
      find_package(MUMPS REQUIRED)
      if(WITH_MPI)
        find_package(SCALAPACK REQUIRED) # parallel version of MUMPS needs ScaLapack+BLACS
      endif(WITH_MPI)
      include_directories(${MUMPS_INCLUDE_DIR})
      set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})
      LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
    endif(WITH_MUMPS)

    if(WITH_TRILINOS)
      find_package(TRILINOS REQUIRED)
      include_directories(${TRILINOS_INCLUDE_DIR})
    endif(WITH_TRILINOS)
    
    # Experimental
    if(WITH_ZOLTAN)
      find_package(ZOLTAN REQUIRED)
      include_directories(${ZOLTAN_INCLUDE_DIR})
    endif(WITH_ZOLTAN)

    if(WITH_SUPERLU)
      find_package(SUPERLU REQUIRED)
      include_directories(${SUPERLU_INCLUDE_DIR})
      if(SUPERLU_MT)
        add_definitions(-DSLU_MT)
      endif(SUPERLU_MT)
    endif(WITH_SUPERLU)

    # SMP.
    if(WITH_OPENMP)
      find_package(OpenMP REQUIRED)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif(WITH_OPENMP)

    if(WITH_EXODUSII)
      find_package(EXODUSII REQUIRED)
      include_directories(${EXODUSII_INCLUDE_DIR})
    endif(WITH_EXODUSII)

    # If using any package that requires MPI (e.g. parallel versions of MUMPS, PETSC).
    if(WITH_MPI)
      if(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH) # If MPI was not defined by the user
        find_package(MPI REQUIRED)
      else(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
        if(NOT EXISTS MPI_INCLUDE_PATH)
          message(FATAL_ERROR "MPI include directory ${MPI_INCLUDE_PATH} not found")
        endif(NOT EXISTS MPI_INCLUDE_PATH)
      endif(NOT MPI_LIBRARIES OR NOT MPI_INCLUDE_PATH)
      include_directories(${MPI_INCLUDE_PATH})
    endif(WITH_MPI)
    
    
    IF(WITH_STACKTRACE)
      IF(MSVC)
        find_package(STACK_WALKER REQUIRED)
        include_directories(${STACK_WALKER_INCLUDE_DIR})
      ENDIF(MSVC)
    ENDIF(WITH_STACKTRACE)

  if(NOT H2D_WITH_GLUT)
    set(NOGLUT  YES)
  endif(NOT H2D_WITH_GLUT)

  # Create file with preprocessor definitions exposing the build settings to the source code.
  configure_file(
    config.h.cmake
    include/config.h
  )
  
  # Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
  # and set its build properties (compile and link flags, installation directories).
  macro(BUILD_LIB LIBRARY BUILD_FLAGS)
    add_library(${LIBRARY} SHARED ${SRC} ${HEADERS})
    set_target_properties(${LIBRARY} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})

    # Link with all third party libraries corresponding to selected configuration
    # options. Linking of optional libraries for the different Hermeses ({1, 2, 3}D)
    # is controled in their respective CMakeLists.txt. Note that for unselected
    # TPL's, the respective variables will be empty.

    target_link_libraries(${LIBRARY}
      ${EXODUSII_LIBRARIES}
      ${HDF5_LIBRARY}
      ${METIS_LIBRARY}
      ${UMFPACK_LIBRARIES}
      ${TRILINOS_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${MUMPS_LIBRARIES}
      ${SUPERLU_LIBRARY}
      ${LAPACK_LIBRARIES}
      ${STACK_WALKER_LIBRARY}
      ${PTHREAD_LIBRARY} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES}
      ${CLAPACK_LIBRARY} ${BLAS_LIBRARY} ${F2C_LIBRARY}
      ${ADDITIONAL_LIBS}
    )

    if(NOT MSVC)
      IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        target_link_libraries(${LIBRARY} "rt")
      ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    endif(NOT MSVC)
  endmacro(BUILD_LIB)
  
  include_directories(include)
  
  add_subdirectory(src/algebra)
  add_subdirectory(src/base)
  add_subdirectory(src/util)
    
  if(WITH_H2D)
    set(HERMES2D "hermes2d")              # Base name of the library.
    
    add_subdirectory(src/2d)

    IF(EXISTS "2d/test_examples")
      if(H2D_WITH_TEST_EXAMPLES)
        add_subdirectory(2d/test_examples)
      endif(H2D_WITH_TEST_EXAMPLES)
    ENDIF(EXISTS "2d/test_examples")
  endif(WITH_H2D)


  #------------------------------------------------------------------------------
  # Report.
  #~~~~~~~~
  message("\nConfiguration results")
  message("---------------------")
  message("Build with UMFPACK: ${WITH_UMFPACK}")
  message("Build with PETSC: ${WITH_PETSC}")
  message("Build with MUMPS: ${WITH_MUMPS}")
  message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
  message("Build with TRILINOS: ${WITH_TRILINOS}")
  message("Build with MPI: ${WITH_MPI}")
  message("Build with OPENMP: ${WITH_OPENMP}")
  message("Build with EXODUSII: ${WITH_EXODUSII}")
  
  message("---------------------")
  message("General info:")
  message("\tBuild Debug version: ${DEBUG}")
  message("\tBuild with tests: ${HERMES_COMMON_WITH_TESTS}")

  message("Build Hermes2D: ${WITH_H2D}")
  if(WITH_H2D)
    message("\tBuild Hermes2D with test examples: ${H2D_WITH_TEST_EXAMPLES}")
  message("---------------------")
    message("\tBuild Hermes2D with GLUT: ${H2D_WITH_GLUT}")
    message("\tBuild Hermes2D with VIEWER_GUI: ${H2D_WITH_VIEWER_GUI}")
  endif(WITH_H2D)
  message("\n")
