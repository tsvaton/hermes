project(algebra)
    
  # Check for stacktrace and math library.
  include(CheckFunctionExists)
  CHECK_FUNCTION_EXISTS(fmemopen HAVE_FMEMOPEN)
  CHECK_FUNCTION_EXISTS(log2 HAVE_LOG2)

	# Additional building steps required by MSVC.
  include(BuildAndInstallScripts)

  #
  # Source files for the Hermes_common library.
  #
  set(SRC
    matrix.cpp
    tables.cpp
    solvers/dp_interface.cpp
    solvers/linear_matrix_solver.cpp
    solvers/nonlinear_solver.cpp
    solvers/newton_solver_nox.cpp
    solvers/epetra.cpp
    solvers/aztecoo_solver.cpp
    solvers/amesos_solver.cpp
    solvers/mumps_solver.cpp
    solvers/superlu_solver.cpp
    solvers/superlu_solver_real.cpp
    solvers/superlu_solver_cplx.cpp
    solvers/petsc_solver.cpp
    solvers/umfpack_solver.cpp
    solvers/precond_ml.cpp
    solvers/precond_ifpack.cpp
  )

  set(HEADERS
    ${CMAKE_HOME_DIRECTORY}/include/algebra/matrix.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/tables.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/dp_interface.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/linear_matrix_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/nonlinear_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/newton_solver_nox.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/epetra.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/aztecoo_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/amesos_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/mumps_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/superlu_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/petsc_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/umfpack_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/precond_ml.h
    ${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/precond_ifpack.h
  )
  
  #
  # Configure building and installation of all versions of hermes_common needed by
  # the selected Hermes configurations.
  #
  if(MSVC)
    BUILD_LIB(algebra ${DEBUG_FLAGS})
    target_link_libraries(algebra util)
    set_target_properties(algebra PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${FLAGS}")
    INSTALL_LIB(algebra)
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${DEBUG_FLAGS}")
    endif(DEBUG)
    BUILD_COMMON_LIB(algebra ${BUILD_FLAGS})
    INSTALL_LIB(algebra)
  endif(MSVC)

  # Install header files.

  file(GLOB INC_COMMON    "${CMAKE_HOME_DIRECTORY}/include/algebra/*.h")
  file(GLOB INC_SOLVER    "${CMAKE_HOME_DIRECTORY}/include/algebra/solvers/*.h")
  install(FILES ${INC_COMMON}    DESTINATION ${TARGET_ROOT}/include/algebra)
  install(FILES ${INC_SOLVER}    DESTINATION ${TARGET_ROOT}/include/algebra/solvers)
