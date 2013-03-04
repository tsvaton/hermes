project(hermes2d)
  # Generating classes for XML mesh parsing according to your version of XSD.
  message("Hermes2D: Generating classes for XML mesh parsing according to your version of XSD.")
  if(WIN32)
    execute_process(COMMAND "${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_win.bat" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/" OUTPUT_QUIET)
  else(WIN32)
    execute_process(COMMAND "chmod 700 ${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_linux.sh" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/")
    execute_process(COMMAND "${PROJECT_SOURCE_DIR}/xml_schemas/xml_schemas_generate_linux.sh" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/xml_schemas/" OUTPUT_QUIET)
  endif(WIN32)

  SET(H2D_XML_SCHEMAS_DIRECTORY "\"${PROJECT_SOURCE_DIR}/xml_schemas\"")
  SET(H2D_PRECALCULATED_FORMS_DIRECTORY "\"${PROJECT_SOURCE_DIR}/precalculated_forms\"")

  if(MSVC)
    set(HERMES2D_RELEASE ${HERMES2D})     # Name of the release version.
    set(HERMES2D_DEBUG ${HERMES2D})       # Name of the debug version.
  else(MSVC)
    if(DEBUG)
      set(HERMES2D_DEBUG "hermes2d-debug")
      set(HERMES2D ${HERMES2D_DEBUG})
      # set(CMAKE_BUILD_TYPE Debug)
    # This does not work with multi-version configurations,
    # because it is linked to the project (i.e. we would have
    # to have separate projects for Debug and Release)
    else(DEBUG)
      set(HERMES2D_RELEASE "hermes2d")
      set(HERMES2D ${HERMES2D_RELEASE})
      # set(CMAKE_BUILD_TYPE Release)
    endif(DEBUG)
  endif(MSVC)

  #
  # Optional libraries and settings for H2D.
  #
  if(H2D_WITH_GLUT)
    if(NOT APPLE)   # See hermes/OSX.cmake for APPLE-specific search procedures.
      find_package(GLUT REQUIRED)
      find_package(GLEW REQUIRED)
    endif(NOT APPLE)
  else(H2D_WITH_GLUT)
    add_definitions(-DNOGLUT)
  endif(H2D_WITH_GLUT)

  if(H2D_WITH_VIEWER_GUI)
    find_package(ANTTWEAKBAR REQUIRED)
    include_directories(${ANTTWEAKBAR_INCLUDE_DIR})
  endif(H2D_WITH_VIEWER_GUI)
  
  # Mesh format.
  if(WITH_HDF5)
    find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIR})
  endif(WITH_HDF5)

  if(WITH_EXODUSII)
    find_package(EXODUSII REQUIRED)
    include_directories(${EXODUSII_INCLUDE_DIR})
  endif(WITH_EXODUSII)
  include_directories(${XSD_INCLUDE_DIR})
  include_directories(${XERCES_INCLUDE_DIR})
  
  #
  # Source files for the Hermes2d library.
  #
  set(SRC
    forms.cpp
    asmlist.cpp
    newton_solver.cpp
    picard_solver.cpp
    linear_solver.cpp
    
    calculation_continuity.cpp

    adapt/adapt.cpp
    adapt/kelly_type_adapt.cpp

    boundary_conditions/essential_boundary_conditions.cpp

    function/transformable.cpp
    function/function.cpp
    function/exact_solution.cpp
    function/solution.cpp
    function/filter.cpp
    function/mesh_function.cpp
    function/solution_h2d_xml.cpp

    api2d.cpp
    mixins2d.cpp
    
    mesh/refmap.cpp
    mesh/curved.cpp
    mesh/refinement_type.cpp
    mesh/element_to_refine.cpp
    mesh/mesh_reader_exodusii.cpp
    mesh/hash.cpp
    mesh/mesh_reader_h2d.cpp
    mesh/mesh_reader_h2d_xml.cpp
    mesh/mesh_reader_h1d_xml.cpp
    mesh/mesh_h2d_xml.cpp
    mesh/mesh_h1d_xml.cpp
    mesh/subdomains_h2d_xml.cpp
    mesh/mesh.cpp
    mesh/traverse.cpp
    mesh/mesh_data.cpp

    quadrature/limit_order.cpp
    quadrature/quad_std.cpp

    refinement_selectors/selector.cpp
    refinement_selectors/order_permutator.cpp
    refinement_selectors/optimum_selector.cpp
    refinement_selectors/proj_based_selector.cpp
    refinement_selectors/l2_proj_based_selector.cpp
    refinement_selectors/h1_proj_based_selector.cpp
    refinement_selectors/hcurl_proj_based_selector.cpp

    shapeset/shapeset.cpp
    shapeset/shapeset_h1_ortho.cpp
    shapeset/shapeset_h1_jacobi.cpp
    shapeset/shapeset_h1_quad.cpp
    shapeset/shapeset_hc_legendre.cpp
    shapeset/shapeset_hc_gradleg.cpp
    shapeset/shapeset_hd_legendre.cpp
    shapeset/shapeset_l2_legendre.cpp
    shapeset/precalc.cpp

    space/space.cpp
    space/space_h1.cpp
    space/space_hcurl.cpp
    space/space_l2.cpp
    space/space_hdiv.cpp
    space/space_h2d_xml.cpp

    views/base_view.cpp
    views/mesh_view.cpp
    views/order_view.cpp
    views/scalar_view.cpp
    views/stream_view.cpp
    views/vector_base_view.cpp
    views/vector_view.cpp
    views/view.cpp
    views/view_data.cpp
    views/view_support.cpp
    views/linearizer.cpp
    views/linearizer_base.cpp
    views/orderizer.cpp
    views/vectorizer.cpp

    weakform/weakform.cpp

    neighbor.cpp
    graph.cpp
    global.cpp
    discrete_problem.cpp
    discrete_problem_linear.cpp
    runge_kutta.cpp
    spline.cpp

    projections/ogprojection.cpp
    projections/ogprojection_nox.cpp
    projections/localprojection.cpp
    
    weakform_library/weakforms_elasticity.cpp
    weakform_library/weakforms_h1.cpp
    weakform_library/weakforms_hcurl.cpp
    weakform_library/weakforms_maxwell.cpp
    weakform_library/weakforms_neutronics.cpp
  )
  
  set(HEADERS
    ${CMAKE_HOME_DIRECTORY}/include/2d/forms.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/asmlist.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/newton_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/picard_solver.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/linear_solver.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/calculation_continuity.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/adapt/adapt.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/adapt/kelly_type_adapt.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/boundary_conditions/essential_boundary_conditions.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/function/transformable.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/function.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/exact_solution.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/solution.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/filter.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/mesh_function.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/function/solution_h2d_xml.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/integrals/l2.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/integrals/h1.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/integrals/hdiv.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/integrals/hcurl.h
    
    ${CMAKE_HOME_DIRECTORY}/include/2d/api2d.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mixins2d.h
    
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/refmap.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/curved.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/refinement_type.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/element_to_refine.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_reader_exodusii.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/hash.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_reader_h2d.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_reader_h2d_xml.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_reader_h1d_xml.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_h2d_xml.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_h1d_xml.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/subdomains_h2d_xml.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/traverse.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/mesh/mesh_data.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/quadrature/limit_order.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/quadrature/quad.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/quadrature/quad_all.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/selector.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/order_permutator.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/optimum_selector.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/proj_based_selector.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/l2_proj_based_selector.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/h1_proj_based_selector.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/hcurl_proj_based_selector.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset_h1_all.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset_h1_quad.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset_hc_all.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset_hd_all.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/shapeset_l2_all.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/precalc.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space_h1.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space_hcurl.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space_l2.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space_hdiv.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/space/space_h2d_xml.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/views/base_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/mesh_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/order_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/scalar_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/stream_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/vector_base_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/vector_view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/view.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/view_support.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/linearizer.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/linearizer_base.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/orderizer.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/views/vectorizer.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform/weakform.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/neighbor.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/graph.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/global.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/discrete_problem.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/discrete_problem_linear.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/runge_kutta.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/spline.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/projections/ogprojection.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/projections/ogprojection_nox.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/projections/localprojection.h

    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/weakforms_elasticity.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/weakforms_h1.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/weakforms_hcurl.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/weakforms_maxwell.h
    ${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/weakforms_neutronics.h
  )
    
  #
  # Build and install the library.
  #
  include(BuildAndInstallScripts)

  # Add a target for a specified version of the Hermes library,
  # and set its build properties (compile and link flags, installation directories). Due to
  # different optional libraries for H1D, H2D and H3D, this macro is dimension-specific.
  macro(BUILD_2D_LIB HERMES_LIB BUILD_FLAGS)
    add_library(${HERMES_LIB} SHARED ${SRC} ${HEADERS})

    set(_FLAGS ${BUILD_FLAGS})

    set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${_FLAGS})

    target_link_libraries(  ${HERMES_LIB}
      util base algebra
      ${GLUT_LIBRARY} ${GLEW_LIBRARY}
      ${ANTTWEAKBAR_LIBRARY}
      ${XSD_LIBRARY}
      ${XERCES_LIBRARY}
      ${LAPACK_LIBRARY}
      ${CLAPACK_LIBRARY} ${BLAS_LIBRARY}
    )
  endmacro(BUILD_2D_LIB)

  # Build and install:
  if(MSVC)
    if(DEBUG)
      set(FLAGS "${DEBUG_FLAGS}")
    elseif(DEBUG)
      set(FLAGS "${RELEASE_FLAGS}")
    endif(DEBUG)

    BUILD_2D_LIB(${HERMES2D} ${FLAGS})
    set_target_properties(${HERMES2D} PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${FLAGS}")
    INSTALL_LIB(${HERMES2D})
  else(MSVC)
    if(DEBUG)
      BUILD_2D_LIB(${HERMES2D_DEBUG} ${HERMES_COMMON_LIB_DEBUG} ${DEBUG_FLAGS})
      INSTALL_LIB(${HERMES2D_DEBUG})
    else(DEBUG)
      BUILD_2D_LIB(${HERMES2D_RELEASE} ${HERMES_COMMON_LIB_RELEASE} ${RELEASE_FLAGS})
      INSTALL_LIB(${HERMES2D_RELEASE})
    endif(DEBUG)
  endif(MSVC)

  # Install header files.

  file(GLOB INC            "${CMAKE_HOME_DIRECTORY}/include/2d/*.h")
  file(GLOB INC_FORM       "${CMAKE_HOME_DIRECTORY}/include/2d/form/*.h")
  file(GLOB INC_SPACE      "${CMAKE_HOME_DIRECTORY}/include/2d/space/*.h")
  file(GLOB INC_MESH       "${CMAKE_HOME_DIRECTORY}/include/2d/mesh/*.h")
  file(GLOB INC_XSD        "xml_schemas/*.xsd")
  file(GLOB INC_BC         "${CMAKE_HOME_DIRECTORY}/include/2d/boundary_conditions/*.h")
  file(GLOB INC_SHAPESET   "${CMAKE_HOME_DIRECTORY}/include/2d/shapeset/*.h")
  file(GLOB INC_FUNCTION   "${CMAKE_HOME_DIRECTORY}/include/2d/function/*.h")
  file(GLOB INC_INTEGRALS  "${CMAKE_HOME_DIRECTORY}/include/2d/integrals/*.h")
  file(GLOB INC_LINEARIZER "${CMAKE_HOME_DIRECTORY}/include/2d/linearizer/*.h")
  file(GLOB INC_ADAPT      "${CMAKE_HOME_DIRECTORY}/include/2d/adapt/*.h")
  file(GLOB INC_PROJECTIONS "${CMAKE_HOME_DIRECTORY}/include/2d/projections/*.h")
  file(GLOB INC_PRECALC     "${CMAKE_HOME_DIRECTORY}/include/2d/precalculated_forms/*.h")
  file(GLOB INC_QUAD       "${CMAKE_HOME_DIRECTORY}/include/2d/quadrature/*.h")
  file(GLOB INC_REF_SEL    "${CMAKE_HOME_DIRECTORY}/include/2d/refinement_selectors/*.h")
  file(GLOB INC_VIEWS      "${CMAKE_HOME_DIRECTORY}/include/2d/views/*.h")
  file(GLOB INC_WEAKFORM   "${CMAKE_HOME_DIRECTORY}/include/2d/weakform/*.h")
  file(GLOB INC_WEAKFORM_LIBRARY "${CMAKE_HOME_DIRECTORY}/include/2d/weakform_library/*.h")
  
  install(FILES ${INC}           DESTINATION ${TARGET_ROOT}/include/2d)
  install(FILES ${INC_SPACE}     DESTINATION ${TARGET_ROOT}/include/2d/space)
  install(FILES ${INC_MESH}      DESTINATION ${TARGET_ROOT}/include/2d/mesh)
  install(FILES ${XSD}           DESTINATION ${TARGET_ROOT}/include/2d/xml_schemas)
  install(FILES ${INC_BC}        DESTINATION ${TARGET_ROOT}/include/2d/boundary_conditions)
  install(FILES ${INC_SHAPESET}  DESTINATION ${TARGET_ROOT}/include/2d/shapeset)
  install(FILES ${INC_FUNCTION}  DESTINATION ${TARGET_ROOT}/include/2d/function)
  install(FILES ${INC_INTEGRALS} DESTINATION ${TARGET_ROOT}/include/2d/integrals)
  install(FILES ${INC_LINEARIZER}DESTINATION ${TARGET_ROOT}/include/2d/linearizer)
  install(FILES ${INC_ADAPT}     DESTINATION ${TARGET_ROOT}/include/2d/adapt)
  install(FILES ${INC_PROJECTIONS}  DESTINATION ${TARGET_ROOT}/include/2d/projections)
  install(FILES ${INC_PRECALC}  DESTINATION ${TARGET_ROOT}/include/2d/precalculated_forms)
  install(FILES ${INC_QUAD}      DESTINATION ${TARGET_ROOT}/include/2d/quadrature)
  install(FILES ${INC_REF_SEL}   DESTINATION ${TARGET_ROOT}/include/2d/refinement_selectors)
  install(FILES ${INC_VIEWS}     DESTINATION ${TARGET_ROOT}/include/2d/views)
  install(FILES ${INC_WEAKFORM}  DESTINATION ${TARGET_ROOT}/include/2d/weakform)
  install(FILES ${INC_WEAKFORM_LIBRARY}  DESTINATION ${TARGET_ROOT}/include/2d/weakform_library)

