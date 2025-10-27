esma_set_this ()

esma_add_subdirectories(Shared GOCART TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM)

set (alldirs
  GEOSpchem_GridComp
  GEOSCHEMchem_GridComp
  H2O_GridComp
  DNA_GridComp
  HEMCO_GridComp
  RRG
  )

set (srcs
  GEOS_ChemGridComp.F90
  GEOS_ChemEnvGridComp.F90
  )

if (BUILD_PYCHEM_INTERFACE)
  list (APPEND srcs
  pyChem/pyChem/interface/cffi_lib/pychem_interface.f90
  pyChem/pyChem/interface/cffi_lib/pychem_interface.c)
endif ()

if (BUILD_PYCHEM_INTERFACE)

  message(STATUS "Building pyChem interface")

  add_definitions(-DPYCHEM_INTEGRATION)
  
  # The Python library creation requires mpiexec/mpirun to run on a
  # compute node. Probably a weird SLURM thing?
  find_package(MPI REQUIRED)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  # Set up some variables in case names change
  set(PYCHEM_INTERFACE_LIBRARY     ${CMAKE_CURRENT_BINARY_DIR}/libpyChem_interface_py${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(PYCHEM_INTERFACE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/pychem_interface_py.h)
  set(PYCHEM_INTERFACE_FLAG_HEADER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/pyChem/pyChem/interface/cffi_lib/pychem_interface.h)
  set(PYCHEM_INTERFACE_SRCS        ${CMAKE_CURRENT_SOURCE_DIR}/pyChem/pyChem/interface/cffi_lib/interface.py)

  # This command creates the shared object library from Python
  add_custom_command(
    OUTPUT ${PYCHEM_INTERFACE_LIBRARY}
    # Note below is essentially:
    #  mpirun -np 1 python file
    # but we use the CMake options as much as we can for flexibility
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYCHEM_INTERFACE_FLAG_HEADER_FILE} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python3_EXECUTABLE} ${PYCHEM_INTERFACE_SRCS}
    BYPRODUCTS ${PYCHEM_INTERFACE_HEADER_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${PYCHEM_INTERFACE_SRCS}
    COMMENT "Building pyChem interface library with Python"
    VERBATIM
  )

  # This creates a target we can use for dependencies and post build
  add_custom_target(generate_pychem_interface_library DEPENDS ${PYCHEM_INTERFACE_LIBRARY})

  # Because of the weird hacking of INTERFACE libraries below, we cannot
  # use the "usual" CMake calls to install() the .so. I think it's because
  # INTERFACE libraries don't actually produce any artifacts as far as
  # CMake is concerned. So we add a POST_BUILD custom command to "install"
  # the library into install/lib
  add_custom_command(TARGET generate_pychem_interface_library
    POST_BUILD
    # We first need to make a lib dir if it doesn't exist. If not, then
    # the next command can copy the script into a *file* called lib because
    # of a race condition (if install/lib/ isn't mkdir'd first)
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
    # Now we copy the file (if different...though not sure if this is useful)
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PYCHEM_INTERFACE_LIBRARY}" ${CMAKE_INSTALL_PREFIX}/lib
    )

  # We use INTERFACE libraries to create a sort of "fake" target library we can use
  # to make libFVdycoreCubed_GridComp.a depend on. It seems to work!
  add_library(pychem_interface_py INTERFACE)

  # The target_include_directories bits were essentially stolen from the esma_add_library
  # code...
  target_include_directories(pychem_interface_py INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # stubs
    # modules and copied *.h, *.inc
    $<BUILD_INTERFACE:${esma_include}/${this}>
    $<INSTALL_INTERFACE:include/${this}>
    )
  target_link_libraries(pychem_interface_py INTERFACE ${PYCHEM_INTERFACE_LIBRARY})

  # This makes sure the library is built first
  add_dependencies(pychem_interface_py generate_pychem_interface_library)

  # This bit is to resolve an issue and Google told me to do this. I'm not
  # sure that the LIBRARY DESTINATION bit actually does anything since
  # this is using INTERFACE
  install(TARGETS pychem_interface_py
    EXPORT ${PROJECT_NAME}-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

endif ()


# Note: GOCART_GridComp is added specifically as it is not longer
#       automatically a dependence through alldirs
if (BUILD_PYMOIST_INTERFACE)
  esma_add_library (${this}
    SRCS ${srcs}
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES pychem_interface_py MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM RRG ESMF::ESMF)
else ()
  esma_add_library (${this}
    SRCS ${srcs}
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM RRG ESMF::ESMF)
endif ()

install(
   FILES GEOS_ChemGridComp.rc ChemEnv_ExtData.rc  ChemEnv.rc ChemEnv_ExtData.yaml
   DESTINATION etc
   )
