esma_set_this ()

set (srcs
   HEMCO_GridCompMod.F90
   )

esma_add_library (
  ${this}
  SRCS ${srcs}
  DEPENDENCIES Chem_Shared MAPL HEMCO esmf NetCDF::NetCDF_Fortran
  )
target_compile_definitions (${this} PRIVATE ESMF_ DEVEL GEOS_FP)

set (acg_headers
  HEMCO_InternalSpec___.h HEMCO_DeclarePointer___.h HEMCO_GetPointer___.h
  )


find_file (generator
  NAME mapl_acg.pl
  PATHS ${esma_include}/MAPL_Base ${esma_etc}/MAPL)


add_custom_command (
  OUTPUT ${acg_headers}
  COMMAND ${generator} ${acg_flags} -v -F ${CMAKE_CURRENT_SOURCE_DIR}/HEMCO_Registry.rc
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  MAIN_DEPENDENCY HEMCO_Registry.rc
#  DEPENDS ${include_MAPL_Base}/mapl_acg.pl
  COMMENT "Generating automated code for HEMCO"
  )
add_custom_target (phony_HEMCO DEPENDS ${acg_headers})
add_dependencies (${this} phony_HEMCO)

file (GLOB resource_files CONFIGURE_DEPENDS "*.rc")

install(
   FILES ${resource_files}
   DESTINATION etc
   )

file (GLOB yaml_files CONFIGURE_DEPENDS "*.yaml")

install(
   FILES ${yaml_files}
   DESTINATION etc
   )
