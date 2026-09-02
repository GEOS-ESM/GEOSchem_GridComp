esma_set_this ()

set (srcs
  DNA_GridCompMod.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)
target_compile_definitions (${this} PRIVATE MAPL_MODE GEOS5)

set (resource_files
   DNA_ExtData.rc
   DNA_ExtData.yaml
   DNA_GridComp.rc
   DNA_Registry.rc
   )

install(
   FILES ${resource_files}
   DESTINATION etc
   )
