esma_set_this ()

set (srcs
  DNA_GridCompMod.F90
  )

include_directories (${include_MAPL_Base})


esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)
target_compile_definitions (${this} PRIVATE MAPL GEOS5)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

set (resource_files
   DNA_ExtData.rc
   DNA_GridComp.rc
   DNA_Registry.rc
   )

install(
   FILES ${resource_files}
   DESTINATION etc
   )
