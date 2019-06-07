esma_set_this ()

set (srcs
  DNA_GridCompMod.F90
  )

include_directories (${include_MAPL_Base})


esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)
target_compile_definitions (${this} PRIVATE MAPL GEOS5)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})


