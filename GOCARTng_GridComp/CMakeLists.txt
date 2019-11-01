esma_set_this ()

set (alldirs
  DUng_GridComp
  SSng_GridComp
#  SUng_GridComp
  BCng_GridComp
  OCng_GridComp
  )

set (srcs
  GOCART_GridCompMod.F90
  )

set (resource_files
   GOCARTng_GridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (dependencies Chem_Base Chem_Shared MAPL_Base GMAO_mpeu)
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES ${dependencies}
  INCLUDES ${INC_ESMF})

