esma_set_this ()

set (srcs
  TR_GravSetMod.F90
  DryDepositionGmiMod.F90
  WetRemovalGmiMod.F90
  TR_GridCompMod.F90
  )

esma_add_library (${this} 
   SRCS ${srcs} 
   DEPENDENCIES MAPL GMAO_mpeu Chem_Shared GEOS_Shared
   INCLUDES ${INC_ESMF})

#esma_generate_gocart_code (${this} -F)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
