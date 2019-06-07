esma_set_this ()

set (srcs
  TR_GravSetMod.F90
  DryDepositionGmiMod.F90
  WetRemovalGmiMod.F90
  TR_GridCompMod.F90
  )

esma_add_library (${this} 
   SRCS ${srcs} 
   DEPENDENCIES MAPL_Base GMAO_mpeu Chem_Shared GEOS_Shared
   INCLUDES ${INC_ESMF} ${INC_NETCDF})

#esma_generate_gocart_code (${this} -F)

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach(rc_file ${rc_files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${rc_file} ${esma_etc}/${rc_file} COPYONLY)
endforeach()
