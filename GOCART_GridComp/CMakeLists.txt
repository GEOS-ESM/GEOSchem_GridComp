esma_set_this ()

set (alldirs
  O3_GridComp
  CO_GridComp
  CO2_GridComp
  CFC_GridComp
  DU_GridComp
  SS_GridComp  
  BC_GridComp
  OC_GridComp
  SU_GridComp
  Rn_GridComp
  CH4_GridComp
  NI_GridComp
  BRC_GridComp
  )

set (srcs
  GOCART_GridCompMod.F90
  Aero_GridCompMod.F90
  )

set (dependencies Chem_Base Chem_Shared MAPL_Base GMAO_mpeu)
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES ${dependencies}
  INCLUDES ${INC_ESMF})

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach(rc_file ${rc_files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${rc_file} ${esma_etc}/${rc_file} COPYONLY)
endforeach()

