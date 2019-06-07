esma_set_this ()

set (srcs
#  m_ana.F90
  LDE_Mod.F90
  GAAS_GridCompMod.F90
  )


esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL_Base Chem_Base Chem_Shared GMAO_mpeu FVdycoreCubed_GridComp)
target_include_directories (${this} PUBLIC ${INC_NETCDF})

esma_generate_gocart_code (${this} -F)

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach(rc_file ${rc_files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${rc_file} ${esma_etc}/${rc_file} COPYONLY)
endforeach()
