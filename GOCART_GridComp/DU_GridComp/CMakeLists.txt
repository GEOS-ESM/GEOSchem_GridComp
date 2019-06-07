esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES Chem_Shared Chem_Base GMAO_mpeu)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

esma_generate_gocart_code (${this} "-B\;-E\;-C\;-N\;GOCART")

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach(rc_file ${rc_files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${rc_file} ${esma_etc}/${rc_file} COPYONLY)
endforeach()
