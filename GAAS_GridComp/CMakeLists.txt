esma_set_this ()

set (srcs
#  m_ana.F90
  LDE_Mod.F90
  GAAS_GridCompMod.F90
  )


esma_add_library (${this} SRCS ${srcs} DEPENDENCIES MAPL Chem_Base Chem_Shared GMAO_mpeu FVdycoreCubed_GridComp)
target_include_directories (${this} PUBLIC ${INC_NETCDF})

esma_generate_gocart_code (${this} -F)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

install(PROGRAMS ana_lde.py DESTINATION bin)

ecbuild_add_executable (TARGET ana_lde.x SOURCES ana_lde.F90  LIBS ${this})

