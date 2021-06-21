esma_set_this ()

set (srcs
  Carbon_GridCompMod.F90
  )

set (resource_files
   Carbon_GridComp.rc
   )

#mapl_acg (${this}   DU2G_StateSpecs.rc 
#          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS 
#          GET_POINTERS DECLARE_POINTERS)

esma_add_library (${this} 
   SRCS ${srcs} 
   DEPENDENCIES MAPL GEOS_Shared
   INCLUDES ${INC_ESMF})


file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

