
#############################################################
# If using in GEOSgcm, uncomment lines below
#############################################################
esma_set_this ()

esma_mepo_style(HEMCO hemco_dir)

set (src_directories
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Shared/Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Shared/NcdfUtil
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Shared/GeosUtil
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Core
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Extensions
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Interfaces/Shared
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Interfaces/Standalone
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Interfaces/MAPL_ESMF
   )

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs CONFIGURE_DEPENDS ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()

# These files contain Fortran MODULE statements that are
# deactivated by CPP statements.  CMake incorrectly determines
# that the file produces the corresponding .mod file, but they do not.
# Essentially a clash between using CPP to deactivate a file vs
# the config system.  The latter is preferable.  TLC
list (REMOVE_ITEM srcs
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Extensions/hcox_tomas_dustdead_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Extensions/hcox_tomas_seasalt_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Extensions/hcox_tomas_jeagle_mod.F90
   )

# The glob above brings in two nigh-identical hco_julday_mod
# modules. GNU make seems to allow this, but ninja does not. So we only
# build one
list (REMOVE_ITEM srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/${hemco_dir}/src/Shared/NcdfUtil/julday_mod.F90
  )

esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran
   )

target_compile_definitions (${this} PRIVATE ESMF_ USE_REAL8 MODEL_GEOS)
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})
