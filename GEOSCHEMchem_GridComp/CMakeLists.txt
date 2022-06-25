esma_set_this ()

set (alldirs
  )

esma_mepo_style(geos-chem geos_chem_path)
set (geos_chem_dir ${CMAKE_CURRENT_SOURCE_DIR}/${geos_chem_path})


set (src_directories
#   ${CMAKE_CURRENT_SOURCE_DIR}
   ${geos_chem_dir}/Headers
   ${geos_chem_dir}/History
   ${geos_chem_dir}/GeosUtil
   ${geos_chem_dir}/NcdfUtil
   ${geos_chem_dir}/ISORROPIA
   ${geos_chem_dir}/KPP/fullchem
   ${geos_chem_dir}/GeosCore
   ${geos_chem_dir}/Interfaces/GCHP
   ${geos_chem_dir}/ObsPack
  )
include_directories (${geos_chem_dir}/Headers)
include_directories (${geos_chem_dir}/Interfaces/GEOS)

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs CONFIGURE_DEPENDS ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()
# Unfortunately, Intel always searches for include files in the directory of the source file, so if we
# leave the grid comp in place it will find the GCHP include files.  So we copy to the bin dir.
configure_file(${geos_chem_dir}/Interfaces/GCHP/Chem_GridCompMod.F90 Chem_GridCompMod.F90 COPYONLY)
list (APPEND srcs ${CMAKE_CURRENT_BINARY_DIR}/Chem_GridCompMod.F90)

list (REMOVE_ITEM srcs
   # this one is moved to HEMCO
   #${geos_chem_dir}/GeosUtil/regrid_a2a_mod.F90
   ${geos_chem_dir}/NcdfUtil/charpak_mod.F90
   #${geos_chem_dir}/Headers/charpak_mod.F90
   #${geos_chem_dir}/NcdfUtil/julday_mod.F90
   #${geos_chem_dir}/GeosUtil/julday_mod.F90

   # not sure about the rest
   ${geos_chem_dir}/GeosCore/main.F90
   ${geos_chem_dir}/GeosCore/tpcore_fvdas_mod.F90
   ${geos_chem_dir}/GeosCore/tpcore_geos57_window_mod.F90
   ${geos_chem_dir}/GeosCore/tpcore_geos5_window_mod.F90
   ${geos_chem_dir}/GeosCore/transport_mod.F90
   ${geos_chem_dir}/GeosCore/rrtmg_rad_transfer_mod.F90
   ${geos_chem_dir}/GeosCore/tomas_mod.F90
   ${geos_chem_dir}/GeosCore/exchange_mod.F90
   ${geos_chem_dir}/GeosCore/get_met_mod.F90
   ${geos_chem_dir}/GeosCore/flexgrid_read_mod.F90
   ${geos_chem_dir}/GeosCore/aero_drydep.F90
   ${geos_chem_dir}/GeosCore/exchange_mod.F90
   ${geos_chem_dir}/GeosCore/YuIMN_Code.F90
   ${geos_chem_dir}/GeosCore/apm_driv_mod.F90
   ${geos_chem_dir}/Interfaces/GCHP/Chem_GridCompMod.F90

   ${geos_chem_dir}/NcdfUtil/julday_mod.F
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   if (CMAKE_BUILD_TYPE MATCHES Release)
      set (CMAKE_Fortran_FLAGS_RELEASE "")
      foreach (src ${srcs})
         if (${src} MATCHES GEOSCHEMchem_GridCompMod.F90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${FOPT2} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
         else ()
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_FLAGS_RELEASE} ${EXTENDED_SOURCE}")
         endif ()
      endforeach ()
   endif ()
endif ()

esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL MAPL_cfio_r4 Chem_Base Chem_Shared HEMCO
   )

# Enable OpenMP
target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

target_compile_definitions(${this} PRIVATE
   ESMF_ UCX EXTERNAL_GRID EXTERNAL_FORCING USE_REAL8 MODEL_GEOS)

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   target_compile_definitions(${this} PRIVATE LINUX_IFORT)
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
   target_compile_definitions(${this} PRIVATE LINUX_PGI)
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   target_compile_definitions(${this} PRIVATE LINUX_GFORTRAN)
endif ()


set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

new_esma_generate_automatic_code (
   ${this} GEOSCHEMchem_Registry.rc
   "GEOSCHEMCHEM_InternalSpec___.h;GEOSCHEMCHEM_DeclarePointer___.h;GEOSCHEMCHEM_GetPointer___.h"
   GEOSCHEMCHEM_History___.rc
   ${include_GEOSCHEMchem_GridComp} ${esma_etc}
   "-v;-F"
   )

file (GLOB resource_files CONFIGURE_DEPENDS "*.rc")

install(
   FILES ${resource_files}
   DESTINATION etc
   )
