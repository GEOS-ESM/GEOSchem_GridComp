esma_set_this ()

set (alldirs
  )

set (src_directories
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/History
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosUtil
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/NcdfUtil
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/ISOROPIA
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/KPP/Standard
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/ObsPack
  )

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/gc_column/Headers)

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()

list (REMOVE_ITEM srcs
   # this one is moved to HEMCO
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosUtil/regrid_a2a_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/NcdfUtil/charpak_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/Headers/charpak_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/NcdfUtil/julday_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosUtil/julday_mod.F

   # not sure about the rest
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/main.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/tpcore_fvdas_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/tpcore_geos57_window_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/tpcore_geos5_window_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/transport_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/rrtmg_rad_transfer_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/tomas_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/exchange_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/get_met_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/flexgrid_read_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/aero_drydep.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/exchange_mod.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GeosCore/YuIMN_Code.F
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/gigc_manualexports_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GIGC_GridCompMod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/Chem_GridCompMod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GEOS_ctmEnvGridComp.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GEOSChem.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/Headers/inquireMod.F90
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
   DEPENDENCIES MAPL MAPL_cfio_r4 Chem_Base Chem_Shared
   )

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
target_link_libraries(${this} PUBLIC HCOI_MAPL_ESMF)

target_compile_definitions(${this} PRIVATE ESMF_ DEVEL NO_BPCH UCX 
   EXTERNAL_GRID EXTERNAL_FORCING GEOS_FP USE_REAL8 USE_TEND NC_DIAG DISCOVER MODEL_GEOS)

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

file (GLOB resource_files "*.rc")

install(
   FILES ${resource_files}
   DESTINATION etc
   )
