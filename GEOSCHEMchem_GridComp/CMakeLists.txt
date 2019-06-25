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
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GIGC_GridCompMod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/Chem_GridCompMod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GEOS_ctmEnvGridComp.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/GIGC/GEOSChem.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/gc_column/Headers/inquireMod.F90
   )

esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL_Base MAPL_cfio_r4 Chem_Base Chem_Shared HEMCO
   )

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

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


set (acg_headers
   GEOSCHEMchem_InternalSpec___.h 
   GEOSCHEMchem_DeclarePointer___.h 
   GEOSCHEMchem_GetPointer___.h
   GEOSCHEMchem_History___.rc
  )

add_custom_command (
  OUTPUT ${acg_headers}
  COMMAND ${esma_include}/MAPL_Base/mapl_acg.pl ${acg_flags} -v -F ${CMAKE_CURRENT_SOURCE_DIR}/GEOSCHEMchem_Registry.rc
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  MAIN_DEPENDENCY GEOSCHEMchem_Registry.rc 
  COMMENT "Generating automated code for GEOSCHEMchem"
  )
add_custom_target (phony_GEOSCHEMchem DEPENDS ${acg_headers})
add_dependencies (${this} phony_GEOSCHEMchem)

file (GLOB resource_files "*.rc")

install(
   FILES ${resource_files}
   DESTINATION etc
   )
