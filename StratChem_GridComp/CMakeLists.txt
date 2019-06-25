esma_set_this ()


set (subsrcs
  SC_GridCompMod.F90
  brpandl.F  brpandlest.F
  calc_qqjk.F
  calcsulf.F  chxpart.F
  density.F
  fallvel.F
  gamma_shi.F getBaseChem.F
  hoxpart.F hoxpartest.F
  interp_s.F
  jcalc4.F
  khet3d.F  krates.F
  localZeroSet.F90
  meso_phot.F
  ncpandl.F ncpandlest.F
  o3colum.F oxpandl.F oxpandlest.F
  part.F partest.F pscs.F
  rdkrate.F
  sadice.F sadnat.F sadsts.F sediment.F solverd2.F solverest.F
  sproic.F  storeBaseChem.F sza.F
  ternary.F
  verifyChemDt.F
  )

set (srcs)
foreach (src ${subsrcs})
  list (APPEND srcs SC_GridComp/${src})
endforeach ()
list (APPEND srcs StratChem_GridCompMod.F90)

# Workaround for several older files that exceed the
# 72 column limit of fixed-format Fortran source code.
foreach (file
    brpandl.F brpandlest.F
    calc_qqjk.F
    chxpart.F getBaseChem.F
    hoxpart.F hoxpartest.F
    interp_s.F khet3d.F
    krates.F meso_phot.F
    ncpandl.F ncpandlest.F
    oxpandl.F oxpandlest.F
    part.F partest.F
    pscs.F rdkrate.F
    sediment.F
    solverd2.F solverest.F
    sproic.F
    storeBaseChem.F
    verifyChemDt.F
)
  if (EXTENDED_SOURCE)
    set_source_files_properties (SC_GridComp/${file} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  endif ()
endforeach ()

include_directories (${include_Chem_Shared})
include_directories (${include_MAPL_Base})
include_directories (${include_MAPL_cfio_r4})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${include_GMAO_mpeu})

set (dependencies Chem_Shared MAPL_Base MAPL_cfio_r4 GMAO_mpeu)
esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
target_include_directories (${this} PUBLIC ${INC_NETCDF})

string (REPLACE "_GridComp" "" name ${this})

set (automatic_headers
  ${name}_ExportSpec___.h
  ${name}_ImportSpec___.h
  ${name}_GetPointer___.h
  )
set (automatic_rc
  ${name}_History___.rc
  )

set (registry SC_GridComp/SC_Registry.rc)

new_esma_generate_automatic_code (
  ${this} ${registry}
  "${automatic_headers}" "${automatic_rc}"
  ${include_GEOSchem_GridComp} ${esma_etc}
  "-N STRATCHEM"
  )


set (resource_files
   SC_GridComp/SC_GridComp.rc
   SC_GridComp/SC_Registry.rc
   SC_GridComp/StratChem_ExtData.rc
   )

install(
   FILES ${resource_files}
   DESTINATION etc
   )
