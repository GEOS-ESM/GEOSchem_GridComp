esma_set_this ()

# These are not "SUBCOMPONENTS" because
# they provide more than a single module/library
esma_add_subdirectories(Shared GOCART)

# This set(rename_) is for the automatic stubber
# If directory is not available, this renames the
# module in the stubbed routine to be
# name which is expected by the parent.
set (rename_TR TR_GridComp)
set (rename_GMI GMIchem_GridComp)
set (rename_StratChem StratChem_GridComp)
set (rename_MAM MAM_GridComp)
set (rename_MATRIX MATRIX_GridComp)
set (rename_CARMA CARMAchem_GridComp)
set (rename_GAAS GAAS_GridComp)
set (rename_ACHEM ACHEM_GridComp)
set (rename_RRG RRG_GridComp)

set (alldirs
  GEOSpchem_GridComp
  GEOSCHEMchem_GridComp
  H2O_GridComp
  DNA_GridComp
  HEMCO_GridComp
  RRG
  TR
  GMI
  StratChem
  MAM
  MATRIX
  CARMA
  GAAS
  ACHEM
  )

set (srcs
  GEOS_ChemGridComp.F90
  GEOS_ChemEnvGridComp.F90
  )

# Note: GOCART_GridComp is added specifically as it is not longer
#       automatically a dependence through alldirs
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp ESMF::ESMF)

install(
   FILES GEOS_ChemGridComp.rc ChemEnv_ExtData.rc  ChemEnv.rc ChemEnv_ExtData.yaml
   DESTINATION etc
   )
