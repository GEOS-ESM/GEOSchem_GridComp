esma_set_this ()

esma_add_subdirectories(Shared GOCART TR GMI StratChem)

set (alldirs
  GEOSpchem_GridComp
  CARMAchem_GridComp
  GEOSCHEMchem_GridComp
  MATRIXchem_GridComp
  MAMchem_GridComp
  GAAS_GridComp
  H2O_GridComp
  GEOSachem_GridComp
  DNA_GridComp
  HEMCO_GridComp
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
  DEPENDENCIES MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp TR GMI StratChem esmf)

install(
   FILES GEOS_ChemGridComp.rc ChemEnv_ExtData.rc  ChemEnv.rc ChemEnv_ExtData.yaml
   DESTINATION etc
   )
