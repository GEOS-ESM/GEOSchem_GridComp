esma_set_this ()

add_subdirectory(Shared)

set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)

set (alldirs
  GEOSpchem_GridComp
  GOCART_GridComp
  StratChem_GridComp
  GMIchem_GridComp
  CARMAchem_GridComp
  GEOSCHEMchem_GridComp 
  MATRIXchem_GridComp
  MAMchem_GridComp 
  GAAS_GridComp
  H2O_GridComp
  TR_GridComp 
  GEOSachem_GridComp
  DNA_GridComp
  @HEMCO_GridComp
  )

set (srcs
  GEOS_ChemGridComp.F90
  GEOS_ChemEnvGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL Chem_Shared Chem_Base
  INCLUDES ${INC_ESMF})

install(
   FILES GEOS_ChemGridComp.rc
   DESTINATION etc
   )
