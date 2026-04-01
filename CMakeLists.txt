esma_set_this ()

# Keep Shared (Chem_Base, Chem_Shared) and GOCART (for GOCART2G_GridComp).
# Remaining subdirs not yet ported to MAPL3 - restore each when ported:
esma_add_subdirectories(Shared GOCART)
# esma_add_subdirectories(Shared GOCART TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM)

# Temporarily bypass GEOSchem top-level library: GEOS_ChemGridComp.F90
# USE-s non-ported components. Restore the block below when all SUBCOMPONENTS
# are ported to MAPL3.

# set (alldirs
#   GEOSpchem_GridComp     # not yet ported to MAPL3 - restore when ported
#   GEOSCHEMchem_GridComp  # not yet ported to MAPL3 - restore when ported
#   H2O_GridComp           # not yet ported to MAPL3 - restore when ported
#   DNA_GridComp           # not yet ported to MAPL3 - restore when ported
#   HEMCO_GridComp         # not yet ported to MAPL3 - restore when ported
#   RRG                    # not yet ported to MAPL3 - restore when ported
#   )

# set (srcs
#   GEOS_ChemGridComp.F90
#   GEOS_ChemEnvGridComp.F90
#   )

# Note: GOCART_GridComp is added specifically as it is not longer
#       automatically a dependence through alldirs
# esma_add_library (${this}
#   SRCS ${srcs}
#   SUBCOMPONENTS ${alldirs}
#   DEPENDENCIES MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM RRG ESMF::ESMF)

# install(
#    FILES GEOS_ChemGridComp.rc ChemEnv_ExtData.rc  ChemEnv.rc ChemEnv_ExtData.yaml
#    DESTINATION etc
#    )
