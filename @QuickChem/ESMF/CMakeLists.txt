set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)
set (BUILD_GEOS_INTERFACE TRUE)

if (UFS_QUICKCHEM)
  esma_add_subdirectories(
    QAerosol_GridComp
    QuickChem_GridComp
    Shared
    UFS
  )
else ()
  esma_add_subdirectories(
    Shared
    QAerosol_GridComp
    QuickChem_GridComp
    HEMCO_GridComp
  )
endif ()
