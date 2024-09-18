esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_PChemGridComp.F90
  DEPENDENCIES MAPL MAPL_cfio_r4 Chem_Base ESMF::ESMF NetCDF::NetCDF_Fortran)
