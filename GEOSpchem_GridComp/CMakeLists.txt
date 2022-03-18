esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_PChemGridComp.F90
  DEPENDENCIES MAPL MAPL_cfio_r4 Chem_Base esmf NetCDF::NetCDF_Fortran)
