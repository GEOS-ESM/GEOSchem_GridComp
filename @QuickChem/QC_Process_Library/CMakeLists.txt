esma_set_this ()

set (srcs
  QuickChem_Process.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES esmf NetCDF::NetCDF_Fortran
  )
