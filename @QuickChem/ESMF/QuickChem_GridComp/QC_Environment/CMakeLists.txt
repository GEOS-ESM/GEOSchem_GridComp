esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES MAPL QChem_Shared2G QC_Process_Library
  )
