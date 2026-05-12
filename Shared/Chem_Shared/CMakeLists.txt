esma_set_this()

set (srcs
  Chem_ConstMod.F90
  Chem_SettlingMod.F90
  DryDepositionMod.F90
  DustEmissionMod.F90
  WetRemovalMod.F90
  SeasaltEmissionMod.F90
  NitrateChemDriverMod.F90
  SulfateChemDriverMod.F90
  ConvectionMod.F90
  SmokeEmissionMod.F90
  VegLaiMod.F90
  GmiPrintError_mod.F90
  GmiTimeControl_mod.F90
  GmiSolar_mod.F90
  GmiResistance.F90
  GmiEmissUtils.F90
  Henrys_Law_cts.F90
  Chem_GroupMod.F90
  Bundle_IncrementMod.F90
  Lightning_mod.F90
  lightning_toolbox_mod.F90
  SZA_from_MAPL_mod.F90
  mod_diag.F90
  diaglist.F90
  Chem_UtilMod.F90
  Chem_StateMod.F90
  )

esma_add_library(${this}
   SRCS ${srcs}
   DEPENDENCIES GMAO_mpeu MAPL Chem_Base GEOS_Shared GMAO_hermes ESMF::ESMF NetCDF::NetCDF_Fortran)
#add_definitions(-DDEBUG)
add_definitions(-DSPMD -DFVCHEM)

file (COPY gmi_drydep_data.h DESTINATION ${include_${this}})
file (COPY gmi_emiss_constants.h DESTINATION ${include_${this}})
file (COPY gmi_phys_constants.h DESTINATION ${include_${this}})
file (COPY gmi_time_constants.h DESTINATION ${include_${this}})

install(
  PROGRAMS bundleParser.py
  DESTINATION bin
)


