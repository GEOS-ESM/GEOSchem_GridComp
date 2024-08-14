esma_set_this ()

set (kpp_gas_dir kpp/gas)
set (srcs
  ${kpp_gas_dir}/kpp_achem_gas_Precision.f90
  ${kpp_gas_dir}/kpp_achem_gas_Parameters.f90
  ${kpp_gas_dir}/kpp_achem_gas_Global.f90
  ${kpp_gas_dir}/kpp_achem_gas_Function.f90
  ${kpp_gas_dir}/kpp_achem_gas_JacobianSP.f90
  ${kpp_gas_dir}/kpp_achem_gas_Jacobian.f90
  ${kpp_gas_dir}/kpp_achem_gas_HessianSP.f90
  ${kpp_gas_dir}/kpp_achem_gas_Hessian.f90
  ${kpp_gas_dir}/kpp_achem_gas_StoichiomSP.f90
  ${kpp_gas_dir}/kpp_achem_gas_Stoichiom.f90
  ${kpp_gas_dir}/kpp_achem_gas_Rates.f90
  ${kpp_gas_dir}/kpp_achem_gas_Monitor.f90
  ${kpp_gas_dir}/kpp_achem_gas_Util.f90
  ${kpp_gas_dir}/kpp_achem_gas_LinearAlgebra.f90
  ${kpp_gas_dir}/kpp_achem_gas_Integrator.f90
  GACL_ConstantsMod.F90
  GACL_EmissionsMod.F90
  GACL_ReactionRatesMod.F90
  GACL_DryDepositionMod.F90
  GEOS_AChemGridCompMod.F90
  )

esma_add_library (
  ${this}
  SRCS ${srcs}
  DEPENDENCIES Chem_Shared MAPL GMAO_mpeu ESMF::ESMF NetCDF::NetCDF_Fortran
  )
target_compile_definitions (${this} PRIVATE MAPL_MODE GEOS5)
set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${PP})

new_esma_generate_automatic_code (
  ${this} GEOSachem_Registry.rc
  "GEOS_AChem_ExportSpec___.h;GEOS_AChem_GetPointer___.h"
  GEOS_AChem_History___.rc
  ${include_GEOSachem_GridComp} ${esma_etc}
  "-f"
  )

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

