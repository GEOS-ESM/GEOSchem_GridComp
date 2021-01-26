esma_set_this ()

set (microphysics_dir microphysics)
set (srcs
  ${microphysics_dir}/infnan.F90
  ${microphysics_dir}/cam_logfile.F90
  ${microphysics_dir}/abortutils.F90
  ${microphysics_dir}/chem_mods.F90
  ${microphysics_dir}/constituents.F90
  ${microphysics_dir}/modal_aero_data.F90
  ${microphysics_dir}/modal_aero_newnuc.F90
  ${microphysics_dir}/module_data_mosaic_kind.F90
  ${microphysics_dir}/module_data_mosaic_aero.F90
  ${microphysics_dir}/module_data_mosaic_asect.F90
  ${microphysics_dir}/module_data_mosaic_asecthp.F90
  ${microphysics_dir}/module_data_mosaic_constants.F90
  ${microphysics_dir}/module_mosaic_ext.F90
  ${microphysics_dir}/module_mosaic_support.F90
  ${microphysics_dir}/module_mosaic_astem.F90
  ${microphysics_dir}/module_data_mosaic_main.F90
  ${microphysics_dir}/module_data_mosaic_gas.F90
  ${microphysics_dir}/module_mosaic_lsode.F90
  ${microphysics_dir}/module_mosaic_box_aerchem.F90
  ${microphysics_dir}/modal_aero_coag.F90
  ${microphysics_dir}/modal_aero_calcsize.F90
  ${microphysics_dir}/modal_aero_amicphys.F90
  ${microphysics_dir}/module_mosaic_init_aerpar.F90
  ${microphysics_dir}/module_mosaic_init.F90
  ${microphysics_dir}/module_mosaic_cam_init.F90
  ${microphysics_dir}/modal_aero_initialize_data.F90
  ${microphysics_dir}/modal_aero_wateruptake.F90
  MAML_SizeMod.F90
  MAML_SettlingMod.F90
  MAML_DryDepositionMod.F90
  MAML_DryRemovalMod.F90
  MAML_WetRemovalMod.F90
  MAML_OpticsTableMod.F90
  MAML_OpticsMod.F90
  MAM3_DataMod.F90
  MAM7_DataMod.F90
  MAM_ComponentsDataMod.F90
  MAM_ConstituentsDataMod.F90
  MAM_BaseMod.F90
  MAM_SizeMod.F90
  MAM_DryRemovalMod.F90
  MAM_WetRemovalMod.F90
  MAM_SeasaltMod.F90
  MAM_DustMod.F90
  MAM_BlackCarbonMod.F90
  MAM_OrganicCarbonMod.F90
  MAM_SulfateMod.F90
  MAMchem_GridCompMod.F90
  #      wetdep.F90
  )


set (dependencies MAPL Chem_Shared GMAO_mpeu)
esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
target_compile_definitions (${this} PRIVATE GEOS5 MODAL_AERO MODAL_AERO_7MODE GEOS5_PORT)

esma_generate_gocart_code (${this} -F)

file (GLOB resource_files CONFIGURE_DEPENDS "*.rc")

set (mamscripts mam_optics_calculator.py mam_optics_calculator.csh)

ecbuild_add_executable(TARGET mam_optics_calculator.xx SOURCES mam_optics_calculator.F90 LIBS ${this})

install(
   FILES ${resource_files}
   DESTINATION etc
   )

install(
   PROGRAMS ${mamscripts}
   DESTINATION bin)

