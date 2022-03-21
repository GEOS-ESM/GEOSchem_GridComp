esma_set_this ()

set (srcs
  microphysics/CONST.F
  microphysics/TRAMP_param.F
  microphysics/TRAMP_config.F90
  microphysics/TRAMP_setup.F
  microphysics/TRAMP_quad.F
  microphysics/TRAMP_coag.F
  microphysics/TRAMP_npf.F
  microphysics/TRAMP_actv.F
  microphysics/TRAMP_diam.F
  microphysics/TRAMP_drv.F
  microphysics/TRAMP_subs.F
  microphysics/TRAMP_depv.F
  microphysics/TRAMP_nomicrophysics.F
  microphysics/TRAMP_isofwd2.F
  microphysics/TRAMP_isorev2.F
  microphysics/TRAMP_isocom2.F
  microphysics/TRAMP_thermo_isorr2.F
  microphysics/TRAMP_matrix.F
  MATRIXchem_GridCompMod.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES Chem_Shared MAPL esmf)
if (EXTENDED_SOURCE)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()
target_compile_definitions (${this} PRIVATE TRACERS_AMP TRACERS_AMP_M1 GEOS5_PORT)

esma_generate_gocart_code (${this} -F)

file (GLOB resource_files CONFIGURE_DEPENDS "*.rc")

install(
   FILES ${resource_files}
   DESTINATION etc
   )
