esma_set_this ()

esma_mepo_style(Cloud-J cloud_j_dir)

set (MODEL_GEOSCHEM TRUE) # Will need to rework to use with GMI instead
set (MAPL_ESMF TRUE)

set (srcs
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_cmn_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_error_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_fjx_osa_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_fjx_sub_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_init_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${cloud_j_dir}/src/Core/cldj_sub_mod.F90
   )

# placeholder for excluding files from build
#list (REMOVE_ITEM srcs
#   ${CMAKE_CURRENT_SOURCE_DIR}${cloud_j_dir}/src/Core/placeholder.F90
#   )

# MAPL, ESMF, and netCDF not used but include here for potential future use
esma_add_library (${this}
   SRCS ${srcs}
   DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran
   )

target_compile_definitions(${this} PRIVATE MAPL_ESMF MODEL_GEOSCHEM)
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})
