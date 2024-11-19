esma_set_this ()

esma_mepo_style(HETP hetp_dir)

set (srcs
   ${CMAKE_CURRENT_SOURCE_DIR}/${hetp_dir}/src/Core/hetp_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/${hetp_dir}/src/Core/mach_hetp_mod.F90
   )

esma_add_library (${this}
   SRCS ${srcs}
   )
