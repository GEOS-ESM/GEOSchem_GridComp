esma_set_this ()

set (srcs
  Chem_RegistryMod.F90
  Runtime_RegistryMod.F90
  Chem_ArrayMod.F90
  Species_ArrayMod.F90
  Chem_BundleMod.F90
  Species_BundleMod.F90
  Chem_Mod.F90
  Chem_InitMod.F90
  Chem_MieTableMod.F90
  Chem_MieMod.F90
  Chem_AodMod.F90
  Chem_SimpleBundleMod.F90
  )
esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu MAPL GMAO_gfio_r4 esmf NetCDF::NetCDF_Fortran)
target_compile_definitions (${this} PRIVATE GEOS5)

ecbuild_add_executable(TARGET gogo.x SOURCES gogo.F90 LIBS ${this})

include_directories (${esma_include}/${this})

if (USE_F2PY)
   find_package(F2PY2)
   if (F2PY2_FOUND)
      add_f2py2_module(MieObs_ SOURCES MieObs_py.F90
         DESTINATION lib/Python2
         LIBRARIES Chem_Base MAPL GMAO_mpeu ${ESMF_LIBRARY}
         INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
         USE_MPI USE_NETCDF
         )
      add_dependencies(MieObs_ ${this})
   endif ()
endif ()


foreach (exe Chem_Aod.x Chem_Aod3d.x ctl_crst.x Chem_BundleToG5rs.x reff_calculator.xx ext_calculator.xx)
  get_filename_component (basename ${exe} NAME_WE) # base file name
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${basename}.F90
    LIBS ${this} OpenMP::OpenMP_Fortran)
endforeach ()

# Copy RC files to build tree
file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

install(PROGRAMS ext_calculator.py reff.py extinction_calculator.csh DESTINATION bin)
install(PROGRAMS mieobs.py DESTINATION lib/Python)
