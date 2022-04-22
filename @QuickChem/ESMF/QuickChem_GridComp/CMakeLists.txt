esma_set_this ()

set (alldirs
  QC_Environment
  OH_GridComp
  XX_GridComp
  )

set (srcs
  QuickChem_GridCompMod.F90
  )

set (resource_files
   QuickChem_GridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (dependencies MAPL QChem_Shared2G QC_Process_Library esmf)
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES ${dependencies}
  )

mapl_acg (${this}   QuickChem_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
