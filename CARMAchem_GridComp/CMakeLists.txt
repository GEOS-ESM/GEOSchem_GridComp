esma_set_this ()

set (src_directories
  .
  CARMA
  CARMA/source/base
  )

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()

set (dependencies Chem_Shared Chem_Base GMAO_mpeu)
esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
target_include_directories (${this} PUBLIC CARMA/source/base)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

new_esma_generate_automatic_code (
  ${this} CARMAchem_Registry.rc
  "CARMA_ExportSpec___.h;CARMA_GetPointer___.h"
  CARMA_History___.rc
  ${include_GEOSchem_GridComp} ${esma_etc}
  -v
  )
