esma_set_this()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/QAerosol_GridComp.F90)

    esma_add_library (${this}
            SRCS QAerosol_GridComp.F90
            DEPENDENCIES QuickChem_GridComp MAPL esmf)

    mapl_acg (${this} QAerosol_StateSpecs.rc
            IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
            GET_POINTERS DECLARE_POINTERS)

endif()
