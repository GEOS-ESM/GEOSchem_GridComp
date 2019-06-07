#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 17Oct2008  da Silva  Changed from distributed to centralized makefile.
#
#-------------------------------------------------------------------------

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS := $(shell basename `pwd`)
NAME := $(patsubst %_GridComp,%,$(THIS))
LIB_THIS := lib$(THIS).a

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

esma_install install: $(ACGS) $(LIB_THIS) 
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p $(LIB_THIS) $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p *.rc        $(ESMAETC)

esma_clean clean:
	$(RM) $(OBJS) *~ *.[aox] *.[Mm][Oo][Dd] *___.rc *___.h

esma_distclean distclean:
	$(RM) $(ACGS) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.rc *___.h

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRC_DIRS := .                            \
            ./gc_column/Headers          \
            ./gc_column/History          \
            ./gc_column/GeosUtil         \
            ./gc_column/NcdfUtil         \
            ./gc_column/ISOROPIA         \
            ./gc_column/KPP/Standard     \
            ./gc_column/GeosCore         \
            ./gc_column/GIGC

INC_DIRS := . ./gc_column/Headers $(INC_MAPL_BASE) \
             $(INC_GMAO_SHARED) $(INC_MFHDF3) 

MOD_DIRS = . $(INC_ESMF) $(INC_MPEU) $(INC_MPI) \
             $(INC_DIRS) $(ESMAINC)/MAPL_Base \
             $(INC_GEOS_SHARED)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

SRCS := $(filter-out ./gc_column/GeosCore/main.F 	      \
		     ./gc_column/GeosCore/tpcore_fvdas_mod.F  \
		     ./gc_column/GeosCore/tpcore_geos57_window_mod.F \
		     ./gc_column/GeosCore/tpcore_geos5_window_mod.F  \
		     ./gc_column/GeosCore/transport_mod.F     \
		     ./gc_column/GeosCore/rrtmg_rad_transfer_mod.F \
		     ./gc_column/GIGC/GIGC_GridCompMod.F90    \
		     ./gc_column/GIGC/Chem_GridCompMod.F90    \
		     ./gc_column/GIGC/GEOS_ctmEnvGridComp.F90    \
		     ./gc_column/GIGC/GEOSChem.F90    \
		     ./gc_column/Headers/inquireMod.F90    \
                     , $(SRCS))

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h \
        $(NAME)_DeclarePointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

USER_FDEFS  = $(D)ESMF_ $(D)DEVEL $(D)NO_BPCH $(D)UCX $(D)EXTERNAL_GRID $(D)EXTERNAL_FORCING $(D)GEOS_FP $(D)USE_REAL8 $(D)USE_TEND $(D)NC_DIAG $(D)DISCOVER $(D)MODEL_GEOS
USER_CDEFS  = $(I)ESMF_ $(I)DEVEL
USER_FFLAGS = $(D)ESMF_ $(BIG_ENDIAN)

ifeq ($(ESMA_FC),ifort)
   USER_FDEFS += $(D)LINUX_IFORT
else
ifeq ($(ESMA_FC),pgfortran)
   USER_FDEFS += $(D)LINUX_PGI
else
ifeq ($(ESMA_FC),gfortran)
   USER_FDEFS += $(D)LINUX_GFORTRAN
endif
endif
endif

USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

ACG_FLAGS += -F

vpath % /usr/include $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 


$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

%.x : $(LIB) %.o
	@echo "$(LIB_MPI)"
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB_THIS) \
	$(LIB_CHEM_BASE) $(LIB_CHEM_SHARED) $(LIB_PILGRIM)\
	$(LIB_MAPL_BASE) $(LIB_CFIO) $(LIB_MPEU) \
	$(LIB_ESMF) $(LIB_SDF) \
	$(LIB_SYS) $(LIB_MPI)

test_rc: $(LIB)
	$(CP) $(ESMAETC)/GEOSCHEMchem_GridComp.rc	  ./Tests

clean_rc:
	-$(RM) -f Tests/ExtData Tests/*_Registry.rc Tests/*_GridComp.rc 

# Reduce the optimization for (only)
# GEOSCHEMchem_GridCompMod.F90 if it compiles too slowly.
# -------------------------------------------------------
GEOSCHEMchem_GridCompMod.o : GEOSCHEMchem_GridCompMod.F90
	$(FC) -c $(patsubst $(FOPT),$(FOPT2),$(F90FLAGS)) $<

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
