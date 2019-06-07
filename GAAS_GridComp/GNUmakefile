#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS := $(shell basename `pwd`)
NAME := GAAS
LIB   = lib$(THIS).a

BINS = ana_lde.py ana_lde.x

esma_install install: $(DEPS) $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.rc        $(ESMAETC)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p $(BINS)     $(ESMABIN)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd] *.x

esma_distclean distclean:
	$(RM) $(ACGS) *~ *.[aoxd] *.[Mm][Oo][Dd] *.x

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


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
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS  = $(wildcard m_ana.F90 LDE_Mod.F90 GAAS_GridCompMod.F90)

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS))) 
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h \
	$(NAME)_ImportSpec___.h $(NAME)_DeclarePointer___.h \
        $(NAME)_History___.rc

ACG_FLAGS += -F

#FOPT = $(FOPT3)

MOD_DIRS = . $(INC_OBS_AOD) $(INC_ODS) $(INC_PSAS) $(INC_HERMES) \
             $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_MPI) $(INC_GEOS_FV3)

INC_DIRS = $(MOD_DIRS)

THIS_GFIO := GMAO_gfio_r4
THIS_CFIO := MAPL_cfio_r4

USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(MOD_DIRS),$(I)$(dir)) 

USER_LDFLAGS = $(OMPFLAG)

FREAL = $(FREAL4)

vpath % $(MOD_DIRS)

$(LIB) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

#LIB_SCI_ = -llapack -lblas

LIB_FVCUBE = $(ESMALIB)/libFVdycoreCubed_GridComp.a\
             $(ESMALIB)/libfvdycore.a \
             $(ESMALIB)/libGFDL_fms.a 

%.x : $(LIB) %.o
	$(FC) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_FVCUBE) \
              $(LIB_CHEM_BASE) $(LIB_HERMES) \
              $(LIB_MAPL_BASE) $(LIB_GMAO_pFIO) $(LIB_FVCUBE) $(LIB_GFIO) \
              $(LIB_CFIO) $(LIB_MPEU) \
              $(LIB_SDF) $(LIB_ESMF) $(LIB_SCI) $(LIB_SYS) $(LIB_MPI)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
