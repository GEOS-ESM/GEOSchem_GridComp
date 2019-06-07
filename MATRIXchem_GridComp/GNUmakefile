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
       ESMADIR := $(PWD)/../../..
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

esma_install install: $(LIB_THIS) 
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

SRC_DIRS := . microphysics

MOD_DIRS = . $(INC_ESMF) $(INC_MPEU) $(INC_MPI) $(INC_GMAO_SHARED) \
             $(INC_MAPL_BASE) $(INC_MFHDF3) $(INC_GEOS_SHARED) $(INC_MAPL_BASE) \
             $(INC_CHEM_BASE) $(INC_CHEM_SHARED) 

INC_DIRS = $(MOD_DIRS)

THIS_GFIO := GMAO_gfio_r4
THIS_CFIO := MAPL_cfio_r4

SRCS := CONST.F \
        TRAMP_param.F \
        TRAMP_config.F90 \
        TRAMP_setup.F \
        TRAMP_quad.F \
        TRAMP_coag.F \
        TRAMP_npf.F \
        TRAMP_actv.F \
        TRAMP_diam.F \
        TRAMP_drv.F \
        TRAMP_subs.F \
        TRAMP_depv.F \
        TRAMP_nomicrophysics.F \
        TRAMP_isofwd2.F \
        TRAMP_isorev2.F \
        TRAMP_isocom2.F \
        TRAMP_thermo_isorr2.F \
        TRAMP_matrix.F \
        MATRIXchem_GridCompMod.F90

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

USER_FDEFS  = $DTRACERS_AMP $(D)TRACERS_AMP_M1 $(D)GEOS5_PORT $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) $(PP)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

ACG_FLAGS += -F

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)


$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

