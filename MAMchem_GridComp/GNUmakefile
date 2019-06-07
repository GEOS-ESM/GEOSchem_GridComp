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


#                  --------------------
#                  User Defined Targets
#                  --------------------

BINS = mam_optics_calculator.xx #test.x
SCRP = mam_optics_calculator.csh mam_optics_calculator.py

DEBUG = #$(D)DEBUG

MICROPHYSICS_DIR = ./microphysics/

SRCS = $(MICROPHYSICS_DIR)/infnan.F90 \
       $(MICROPHYSICS_DIR)/cam_logfile.F90 \
       $(MICROPHYSICS_DIR)/abortutils.F90 \
       $(MICROPHYSICS_DIR)/chem_mods.F90 \
       $(MICROPHYSICS_DIR)/constituents.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_data.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_newnuc.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_kind.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_aero.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_asect.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_asecthp.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_constants.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_ext.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_support.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_astem.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_main.F90 \
       $(MICROPHYSICS_DIR)/module_data_mosaic_gas.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_lsode.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_box_aerchem.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_coag.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_calcsize.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_amicphys.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_init_aerpar.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_init.F90 \
       $(MICROPHYSICS_DIR)/module_mosaic_cam_init.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_initialize_data.F90 \
       $(MICROPHYSICS_DIR)/modal_aero_wateruptake.F90 \
       MAML_SizeMod.F90 \
       MAML_SettlingMod.F90 \
       MAML_DryDepositionMod.F90 \
       MAML_DryRemovalMod.F90 \
       MAML_WetRemovalMod.F90 \
       MAML_OpticsTableMod.F90 \
       MAML_OpticsMod.F90 \
       MAM3_DataMod.F90 \
       MAM7_DataMod.F90 \
       MAM_ComponentsDataMod.F90 \
       MAM_ConstituentsDataMod.F90 \
       MAM_BaseMod.F90 \
       MAM_SizeMod.F90 \
       MAM_DryRemovalMod.F90 \
       MAM_WetRemovalMod.F90 \
       MAM_SeasaltMod.F90 \
       MAM_DustMod.F90 \
       MAM_BlackCarbonMod.F90 \
       MAM_OrganicCarbonMod.F90 \
       MAM_SulfateMod.F90 \
       MAMchem_GridCompMod.F90
#      wetdep.F90

MAMchem_GridCompMod.o: modal_aero_initialize_data.o

MOD_DIRS = . $(INC_ESMF) $(INC_GMAO_SHARED)

INC_DIRS = $(MOD_DIRS)

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS)))) \
        $(notdir $(addsuffix .d, $(basename $(BINS))))

FREAL = $(FREAL4)
THIS_CFIO = MAPL_cfio_r4

USER_FDEFS  = $(D)GEOS5 $(D)MODAL_AERO $(D)MODAL_AERO_7MODE $(D)GEOS5_PORT $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

USER_LDFLAGS = $(OMPFLAG)

ACG_FLAGS += -F

vpath % $(MICROPHYSICS_DIR) $(INC_DIRS) $(MOD_DIRS)

LIBS   = $(LIB) $(LIB_MAPL_BASE) \
         $(LIB_EU) $(LIB_CFIO) $(LIB_GFIO) \
         $(LIB_MPEU) $(LIB_MFHDF3) $(LIB_SDF) \
         $(LIB_ESMF) $(LIB_MPI) $(LIB_SYS)

MPLIBS = $(LIB) $(LIB_MAPL_BASE) $(LIB_GMAO_pFIO) \
         $(LIB_MPEU) $(LIB_CFIO) $(LIB_GFIO) \
         $(LIB_MPEU) $(LIB_MFHDF3) $(LIB_SDF) \
         $(LIB_ESMF) $(LIB_MPI) $(LIB_SYS)


$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)


$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

# Standard ESMA targets
# ---------------------
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

esma_install install: $(LIB_THIS) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p $(LIB_THIS) $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p *.rc        $(ESMAETC)
	$(CP) -p $(BINS) $(SCRP) $(ESMABIN)

esma_clean clean:
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.rc *___.h

esma_distclean distclean:
	$(RM) $(ACGS) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.rc *___.h

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


# Generic single PE apps
# ----------------------
%.x :  $(LIB_THIS) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB_THIS) $(LIBS)

# Generic MPI apps
# ----------------
%.xx : $(LIB_THIS) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB_THIS) $(MPLIBS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

