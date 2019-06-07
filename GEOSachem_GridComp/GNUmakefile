#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
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
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.rc *___.h

esma_distclean distclean:
	$(RM) $(ACGS) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.rc *___.h

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

DEBUG = #$(D)DEBUG

KPP_GAS_DIR = kpp/gas/
KPP_GAS_SRC = $(KPP_GAS_DIR)/kpp_achem_gas_Precision.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Parameters.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Global.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Function.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_JacobianSP.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Jacobian.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_HessianSP.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Hessian.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_StoichiomSP.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Stoichiom.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Rates.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Monitor.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Util.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_LinearAlgebra.f90 \
              $(KPP_GAS_DIR)/kpp_achem_gas_Integrator.f90

KPP_SRCS = $(KPP_GAS_SRC) 

GACL_SRCS = GACL_ConstantsMod.F90 \
            GACL_EmissionsMod.F90 \
            GACL_ReactionRatesMod.F90 \
            GACL_DryDepositionMod.F90

SRCS = $(KPP_SRCS) $(GACL_SRCS) GEOS_AChemGridCompMod.F90


MOD_DIRS = . $(INC_ESMF) $(INC_GMAO_SHARED)

INC_DIRS = $(MOD_DIRS)

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

USER_FDEFS  = $(D)MAPL $(D)GEOS5 $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

USER_FFLAGS += $(PP)

ACG_FLAGS += -F

vpath % $(KPP_GAS_DIR) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)


$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

GEOS_AChemGridCompMod.o: kpp_achem_gas_Integrator.o

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

