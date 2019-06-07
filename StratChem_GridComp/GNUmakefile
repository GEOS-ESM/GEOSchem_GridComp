#
# Makefile for ESMA components.
#
#
# REVISION HISTORY:
#
# 17Oct2008  da Silva  Changed from distributed to centralized makefile.
#
#-------------------------------------------------------------------------

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO specific macros

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS := $(shell basename `pwd`)
#CMP := $(patsubst %_GridComp,%,$(THIS))
#NAME := $(shell echo $(CMP) | tr a-z A-Z)
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
	$(CP) -p SC_GridComp/*.rc        $(ESMAETC)
	$(CP) -p *.h *.mod   $(INC_CHEM)

esma_clean clean:
	$(RM) $(ACGS) $(OBJS) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) $(ACGS) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex

#                  --------------------
#                  User Defined Targets
#                  --------------------

ifeq ($(DOING_GEOS5),TRUE)
  DGEOS5 = $(D)GEOS5
endif
DEBUG = -DDEBUG

SRC_DIRS := . SC_GridComp

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_MPI)

MOD_DIRS = . $(INC_DIRS)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h \
	$(NAME)_ImportSpec___.h \
	$(NAME)_GetPointer___.h \
	$(NAME)_History___.rc

DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

# -------------------------------------------------------------------
# To run the reduced SC equation set add the following to USER_FDEFS:
#  -DREDUCED
# See usage notes at the head of SC_GridCompMod.F90.
# -------------------------------------------------------------------

USER_FDEFS  = $(DGEOS5) $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)


$(LIB_THIS) lib : $(ACGS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)

$(ACGS) : SC_GridComp/SC_Registry.rc $(ACG)
	@$(ACG) -N STRATCHEM SC_GridComp/SC_Registry.rc

#$(DEPS) : $(ACGS)

SC_GridCompMod.o : SC_GridComp/SC_GridCompMod.F90
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
