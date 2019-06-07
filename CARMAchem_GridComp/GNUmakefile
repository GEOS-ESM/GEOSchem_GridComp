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
	$(RM) $(OBJS) *~ *.[aox] *.[Mm][Oo][Dd] *___* *.d

esma_distclean distclean:
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo No documentation here
#	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

ifeq ($(DOING_GEOS5),TRUE)
  DGEOS5 = $(D)GEOS5
endif
#DEBUG = -DDEBUG

SRC_DIRS := . \
            CARMA \
            CARMA/source/base

INC_DIRS := . $(INC_ESMF) $(INC_GMAO_SHARED) \
            CARMA/source/base 
 
MOD_DIRS = . $(INC_DIRS)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

BIG_ENDIAN =
THIS_GFIO = GMAO_gfio_r4
FREAL = $(FREAL4)
#FOPT = $(FOPT3)

ifeq ($(ESMA_REAL),$(FREAL8))
        THIS_GFIO = GMAO_gfio_r8
        THIS_CFIO = MAPL_cfio_r8
        FREAL = $(FREAL8)
else
        THIS_GFIO = GMAO_gfio_r4
        THIS_CFIO = MAPL_cfio_r4
        FREAL = $(FREAL4)
endif

# MAT Intel 18 + init snan has a ICE with miess. Remove for now
ifeq ($(ESMA_FC),ifort)
   ifeq ("$(BOPT)","g")
      miess.o : FOPT = $(FOPTG) -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -warn unused 
   endif
endif


USER_FDEFS  = $(DGEOS5) $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)

$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

%.x : $(LIB_THIS) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB_THIS) \
	      $(LIB_CHEM_BASE) $(LIB_CHEM_SHARED) $(LIB_PILGRIM)\
              $(LIB_MAPL_BASE) \
              $(LIB_CFIO) $(LIB_GFIO) $(LIB_MPEU) \
              $(LIB_ESMF) $(LIB_SDF) \
              $(LIB_SYS) $(LIB_MPI)

test_rc: $(LIB)
	$(CP) $(ESMAETC)/CARMAchem_Registry.rc       ./Tests

clean_rc:
	-$(RM) -f Tests/ExtData Tests/*_Registry.rc 

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

