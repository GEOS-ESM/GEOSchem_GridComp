#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

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
NAME = GEOSpchem
THIS = $(NAME)_GridComp
LIB  = lib$(THIS).a
esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS = GEOS_PChemGridComp.F90
OBJS = $(SRCS:.F90=.o)

INC_DIRS = . $(INC_GMAO_SHARED)
MOD_DIRS = . $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))
USER_FFLAGS = $(BIG_ENDIAN) $(USER_FMODS)

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)


esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.F90 > $(ESMADOC)/$(THIS).tex

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

