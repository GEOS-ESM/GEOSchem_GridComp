#
# Makefile for ESMA components.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
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
	@echo "         SITE = $(SITE) "


THIS := $(shell basename `pwd`)
LIB   = $(LIB_THIS)

#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ALLDIRS = GEOSpchem_GridComp GOCART_GridComp \
          StratChem_GridComp GMIchem_GridComp \
          CARMAchem_GridComp GEOSCHEMchem_GridComp \
          MATRIXchem_GridComp MAMchem_GridComp \
          GAAS_GridComp H2O_GridComp TR_GridComp \
          GEOSachem_GridComp DNA_GridComp HEMCO_GridComp

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

export ESMADIR BASEDIR ARCH SITE FREAL

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
	$(CP) -p *.rc           $(ESMAETC)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.[fF]*

local_esma_doc local_doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex
	#@echo "Target $@ not implemented yet in `pwd`"

#                  --------------------
#                  User Defined Targets
#                  --------------------

# If children not present, create stubs
# -------------------------------------
ifeq ( $(wildcard GMIchem_GridComp), $(null) )
        STUBS += GMIchem_GridCompMod___.F90
endif
ifeq ( $(wildcard GOCART_GridComp), $(null) )
        STUBS += GOCART_GridCompMod___.F90
endif
ifeq ( $(wildcard GAAS_GridComp), $(null) )
        STUBS += GAAS_GridCompMod___.F90
endif
ifeq ( $(wildcard H2O_GridComp), $(null) )
        STUBS += H2O_GridCompMod___.F90
endif
ifeq ( $(wildcard StratChem_GridComp), $(null) )
        STUBS += StratChem_GridCompMod___.F90
endif
ifeq ( $(wildcard CARMAchem_GridComp), $(null) )
        STUBS += CARMAchem_GridCompMod___.F90
endif
ifeq ( $(wildcard GEOSCHEMchem_GridComp), $(null) )
        STUBS += GEOSCHEMchem_GridCompMod___.F90
endif
ifeq ( $(wildcard MATRIXchem_GridComp), $(null) )
        STUBS += MATRIXchem_GridCompMod___.F90
endif
ifeq ( $(wildcard MAMchem_GridComp), $(null) )
        STUBS += MAMchem_GridCompMod___.F90
endif
ifeq ( $(wildcard GEOSachem_GridComp), $(null) )
        STUBS += GEOS_AChemGridCompMod___.F90
endif
ifeq ( $(wildcard TR_GridComp), $(null) )
        STUBS += TR_GridCompMod___.F90
endif
ifeq ( $(wildcard DNA_GridComp), $(null) )
        STUBS += DNA_GridCompMod___.F90
endif
ifeq ( $(wildcard HEMCO_GridComp), $(null) )
        STUBS += HEMCO_GridCompMod___.F90
endif

SRCS := GEOS_ChemEnvGridComp.F90 GEOS_ChemGridComp.F90 $(STUBS)
OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

THIS_CFIO = MAPL_cfio_r4

INC_DIRS = . $(INC_ESMF) $(INC_GMAO_SHARED) $(INC_GEOS_CHEM)

MOD_DIRS = . $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

# Stubs
# -----
%___.F90: $(STUB) 
	$(STUB) $* > $@ 

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------
  TR_GridComp_install : GMIchem_GridComp_install

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

