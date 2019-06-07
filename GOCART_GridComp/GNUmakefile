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
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO specific macros

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

# Handle stubs for NCEP implementarion of GOCART Stand alone
# ----------------------------------------------------------
  ifeq ($(GOCART_MODE),stub)
     ALLDIRS = 
     SRCS := GOCART_StubCompMod.F90
  else
     ALLDIRS = O3_GridComp CO_GridComp CO2_GridComp CFC_GridComp \
               DU_GridComp SS_GridComp  BC_GridComp  OC_GridComp \
               SU_GridComp Rn_GridComp CH4_GridComp  NI_GridComp \
               BRC_GridComp
     SRCS := GOCART_GridCompMod.F90 Aero_GridCompMod.F90  
  endif

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

DEBUG = # $(D)DEBUG

ifeq ($(ESMA_REAL),$(FREAL8))
	FREAL = $(FREAL8)
	THIS_CFIO = MAPL_cfio_r8
else
	FREAL = $(FREAL4)
	THIS_CFIO = MAPL_cfio_r4
endif

export ESMADIR BASEDIR ARCH SITE FREAL DEBUG

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) -p $(INC_CHEM)
	$(CP) -p *.mod          $(INC_CHEM)
	$(CP) -p GOCARTdata_AerRegistry.rc $(ESMAETC)
	$(CP) -p GOCARTdata_ExtData.rc     $(ESMAETC)

local_esma_clean local_clean:
	-$(RM) -f *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean: clean_rc
	-$(RM) -f *~ *.[aoxd] *.[Mm][Oo][Dd] *___.[fF]* 

local_esma_doc local_doc:
	#@echo "Target $@ not implemented yet in `pwd`"
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex

#                  --------------------
#                  User Defined Targets
#                  --------------------

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

BIG_ENDIAN =
THIS_GFIO = GMAO_gfio_r4
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

INC_DIRS = . $(INC_ESMF) $(INC_GMAO_SHARED) $(INC_CHEM_GOCART)

MOD_DIRS = . $(INC_DIRS) 

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) $(DEBUG)

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

#                      ---------
#                      Test Code
#                      ---------

LIB_CHILDREN = $(wildcard \
			  $(ESMALIB)/libCFC_GridComp.a\
			  $(ESMALIB)/libRn_GridComp.a\
			  $(ESMALIB)/libBC_GridComp.a\
                          $(ESMALIB)/libCH4_GridComp.a\
                          $(ESMALIB)/libCO2_GridComp.a\
                          $(ESMALIB)/libCO_GridComp.a\
                          $(ESMALIB)/libDU_GridComp.a\
                          $(ESMALIB)/libO3_GridComp.a\
                          $(ESMALIB)/libOC_GridComp.a\
                          $(ESMALIB)/libBRC_GridComp.a\
                          $(ESMALIB)/libNI_GridComp.a\
                          $(ESMALIB)/libSS_GridComp.a\
                          $(ESMALIB)/libSU_GridComp.a)
%.x : $(LIB) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_CHILDREN) \
	      $(LIB_CHEM_BASE) $(LIB_CHEM_SHARED) $(LIB_PILGRIM)\
              $(LIB_MAPL_BASE) $(LIB_CFIO) $(LIB_GFIO) $(LIB_MPEU) \
              $(LIB_ESMF) $(LIB_SDF) \
              $(LIB_SYS) $(LIB_MPI)

test_rc: $(LIB)
	$(CP) $(ESMAETC)/PIESA/Chem_Registry.rc       ./Tests
	$(CP) $(ESMAETC)/PIESA/Aod_Registry.rc  ./Tests
	$(CP) $(ESMAETC)/PIESA/??_GridComp.rc         ./Tests
	$(CP) $(ESMAETC)/PIESA/???_GridComp.rc        ./Tests

clean_rc:
	-$(RM) -f Tests/ExtData Tests/*_Registry.rc Tests/*_GridComp.rc 

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

