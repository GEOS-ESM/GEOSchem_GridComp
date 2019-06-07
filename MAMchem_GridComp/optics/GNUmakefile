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

#                  --------------------
#                  User Defined Targets
#                  --------------------

PY_MODS = gads.py mam7-lut.py 
PY_EXTS = optics_.$(F2PYEXT)


BINS = 
SCRP = 

DEBUG = #$(D)DEBUG

MIEV_DIR = ./miev/
MIEV_SRCS = $(MIEV_DIR)/ErrPack.f \
            $(MIEV_DIR)/MIEV0.F

SRCS = $(MIEV_SRCS) mie.F90

UT_DIR = ./tests/
UT_SRCS = $(UT_DIR)/ut_mie.F90
UT_BINS := $(notdir $(addsuffix .x, $(basename $(UT_SRCS))))

MOD_DIRS = . 

INC_DIRS = $(MOD_DIRS)

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS)))) \
        $(notdir $(addsuffix .d, $(basename $(BINS))))

# the Mie calculations require double precission
FREAL = $(FREAL8)

USER_FDEFS  = 
USER_FFLAGS = $(OMPFLAG) $(BIG_ENDIAN)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(MIEV_DIR) $(UT_DIR) $(INC_DIRS) $(MOD_DIRS)

LIBS = 
MPLIBS = 

all: python


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
	@echo "         THIS = $(THIS)"

esma_clean clean:
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *.so *.pyf *.pyc

esma_distclean distclean:
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd] *.so *.pyf *.pyc

esma_doc doc:
	@echo "No docs."


python: $(PY_EXTS)

optics_.$(F2PYEXT): $(OBJS)
	f2py --overwrite-signature -m optics_ -h mie.pyf mie.F90
	# change the variable declarations from single to double precision
	$(SED) -i -e 's/real/real(kind=8)/' -e 's/complex/complex(kind=8)/' mie.pyf
	$(F2PY) -c mie.pyf -m optics_ $(M). \
                 $(OBJS) -liomp5 \
                 only: scattering_lognormal :


tests: $(UT_BINS)


# Generic single PE apps
# ----------------------
%.x :  $(OBJS) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(OBJS) -liomp5

# Generic MPI apps
# ----------------
%.xx : $(OBJS) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(MPLIBS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

