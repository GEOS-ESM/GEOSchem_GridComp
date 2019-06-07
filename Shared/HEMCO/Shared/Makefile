#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/ext subdirectory)
#
# !DESCRIPTION: This Makefile builds external routines used by standalone
#  HEMCO into library file libEXT.a.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#
#   make -jN TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
#  27 Jul 2014 - C. Keller   - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT :=..
BIN  :=$(ROOT)/bin
DOC  :=$(ROOT)/doc
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod
RUN  :=$(ROOT)/run

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source code files to compile
# (i.e. everything ending in either .F or .F90)
SRC  := $(wildcard *.F) $(wildcard *.F90)

# Replace .F and .F90 extensions with *.o
TMP  := $(SRC:.F=.o)
OBJ  := $(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug

all: lib

lib: $(OBJ)
	$(AR) crs libEXT.a $(OBJ)
	mv libEXT.a $(LIB)

clean:
	rm -f *.o *.mod *.a *.x

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "BIN     : $(BIN)"
	@echo "DOC     : $(DOC)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "RUN     : $(RUN)"
	@echo "F90     : $(F90)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

charpak_mod.o          : charpak_mod.F

henry_mod.o            : henry_mod.F \
                         precision_mod.o

inquireMod.o           : inquireMod.F90

julday_mod.o           : julday_mod.F \
                         precision_mod.o

precision_mod.o        : precision_mod.F

regrid_a2a_mod.o       : regrid_a2a_mod.F90 \
                         precision_mod.o

#EOC
