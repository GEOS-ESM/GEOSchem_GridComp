#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 17Oct2008  da Silva  Changed from distributed to centralized makefile.
# 12Apr2011  Nielsen   Dependencies for automatic code generation.
#  6Mar2012  Nielsen   Ganymed-1_0_UNSTABLE. ESMF bundle for tendency exports.
#
#-------------------------------------------------------------------------

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../../../..
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
CMP := $(patsubst %_GridComp,%,$(THIS))
NAME := $(shell echo $(CMP) | tr a-z A-Z)
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
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMAETC)/CCMI_REF-C1 $(ESMAETC)/CCMI_REF-C2
	$(MKDIR) $(ESMAETC)/CCMI_REF-C1/chemistry_files $(ESMAETC)/CCMI_REF-C2/chemistry_files
	$(CP) -p $(LIB_THIS) $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	$(CP) -p *.rc                           $(ESMAETC)
	$(CP) -p GMI_GridComp/*.rc              $(ESMAETC)
	$(CP) -p GMI_GridComp/CCMI_REF-C1/*.rc  $(ESMAETC)/CCMI_REF-C1/chemistry_files
	$(CP) -p GMI_GridComp/CCMI_REF-C2/*.rc  $(ESMAETC)/CCMI_REF-C2/chemistry_files


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

SRC_DIRS := . \
            GMI_GridComp \
            GMI_GridComp/GmiChemistry\
            GMI_GridComp/GmiChemistry/AerosolDust\
            GMI_GridComp/GmiChemistry/ioChemistry\
            GMI_GridComp/GmiChemistry/photolysis/fastJX53c_ref\
            GMI_GridComp/GmiChemistry/photolysis/fast_JX\
            GMI_GridComp/GmiChemistry/photolysis/fast_JX65 \
            GMI_GridComp/GmiChemistry/photolysis/fastj\
            GMI_GridComp/GmiChemistry/photolysis/fast_JX53b\
            GMI_GridComp/GmiChemistry/photolysis/fast_JX53c\
            GMI_GridComp/GmiChemistry/photolysis/lookup\
            GMI_GridComp/GmiChemistry/photolysis/utils\
            GMI_GridComp/GmiChemistry/sad\
            GMI_GridComp/GmiChemistry/smv2chem\
            GMI_GridComp/GmiChemistry/stratTropMech\
            GMI_GridComp/GmiChemistry/sulfur\
            GMI_GridComp/GmiDeposition\
            GMI_GridComp/GmiEmission\
            GMI_GridComp/GmiEmission/GOCARTroutines\
            GMI_GridComp/GmiEmission/gsfc\
            GMI_GridComp/GmiEmission/harvard\
            GMI_GridComp/GmiEmission/ioEmission\
            GMI_GridComp/GmiEmission/lightning\
            GMI_GridComp/GmiEmission/llnl\
            GMI_GridComp/GmiShared/GmiESMF\
            GMI_GridComp/GmiShared/GmiIOutilities\
            GMI_GridComp/GmiShared/GmiSupportingModules\
            GMI_GridComp/GmiSpeciesConcentration/spcConcentrationMethod

INC_DIRS := . $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_GEOS_SDYN) \
            GMI_GridComp/GmiChemistry/include\
            GMI_GridComp/GmiChemistry/photolysis/include\
            GMI_GridComp/GmiChemistry/stratTropMech\
            GMI_GridComp/GmiShared/GmiInclude \
            GMI_GridComp/GmiChemistry/photolysis/fast_JX \
            GMI_GridComp/GmiChemistry/photolysis/fast_JX65 \
            GMI_GridComp/GmiChemistry/photolysis/fast_JX53b \
            GMI_GridComp/GmiChemistry/photolysis/fastJX53c_ref \
            GMI_GridComp/GmiChemistry/photolysis/fastj 
 
MOD_DIRS = . $(INC_DIRS)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h  \
	$(NAME)_ImportSpec___.h \
	$(NAME)_InternalSpec___.h \
	$(NAME)_History___.rc \
	Deposition_ExportSpec___.h \
	Deposition_Registry___.rc \
	Reactions_ExportSpec___.h \
	Reactions_Registry___.rc \
	Tendency_ExportSpec___.h \
	Tendency_Registry___.rc \
	GMI_GridComp/SCAV_FillExports___.h \
	GMI_GridComp/DD_FillExports___.h \
	GMI_GridComp/WD_FillExports___.h \
	GMI_GridComp/QK_FillExports___.h \
	GMI_GridComp/QJ_FillExports___.h \
	GMI_GridComp/QQK_FillExports___.h \
	GMI_GridComp/QQJ_FillExports___.h \
	GMI_GridComp/Deposition_DeclarePointer___.h \
	GMI_GridComp/Deposition_GetPointer___.h \
	GMI_GridComp/Reactions_DeclarePointer___.h \
	GMI_GridComp/Reactions_GetPointer___.h

DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

USER_FDEFS  = $(DGEOS5) $(DEBUG)
USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : GMI_GridComp/SCAV_FillExports___.h $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_THIS)

# 3-step automatic code generation.
# ---------------------------------
Tendency_ExportSpec___.h : GMI_GridComp/GMI_Registry.rc $(AGC)
	./gmi_acg.pl -R
	@echo " "
	$(ACG) -N $(NAME) Deposition_Registry___.rc
	/bin/mv -f GMICHEM_DeclarePointer___.h GMI_GridComp/Deposition_DeclarePointer___.h
	/bin/mv -f GMICHEM_GetPointer___.h GMI_GridComp/Deposition_GetPointer___.h
	/bin/mv -f GMICHEM_ExportSpec___.h Deposition_ExportSpec___.h
	$(ACG) -N $(NAME) Reactions_Registry___.rc
	/bin/mv -f GMICHEM_DeclarePointer___.h GMI_GridComp/Reactions_DeclarePointer___.h
	/bin/mv -f GMICHEM_GetPointer___.h GMI_GridComp/Reactions_GetPointer___.h
	/bin/mv -f GMICHEM_ExportSpec___.h Reactions_ExportSpec___.h
	$(ACG) -N $(NAME) GMI_Tendency_Registry___.rc
	/bin/rm -f GMICHEM_DeclarePointer___.h
	/bin/rm -f GMICHEM_GetPointer___.h
	/bin/rm -f GMI_Tendency_Registry___.rc
	/bin/mv -f GMICHEM_ExportSpec___.h Tendency_ExportSpec___.h
	/bin/cp -f GMI_GridComp/GmiChemistry/stratTropMech/setkin_chem_mech.txt setkin_chem_mech.txt___.rc

GMICHEM_ExportSpec___.h : Tendency_ExportSpec___.h
	@echo " "
	@echo "Building GMICHEM import and export specs ..."
	$(ACG) -N $(NAME) GMI_GridComp/GMI_Registry.rc
	$(RM) GMICHEM_GetPointer___.h  GMICHEM_DeclarePointer___.h

GMI_GridComp/SCAV_FillExports___.h : GMICHEM_ExportSpec___.h
	./gmi_acg.pl -v
	/bin/mv -f GMI_GridComp/Deposition_GetPointer2___.h GMI_GridComp/Deposition_GetPointer___.h
	/bin/mv -f GMI_GridComp/Reactions_GetPointer2___.h GMI_GridComp/Reactions_GetPointer___.h

# Reduce the optimization for (only)
# GMIchem_GridCompMod.F90 if it cannot compile interactively.
# -----------------------------------------------------------
GMIchem_GridCompMod.o : GMIchem_GridCompMod.F90
	$(FC) -c $(patsubst $(FOPT),$(FOPT2),$(F90FLAGS)) $<

GmiChem_GridCompClassMod.o : GMI_GridComp/GmiChem_GridCompClassMod.F90
	$(FC) -c $(patsubst $(FOPT),$(FOPT1),$(F90FLAGS)) $<

#$(DEPS) : $(ACGS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

