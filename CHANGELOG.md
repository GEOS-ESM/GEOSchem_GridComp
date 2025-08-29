# Changelog

All notable changes to this project will be documented in this file
after October 2021. For earlier changes, consult file ChangeLog.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added

- Added the capability to compute LFC in Lightning module, needed for CTM

### Removed
### Changed

- Slight refactoring of LOPEZ lightning scheme, minor numerical difference
- Renamed lightning option 'usePreconCape' to 'useImportedCape' in order to cover both GCM and CTM

### Fixed

- ChemEnv now updates the precip exports (total, conv, non-conv) as part of Run1 (not just Run2).

### Deprecated


## [1.16.0 - 2025-02-18]

### Added

- Added the RRG child component as a refactor of an instance-based system for GHG modeling
- Added the ability to specify RATS providers individually, including for CO2 (via CO2_PROVIDER)
- Added SO4REFF connectivity from CARMA to GMI
- Added connectivity (OH, H2O2, NO3) from GMI to GOCART

### Removed

- Removed `MAMchem_GridComp`, `MATRIXchem_GridComp`, `CARMAchem_GridComp`, `GAAS_GridComp`, and `GEOSachem_GridComp` as these are now in separate repos
   - `MAMchem_GridComp` → [MAM](https://github.com/GEOS-ESM/MAM)
   - `MATRIXchem_GridComp` → [MATRIX](https://github.com/GEOS-ESM/MATRIX)
   - `CARMAchem_GridComp` → [CARMA](https://github.com/GEOS-ESM/CARMA)
   - `GAAS_GridComp` → [GAAS](https://github.com/GEOS-ESM/GAAS)
   - `GEOSachem_GridComp` → [ACHEM](https://github.com/GEOS-ESM/ACHEM)


## [1.15.0] - 2025-01-16

### Changed

- The file path was changed for anthropogenic CO emissions that are used by achem. Note that the previous version of the emissions have an incorrect seasonal cycle.
- Update ESMF CMake target to `ESMF::ESMF`
- Overhauled the Lopez lightning scheme, and made it the default scheme; note that lightning is used by GMI for computing NOx emissions; PCHEM does not use lightning

### Fixed

- Updated GAAS_Gridcomp_Extdata.yaml in AMIP/ to avoid the model to crash when GAAS is turned on and AMIP emissions chosen.


## [1.14.0] - 2024-06-21

### Added

- Connectivity from GMI to ACHEM (4 fields), requires GMI v1.2.0 or later to run
- Slight improvement for lightning flash rate calculation (LOPEZ and MOIST schemes). See the option UsePreconCape in ChemEnv.rc . This involves new imports from MOIST: CAPE, BYNCY and INHB. **NOTE** THIS REQUIRES GEOSgcm_GridComp develop branch (as of 12/12/23).
- Added a flag for 'strict' child timing, intended to reduce the timing bias against child GC's that employ 'gather' calls. Such calls are occasionally necessary, but can cause timers to attribute excessive time to a child, time that is actually the synchronization lag time that would eventually be spent -somewhere- in the program, but which gets attributed to the child with a 'gather' or barrier call. The new flag is for timing tests only.

### Changed

- Update CI to use v2 orb

### Fixed

- Fixed CARMA to fix radiation callback
- Fixed code in CARMA to properly check whether GMI or GOCART are providing sulfur inputs
- Fixed CARMA/GOCART2G sulfate production tendency term
- Fix a bug in GAAS where it gets the VM (global instead of the correct current)
- Fix an issue in GAAS where the `aod_?` fields were not declared as `MAPL_RestartSkip` in the Registry file.

## [1.13.1] - 2023-04-24

### Added

- Added CO2 connectivity in GEOS_ChemGridComp for GOCART-GEOS-Chem coupling
- Added a wrapper routine for the MAPL Solar Zenith Angle call, in Chem_Shared
- Added connectivity from GOCART2G aerosols to GMI chem

### Removed

- Removed the GMI routines which computed Solar Zenith Angle, in Chem_Shared; but in a later commit, this was added back temporarily, so that older versions of GMI and TR don't complain.

- Removed parallel read of PChem species file. This parallel read was causing issues at NAS at large node count, so now we just do a read-on-root followed by a broadcast

### Fixed

- Fixed a bug that had prevented GMI running with HEMCO

## [1.13.0] - 2023-03-01

### Added

- Added connectivities in GEOS_ChemGridComp for CARMA-GMI coupling
- Added STS growth in CARMA microphysics

### Changed

- For OPS configuration: removal of links, change of QFED paths from vNRT/ to v2.5r1-nrt/ (note after November 2021, files are v2.6r1) (0 diff)
- For AMIP configuration: update of QFED from v2.5r1 to v2.6r1 (most recent collection, may have small diff)

## [1.12.0] - 2023-01-18

### Removed

- Removed `TR_GridComp`, `GMI_GridComp`, and `StratChem_GridComp` as these are now in separate repos
   - `TR_GridComp` → [TR](https://github.com/GEOS-ESM/TR)
   - `GMIchem_GridComp` → [GMI](https://github.com/GEOS-ESM/GMI)
   - `StratChem_GridComp` → [StratChem](https://github.com/GEOS-ESM/StratChem)

## [1.11.0] - 2023-01-04

### Added

- GMI now exports stOX_loss (stratospheric OX tracer loss), customized for the specific chemical mechanism being run.

### Removed

- GMI lbssad_opt allowed unsupported options, eliminated
- GMI h2oclim_opt had only one viable option (3), removed h2oclim_opt and various unused arrays

### Changed

- Instead of importing a set of QQK diagnostic fields for chemical loss of stOX, TR now imports a single field: stOX_loss
- For OPS configuration: removal of links, change of QFED paths from vNRT/ to v2.5r1-nrt/
- For AMIP configuration: update of QFED from v2.5r1 to v2.6r1 (most recent collection)
- Moved to GitHub Actions for label enforcement
- Update CircleCI to use Baselibs 7.7.0
- For GMI, clean-up related to exports HNO3GASsad, HNO3CONDsad and jNO2val.

### Fixed

- Modified TR to only import stOX_loss if loss_species == OX; without this, a GMI _ASSERT may exit the program needlessly.
- Fixed the long_names for GMI chemical species
- Fixed the long_names for StratChem chemical species
- Updated GAAS to now work again after it was changed to use ExtData, only works with ExtData2G
  - **NOTE 1**: This requires MAPL 2.32 or higher to build as a new procedure had to be created for this to work.
  - **NOTE 2**: As noted above, GAAS will now *only* work with ExtData2G
- Fixed the imports for LAI_FRAC and VEG_FRAC in TR and GMI, for simulations with vertical resolution other than 72.

## [1.10.4] - 2022-11-08

### Removed

- Deleted the old CCMI resource files under GMI.

### Changed

- StratChem is no longer part of Chem_Registry.rc; instead the entries are in SC_Mech_Registry.rc, and there are two versions of the resource file - one for the Full mechanism and one for the Reduced mechanism.  Choice between the two mechanisms remains a compile time option.

### Fixed

- Fix 1d vector of latitudes in GCR emissions, pre-cubesphere leftover

## [1.10.3] - 2022-10-27

### Added

- Capability of multiple GMI mechanisms, choose at compile time
- New GMI mechanism: StratTrop_HFC_S
- Add point emission capability for GMI chemical species in GMIchem

## [1.10.2] - 2022-09-22

### Added

- Species_Bundle and Species_Array modules are simpler versions of the Chemistry counterparts

### Changed

- Updates to emissions from galactic cosmic rays in GMI
- Minor improvement to Runtime_Registry module.
- Broke away the GMI contents from Chem_Registry.rc, into a separate file

## [1.10.1] - 2022-08-30

### Added

- Added new Runtime_Registry module in Shared/Chem_Base, for Chem children to use instead of Chem_Registry.

### Removed

- Removed diagnostic messages for GMI isoprene emissions
- Removed code related to TR from Chem_Registry; now handled with Runtime_Registry

### Changed

- GMI photolysis now uses "random cloud overlap" instead of "maximal overlap"
- GMI now enforces a floor value for transported species
- TR now uses Runtime_Registry instead of Chem_Registry

### Fixed

- Fixed typo in yaml files

## [1.10.0] - 2022-08-16

### Added

- Add YAML validator GitHub Action
  - This action makes sure all YAML files are valid (to a relaxed standard)
- Added flag to control whether GMI feeds back QV value to rest of model

### Fixed

- Fixed Dry Deposition in GMI
- Fixed small memory leak in GMI

### Changed

- Changed GMI from internal SZA calculation to using MAPL SZA.
- Changed GMI to use (CMIP6) emissions and boundary conditions from CCMI REF-D1
- Moved external data files (like emissions) from personal space to GMAO shared space
- Improved diagnostic print statements

## [1.9.6] - 2022-08-04

### Fixed

- Updated CI to work with latest GEOSgcm
- Added QuickChem repo to `.gitignore`

## [1.9.5] - 2022-06-22

### Fixed

- Fix YAML typo in `GEOSachem_GridComp/GEOSachem_ExtData.yaml`

## [1.9.4] - 2022-05-31

### Fixed

- Fixed ExtData2G YAML files for OPS and AMIP emissions to handle magic date and removed input file for AMIP.20C

## [1.9.3] - 2022-05-25

### Added

- ChemEnv now provides several more Overpass diagnostics.  It also now has a (thread safe) internal state.
- Add ExtData2G YAML Files

### Changed

- Several CHEM children no longer use TPREC, removed Connectivity
- Updates to the ConvectionMod in Chem_Shared: reduced dust scavenging by 80%; added temperature-dependent scavenging of Pb, Be species
- Update to the WetRemovalMod in Chem_Shared: added temperature-dependent snow scavenging of Pb, Be species
- Update CircleCI to v1 orb
- More updates to CMake for spack

### Fixed

- Fixed bug related to recent Precip change, to satisfy GMI import
- Corrected the connectivity in GEOS_ChemGridComp.F90 for the pSO2_OCS; destination is now GOCART2G instead of GOCART. This is zero-diff except when running the StratChem/ACHEM/GOCART2G OCS-produced SO2 mechanism (not presently used in FP).

## [1.9.2] - 2022-04-29

### Added

- Added `AMIP.20C` directory to Chem_Base

### Changed

- Cleaned up `AMIP` directory in Chem_Base

### Fixed

- Fixed bug in GEOS_Achem, variables being provided via ExtData did not have restart skip

## [1.9.1] - 2022-03-18

### Changed

- Updates to CMake to support Spack
- Moved CircleCI to use circleci-tools orb
- Code was modified to receive Observed Precip data from Surface GridComp Exports (rather than an ad-hoc READ from within Chem). Therefore, Precip into Chem is based on whatever method is used in Surface.

### Added

- Added Changelog Enforcer GitHub Action

## [1.9.0] - 2022-03-15

### Added

- Added files for new ACHEM scenario: "20th century AMIP"

### Changed

- Modified version of QFED emissions for the ACHEM AMIP scenario.
- Modified GMI to support running without aerosols, since the coupling between GMI and GOCART2G is not working yet.

## [1.8.1] - 2022-02-22

### Fixed

- Fix bug with GNU and Achem Finalize when `aqueous_chemistry: .false.`

## [1.8.0] - 2022-02-14

### Changed

- Renamed HEMCO_GridComp/HEMCOgeosfp_*.rc as
  HEMCO_GridComp/HEMCOgocart2g_*.rc as these are not specific of
  GEOS-FP but rather of any system thans runs with GOCART-2G.
- Pathnames in HEMCO*_Config.rc have been changed to /dev/nulls as files are piped in through ExtData.
- Update rc files (109-->108) to reflect the fact that MEGAN is now a single extension.
- Revised GAAS for compatibility with GOCART-2G
- Chem_GridComp exports both AERO and AERO_ACI states. When AERO PROVIDER
  is GOCART-2G, these 2 states point to the same AERO state that GOCART-2G exports.

### Removed

- HEMCO_GridComp/HEMCOgocart_*.rc as they relate to legacy GOCART
  aerosols which have been removed.

