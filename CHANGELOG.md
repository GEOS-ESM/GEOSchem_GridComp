# Changelog

All notable changes to this project will be documented in this file
after October 2021. For earlier changes, consult file ChangeLog.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
### Removed
### Changed
### Fixed

## [1.21.3] - 2023-10-03

### Fixed

- Fix an issue in GAAS where the `aod_?` fields were not declared as `MAPL_RestartSkip` in the Registry file.

## [1.21.2] - 2023-04-07

### Removed

- Completely remove parallel read of PChem species file.

## [1.21.1] - 2023-04-06

### Removed

- Removed parallel read of PChem species file. This parallel read was causing issues at NAS at large node count, so now we just do a
  read-on-root followed by a broadcast

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

