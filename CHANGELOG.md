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

### Fixed

- Fixed Dry Deposition in GMI

### Changed

- Changed GMI from internal SZA calculation to using MAPL SZA.

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

