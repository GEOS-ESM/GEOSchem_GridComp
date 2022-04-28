# Changelog

All notable changes to this project will be documented in this file
after October 2021. For earlier changes, consult file ChangeLog.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- This is a one-line bug-fix that corrects the connectivity in GEOS_ChemGridComp.F90 for the pSO2_OCS. It was incorrectly pointing to GOCART as destination where it should have been GOCART2G. THis is zero-diff except when running the StratChem/ACHEM/GOCART2G OCS-produced SO2 mechanism (not presently used in FP)

## [1.9.1] - 2022-03-18

### Changed

- Updates to CMake to support Spack
- Moved CircleCI to use circleci-tools orb
- Code was modified to receive Observed Precip data from Surface GridComp Exports (rather than an ad-hoc READ from within Chem). Therefore, Precip into Chem is based on whatever method is used in Surface.

### Removed

### Added

- Added Changelog Enforcer GitHub Action

### Fixed

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

