# Changelog

All notable changes to this project will be documented in this file
after October 2021. For earlier changes, consult file ChangeLog.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Updates to CMake to support Spack

### Removed

### Added

### Fixed

## [1.8.1] - 2021-02-22

### Fixed

- Fix bug with GNU and Achem Finalize when `aqueous_chemistry: .false.`

## [1.8.0] - 2021-02-14

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

### Added

### Fixed

