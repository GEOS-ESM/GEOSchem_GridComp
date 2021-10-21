# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed 
- Renamed HEMCO_GridComp/HEMCOgeosfp_*.rc as
  HEMCO_GridComp/HEMCOgocart2g_*.rc as these are not specific of
  GEOS-FP but rather of any system thans runs with GOCART-2G.
- Pathnames in HEMCO*_Config.rc have been changed to /dev/nulls as files are piped in through ExtData.
- Update rc files (109-->108) to reflect the fact that MEGAN is now a
  single extension.
- Revised GAAS for compatibility with GOCART-2G

### Removed 
- HEMCO_GridComp/HEMCOgocart_*.rc as they relate to legacy GOCART 
  aerosols which have been removed. 

### Added

### Fixed

## [1.7.0] - 2021-10-18

### Changed

- Wet deposition for CO, CO2 and CH4 is disabled, while a bug-fix is developed. (non-zero-diff)
  - Disabling wet deposition for CO2, CO, and CH4 due to conservation and convective transport issues. The corrected units for the Henry's Law coefficients used in GF led to excessive wet deposition of CO2 and possibly other constituents in the x0046 experiments. This code change has been adopted as a temporary fix in ongoing x-experiments.
- Lightning module has been re-located and revamped.
  - The PR (feature/mmanyin/new_lightning_options) is one of three that must be done together (GEOSgcm_GridComp, GMAO_Shared and GEOSchem_GridComp). The result is zero-diff for PCHEM simulations. It is non-zero-diff for simulations that use Lightning Flash Rate, such as GMI.
  This set of PR's is focused on the Lightning Module, which now computes both Flash Rate and NOx Emissions. The module has been moved from GMAO_Shared to Chem_Shared, since its output is used only within Chemistry code. The lightning module has also been significantly re-worked. Different Flash Rate algorithms are now available.
- HEMCO config files for GOCART and GMI have been updated for v2.2
  - The MEGAN configuration inputs needed by the HEMCO configuration file were changed between HEMCO v2.1 and HEMCO v2.2. In HEMCO v2.2, MEGAN is just one single extension and the configuration file needs to reflect that. C. Keller changed the configuration files for both GOCART and GMI so that they work properly with HEMCO v2.2 and higher.
  This change is zero-diff for default runs (where HEMCO is not being used) but result in non-zero changes if HEMCO is enabled as it now computes MEGAN emissions properly (again). When using the old configuration files with HEMCO v2.2, MEGAN emissions are all zero (unfortunately, this is a silent bug). The updated configuration files fix this problem.title:	HEMCO and GEOS-Chem are now in their own repo's
- GMI and TR have cleaner units handling (non-zero-diff, minor change)

## [1.6.0] - 2021-05-17

### Changed

- GEOS-Chem and HEMCO now exist in separate repositories, as specified by GEOSgcm_GridComp.
The corresponding code from this repository has been removed.  The new code is more up-to-date, so running GEOS-Chem or HEMCO is non-zero-diff.
But all other cases are zero-diff.
- Some changes to accommodate python 2 & 3, along with other bookkeeping. 

### Fixed

- GMI is now set to import RI and RL only when running Cloud-J. This is a work-around for a CTM issue.

## [1.5.0] - 2021-04-15

### Fixed 

- Fixes the problem of negative PCHEM tracers.
  - The interpolation of Production and Loss fields in PCHEM has been changed, to prevent extrapolation (for example in gridboxes where the midpoint pressures are greater than the maximum pressure in the prod/loss climatology). (non-zero-diff)

## [1.4.4] - 2021-03-17

### Removed

- GOCART is no longer included in this repository.
When GEOSgcm includes this release and the new GOCART repo, the result is Zero
Diff.

## [1.4.3] - 2021-03-16

### Changed

- Added Python Transition Team Codeowner.
- Silenced PRINT statements related to OVP.
- Zero diff

## [1.4.2] - 2021-02-18

- GMI has new overpass (OVP) diagnostics, and a fix for H2O_TEND
- PCHEM has a new option for water vapor production & loss (off by default)
- ChemEnv has new overpass diagnostics
- Infrastructure updates pertaining to Circle CI, docker authentication, etc

## [1.4.1] - 2020-08-17

### Changed

* Convert PChem to use nf90 interfaces
* Use xlarge resource for CircleCI
* Use Ubuntu20 GCC10 image

## [1.4.0] - 2020-07-07

### Fixed

- In the shared wet-removal routine, replaced "greater than zero" with "greater than tiny", to avoid usage of extremely small numbers; this is technically non-zero-diff for aerosols, and impact should be minimal.
  - NOTE: cases of this have been seen when using GCC compiler, but not (yet) for Intel compiler

### Changed 

* Revision to the GMI FastJ photolysis, mostly removing outdated versions, and clarifying an index range.
* Minor housekeeping.

## [1.3.5] - 2020-06-09

### Changed

- Keeping CMake files in line with the current practices.

## [1.3.4] - 2020-06-04

### Fixed

* Fixed memory leak in TR
* Photolysis in GMI no longer reports spurious FAIL

### Changed

* New photolysis option in GMI: CloudJ
* Removed extraneous files under GMI
* Replaced rcEsmfReadLogical with ESMF_ConfigGetAttribute in GMI (this allows GMI stubbing)
* New GMI boundary condition file
* 3D emissions diagnostic in TR is now kg/m2/s instead of kg/m2
* Several TR diagnostics now account for moist air
* MAPL 2.1.1 needs 6.0.11
* Update CI to Baselibs 6.0.12

## [1.3.3] - 2020-04-07

### Changed

- Enable use of online lightning flash rates in GEOSCHEMchem.
  - Requires GEOSgcm_GridComp PR #187:  Online calculation of lightning flash rate for GEOSCHEMchem

## [1.3.2] - 2020-03-27

### Changed

- GOCART update for 2-moment tuning

## [1.3.1] - 2020-03-25

### Fixed

- Fixes for CMake using F2Py and LaTeX. This release should be used along with these Pull Requests:
  - ESMA_cmake#62
  - GMAO_Shared#84
  - UMD_Etc#11

## [1.3.0] - 2020-03-25

### Changed

- Upgrade to GMI. GMI can now use HEMCO for isoprene emissions, and can save time if RAS_NO_NEG is in effect.

### Fixed

- A fix in GAAS (in REPLAY mode) now corrects forecasts (This requires GEOSgcm_GridComp PR #241.)

## [1.2.1] - 2020-03-20

### Fixed

- Fixed compilation of `MieObs_.so`. The shared object now compiles properly.

## [1.2.0] - 2020-02-12

### Changed

- Changes needed for MAPL 2 compatibility
  - This release requires MAPL2.0 versions of GEOSgcm_App, GEOSgcm_GridComp,
    FVdycoreCubed_GridComp, fvdycore  and  GMAO_Shared.  Numeric results should
    not change significantly, but with changes affecting more than 100 files
    non-zero-diff is possible

## [1.1.0] - 2020-01-28

### Changed

- Catch up with FPP version
  - Transport updates from f525land_fpp and Jason-3_4 (from CVS) This release is meant to be used in conjunction with updates to other repo's -- see GEOSgcm_GridComp PR#190, GEOSgcm_App PR#84, and GMAO_Shared  PR#70.

## [1.0.6] - 2020-01-14

### Fixed

- GAAS needed an alarm, to work in the context of the latest IAU implementation.  This release must be used in conjunction with pull requests GEOSgcm_App#72 and GEOSgcm_GridComp#168.

## [1.0.5] - 2019-12-10

### Fixed

- Enforcement of `GOCART_DT == HEARTBEAT_DT`  (to avoid problems w/ non-aerosol GOCART species)

### Changed

- Minor updates to GEOS-Chem
- Updated Python code in MAM/optics

## [1.0.4] - 2019-10-22

### Fixed

- Fixes an issue involving array index parameters, set to zero to indicate "unused",  and an unused section of code under GMI that included those parameters in array expressions.  This could not compile with "array bounds checking" enabled.  Our solution was to delete the usused section of code for
now. This addresses #30 and #34 .

## [1.0.3] - 2019-10-18

### Changed

- New approach for increment diagnostics 
  - Tendency diagnostics are now available for any CHEM species, with respect to the following processes: advection, convection, turbulence, chemistry and water vapor rescaling.
- GEOS-Chem has been updated to Harvard version 12.4.0
- GMI now uses parameterized ship emissions for NOx

## [1.0.2] - 2019-08-15

### Fixed

- Fixes gmichem_setup and stratchem_setup files. Also allows for stratchem reduced mechanism feature. These were broken during the transition to Git and cmake.

## [1.0.1] - 2019-07-26

### Fixed

- Bug fix for TR GridCompfixe that now allows running at NAS and on GMAO desktops

## [1.0.0] - 2019-07-23

### Added

- Initial Release with Semantic Versioning
  - Repository split from CVS ESMA Repository.  Equivalent to `cvs/GEOSadas-5_25_0` release. 
