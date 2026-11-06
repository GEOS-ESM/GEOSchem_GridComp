# GEOSchem_GridComp (Chemistry and aerosols): architecture

**Audience:** coding agents and developers working on the Chemistry parent gridded component, its in-repository children and its nested child repositories.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSchem_GridComp` @ `d6ad829` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; key `GEOSchem_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to this repository's root unless they start with `src/` (a fixture path).
Parent: [Physics](../../ARCHITECTURE.md) (the component that runs Chemistry), [GEOSgcm_GridComp overview](../../../../docs/ARCHITECTURE.md), [fixture overview](../../../../../../../docs/ARCHITECTURE.md).
Sibling: [Radiation](../../@GEOSradiation_GridComp/docs/ARCHITECTURE.md).
Child repository with its own document: [GOCART](../@GOCART/docs/ARCHITECTURE.md).
Decision log: [decisions.md](decisions.md). Conventions for agents: [../AGENTS.md](../AGENTS.md).

## 1. What this repository is

It holds the Chemistry parent gridded component (a gridded component is a MAPL/ESMF unit with imports, exports and run methods), registered under the name `CHEMISTRY` in Physics (`../GEOS_PhysicsGridComp.F90`, SetServices, `MAPL_AddChild`). The parent has no chemistry of its own. It adds optional child components, wires them to each other, selects which child provides the radiatively active gases and the aerosol state, and runs the children in two phases (`GEOS_ChemGridComp.F90`, SetServices, Run1, Run2).
Most children live in nested repositories (section 2.3), listed in this repository's `.gitignore` and in the fixture's `components.yaml`.

## 2. Module map

### 2.1 Files at the repository root
| File | Content |
|---|---|
| `GEOS_ChemGridComp.F90` | module `GEOS_ChemGridCompMod`: the `CHEMISTRY` parent (section 3) |
| `GEOS_ChemEnvGridComp.F90` | module `GEOS_ChemEnvGridCompMod`: the always-present `CHEMENV` child (section 4) |
| `GEOS_ChemGridComp.rc` | the `ENABLE_<name>` child switches and `strict_child_timing:` (section 3.2) |
| `ChemEnv.rc` | lightning settings read by `CHEMENV` (section 4) |
| `ChemEnv_ExtData.rc`, `ChemEnv_ExtData.yaml` | ExtData (MAPL external-data component) entries for the lightning inputs `MCOR`, `RATIO_LOCAL`, `MIDLAT_ADJ` |
| `CMakeLists.txt` | build (section 5) |
| `CHANGELOG.md`, `ChangeLog` | change logs (the older `ChangeLog` predates October 2021, per the `CHANGELOG.md` header) |

### 2.2 Directories inside this repository
| Directory | Content |
|---|---|
| `GEOSpchem_GridComp/` | `GEOS_PChemGridComp.F90`, module `GEOS_PChemGridCompMod`: PCHEM, parameterized chemistry. Species are predicted from tabled production and loss or specified from zonally symmetric values, per the module comment. When `RATS_PROVIDER` is `PCHEM` it adds the radiatively active gases and the tendencies `OX_TEND`, `H2O_TEND` as exports (`GEOS_PChemGridComp.F90`, SetServices, `Doing_RATs`). If the Chemistry registry has `doing_PC` false, it registers only the generic (empty) services (`GEOS_PChemGridComp.F90`, SetServices). With `ALLOW_PCHEM_AINC_UPDATE:` /= 0 it registers a second run method, routine `AINC_UPDATE` |
| `HEMCO_GridComp/` | `HEMCO_GridCompMod.F90`: MAPL wrapper of HEMCO (emissions). Run1 does the work, Run2 is empty. Instances are listed by `HEMCO_Instances:` and `HEMCO_CONFIG--<nnn>:` in `HEMCO_GridComp.rc`; the shipped configurations are `HEMCOgocart2g_*` and `HEMCOgmi_*` |
| `DNA_GridComp/` | `DNA_GridCompMod.F90`: imports per-band aerosol extinction, single-scattering albedo and asymmetry (`AEROSOL_EXT_BAND01`..`18` and the matching `_SSA_`, `_ASY_` names) and fills an `AERO` state for Radiation. |
| `GEOSCHEMchem_GridComp/` | build wrapper of GEOS-Chem: compiles sources from the nested `@geos-chem` repository and copies `Interfaces/GCHP/Chem_GridCompMod.F90` (module `GEOSCHEMchem_GridCompMod`) into the build tree (`GEOSCHEMchem_GridComp/CMakeLists.txt`). `GEOSCHEMchem_Registry.rc` and the Perl scripts `gcIncAft`, `gcUtImp`, `gcUtInt`, `gcIncBef` generate include files for the interface (`gcIncAft` header comment; `GEOSCHEMchem_GridComp/CMakeLists.txt`, `new_esma_generate_automatic_code`). In GEOS the child is named `GEOSCHEMCHEM` and registers Run1 and Run2 (`@geos-chem/Interfaces/GCHP/Chem_GridCompMod.F90`, SetServices, `MODEL_GEOS`) |
| `Shared/Chem_Base/` | library `Chem_Base`: the chemistry registry (`Chem_RegistryMod.F90`, reads `Chem_Registry.rc`), `Chem_Mod.F90`, Mie tables and AOD tools, and the executables listed in its `CMakeLists.txt` |
| `Shared/Chem_Shared/` | library `Chem_Shared`: helpers shared by children and by Physics/Agcm, including `Bundle_IncrementMod.F90` (tracer increments), `Chem_GroupMod.F90` (chemical families transported by the Agcm), `Lightning_mod.F90`, deposition, emission and wet-removal modules |
| `Shared/HEMCO/` | `CMakeLists.txt` that compiles the HEMCO sources from the nested `@HEMCO` repository (globbed from `src/Core`, `src/Extensions`, `src/Interfaces/...`) |

### 2.3 Nested child repositories
Each is a separate repository in `components.yaml` at `.../@GEOSchem_GridComp/@<name>` (the `GEOSchem_GridComp` `.gitignore` excludes the `@<name>`, `<name>` and `<name>@` spellings). Descriptions come from each module's header comment (*Unverified:* against behaviour).

| Repository | Module (SetServices) | Child name in Chemistry | What it is |
|---|---|---|---|
| `@GOCART` | `GOCART2G_GridCompMod`, `GOCART_GridCompMod` | `GOCART2G`; `GOCART`, `GOCART.data` | aerosol model; GOCART-2G holds the instance children dust, sea salt, sulfate, carbonaceous, nitrate; GOCART (first generation) also holds CO, CO2, CFC, ozone, radon; see [its document](../@GOCART/docs/ARCHITECTURE.md) |
| `@ACHEM` | `ACHEM_GridCompMod` | `ACHEM` | gas and aqueous phase chemistry (OCS, sulfur, secondary organic aerosol precursors); see [its document](../@ACHEM/docs/ARCHITECTURE.md) |
| `@GMI` | `GMIchem_GridCompMod` | `GMICHEM` | Global Modeling Initiative combined troposphere/stratosphere chemistry; see [its document](../@GMI/docs/ARCHITECTURE.md) |
| `@StratChem` | `StratChem_GridCompMod` | `STRATCHEM` | stratospheric chemistry; see [its document](../@StratChem/docs/ARCHITECTURE.md) |
| `@CARMA` | `CARMAchem_GridCompMod` | `CARMA` | aerosol and cloud microphysics; see [its document](../@CARMA/docs/ARCHITECTURE.md) |
| `@MAM` | `MAM_GridCompMod` | `MAM`, `MAM.data` | modal aerosol microphysics; see [its document](../@MAM/docs/ARCHITECTURE.md) |
| `@MATRIX` | `MATRIX_GridCompMod` | `MATRIX` | aerosol microphysics; see [its document](../@MATRIX/docs/ARCHITECTURE.md) |
| `@GAAS` | `GAAS_GridCompMod` | `GAAS` | GEOS aerosol analysis system; see [its document](../@GAAS/docs/ARCHITECTURE.md) |
| `@TR` | `TR_GridCompMod` | `TR` | passive tracers (many `TR_GridComp---<tracer>.rc` files); see [its document](../@TR/docs/ARCHITECTURE.md) |
| `@RRG` | `RRG_GridCompMod` (file `RRGMod.F90`) | `RRG` | carbon gases (CO, CO2, CH4) production and loss; see [its document](../@RRG/docs/ARCHITECTURE.md) |
| `@QuickChem` | `QuickChem_GridCompMod` | none | not wired in: no `MAPL_AddChild`, not in `CMakeLists.txt` (wave-1 documents record the same); see [its document](../@QuickChem/docs/ARCHITECTURE.md) |
| `@HEMCO` (under `Shared/HEMCO`), `@geos-chem` (under `GEOSCHEMchem_GridComp`) | none | none | upstream codes compiled into `HEMCO_GridComp` and `GEOSCHEMchem_GridComp` |

`H2O` is added as a child when `ENABLE_H2O` is true (`GEOS_ChemGridComp.F90`, SetServices) but no `H2O_GridComp` directory exists here. The build generates a stub for it (`CMakeLists.txt`, `alldirs`; fixture `@cmake/esma_support/esma_add_library.cmake`, `esma_create_stub_component`).

## 3. The Chemistry parent (`GEOS_ChemGridComp.F90`)

### 3.1 Place in the hierarchy
- Parent: Physics adds it as `CHEMISTRY` with `AChemSetServices` (the local name of this module's `SetServices`) (`../GEOS_PhysicsGridComp.F90`, SetServices). Physics runs it in two phases ([Physics document](../../ARCHITECTURE.md), run order).
- Children: added with `MAPL_AddChild` in SetServices, in this order: `CHEMENV` (always), then if enabled `HEMCO`, `PCHEM`, `ACHEM`, `GOCART`, `GOCART.data`, `GOCART2G`, `RRG`, `GAAS`, `H2O`, `STRATCHEM`, `GMICHEM`, `CARMA`, `GEOSCHEMCHEM`, `MATRIX`, `MAM`, `MAM.data`, `TR`, `DNA` (`GEOS_ChemGridComp.F90`, SetServices). The run loops go over `GCS(1:NCHLD)` (Run1, Run2), so *Inferred:* children run in the order they were added.
- The twin names `GOCART`/`GOCART.data` and `MAM`/`MAM.data` use the same `SetServices`. The suffix `.data` selects data-driven (ExtData-fed) behaviour inside `GOCART`: SetServices stores `isDataDrivenGC_` (true when the component name ends in `.data`) in `state%data_driven`, which later selects the registry file `GOCARTdata_AerRegistry.rc` and other branches (`@GOCART/ESMF/GOCART_GridComp/GOCART_GridCompMod.F90`, SetServices, `isDataDrivenGC_`). `MAM` stores the same flag in `self%data_driven` but nothing reads it, so `MAM` and `MAM.data` run the same code (`@MAM/MAM_GridCompMod.F90`, SetServices).

### 3.2 Which children exist: `ENABLE_<name>` keys
The parent loads a private config file, `GEOS_ChemGridComp.rc` by default, or the file named by `GEOS_ChemGridComp_RC_File:` in the main config (`GEOS_ChemGridComp.F90`, SetServices). Every `ENABLE_*` key defaults to `.FALSE.` in code. The shipped `GEOS_ChemGridComp.rc` sets `.TRUE.` for `ENABLE_PCHEM`, `ENABLE_GOCART`, `ENABLE_GOCART_DATA`, `ENABLE_GOCART2G`, `ENABLE_ACHEM`, `ENABLE_TR`, `ENABLE_HEMCO`.

| Key | Child |
|---|---|
| `ENABLE_PCHEM`, `ENABLE_ACHEM`, `ENABLE_HEMCO`, `ENABLE_TR`, `ENABLE_DNA`, `ENABLE_RRG`, `ENABLE_GAAS`, `ENABLE_H2O`, `ENABLE_CARMA`, `ENABLE_MATRIX` | the child of the same name |
| `ENABLE_GOCART`, `ENABLE_GOCART_DATA`, `ENABLE_GOCART2G` | `GOCART`, `GOCART.data`, `GOCART2G` |
| `ENABLE_STRATCHEM`, `ENABLE_GMICHEM`, `ENABLE_GEOSCHEM` | `STRATCHEM`, `GMICHEM`, `GEOSCHEMCHEM` |
| `ENABLE_MAM`, `ENABLE_MAM_DATA` | `MAM`, `MAM.data` |

The setup scripts edit this file by `awk` after copying it to the experiment directory: PCHEM on when `RATS_PROVIDER` is PCHEM, TR on, `ENABLE_GOCART_DATA` off when `AERO_PROVIDER` is GOCART2G, and further switches in `gmichem_setup`, `stratchem_setup`, `geoschemchem_setup` (`src/Applications/@GEOSgcm_App/gcm_setup`, `gmichem_setup`, `stratchem_setup`, `geoschemchem_setup`). `construct_extdata_yaml_list.py` reads the same keys to choose the child ExtData files and always adds `ChemEnv_ExtData.yaml` (`src/Applications/@GEOSgcm_App/construct_extdata_yaml_list.py`).
`GAAS` requires `GOCART2G`: SetServices asserts it and raises an error again where the `AERO` connection to `GAAS` is made (`GEOS_ChemGridComp.F90`, SetServices). The choice is read once, in SetServices.

### 3.3 Entry points
- **SetServices**: registers `Init`, `Run1`, `Run2`; loads the chemistry registry (`Chem_RegistryCreate`, default file `Chem_Registry.rc`); adds children; declares the parent's own specs and the promoted exports (3.4); makes the sibling connections (3.5); calls `MAPL_GenericSetServices` last. There is no Finalize entry point.
- **Init**: `MAPL_GenericInitialize` (initializes the children), then fills the internal bundle `CHEM_TRACERS` from the children's friendly fields, then initializes the increment bundles (3.6). If no aerosol provider is selected it sets two attributes on the exported `AERO` state to false (`GEOS_ChemGridComp.F90`, Init).
- **Run1** (Physics calls it as phase 1, after Surface phase 1): acts only when the component's run alarm is ringing. It does not clear the alarm. It skips the whole phase, with a printed message, when the `TH` import sums to exactly zero. Otherwise it resets the increment bundle `CHEMTRIr1`, calls phase 1 of every child that has more than one run phase (`NumRunPhases > 1`), then computes the phase-1 increments.
- **Run2** (phase 2, after Turbulence phase 2 and before Radiation): clears the run alarm, resets `CHEMTRI`, calls phase 2 of every multi-phase child and phase 1 of every single-phase child, then computes `CHEMTRI`.
- Both phases time each child (`MAPL_TimerOn` with the child name) and optionally put a barrier before and after each child when `strict_child_timing:` is true (`GEOS_ChemGridComp.F90`, Run1, Run2).
- The run alarm is the MAPL generic one made in `MAPL_GenericInitialize` from `DT:` with default `RUN_DT:` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `handle_clock_and_main_alarm`). The component-name prefix rule for resource keys (`src/Shared/@MAPL/base/MAPL_Resource.F90`, `get_labels_with_prefix`) means *Inferred:* the key `CHEMISTRY_DT:` sets the parent interval. The templates set only child intervals: `GOCART2G_DT`, `GOCART_DT`, `HEMCO_DT`, `ACHEM_DT`, `PCHEM_DT` from `@CHEM_DT` (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`). Several children test the same alarm and return early (`HEMCO_GridCompMod.F90`, Run1; `@GOCART/ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, Run1, Run2).
- Phases per child (from `MAPL_GridCompSetEntryPoint`): two phases in `CHEMENV`, `HEMCO` (phase 2 empty), `GOCART`, `GOCART2G`, `GMICHEM`, `RRG`, `GEOSCHEMCHEM`; one phase in `STRATCHEM`, `ACHEM`, `TR`, `MAM`, `MATRIX`, `CARMA`, `GAAS`, `DNA`; `PCHEM` has one, or two with `ALLOW_PCHEM_AINC_UPDATE` /= 0. *Inferred:* phase numbers follow registration order, so `AINC_UPDATE` is then phase 2 and `PCHEM` Run moves from Run2 to Run1 of the parent.

### 3.4 Providers: who supplies what to Physics
The parent re-exports child exports as its own with `MAPL_AddExportSpec(..., CHILD_ID=...)`. The provider is chosen by name in the main config (`AGCM.rc`), resolved by `GetProvider_`, which accepts only names of enabled children (`PCHEM`, `GOCART.data`, `GOCART`, `GOCART2G`, `GAAS`, `H2O`, `STRATCHEM`, `GMICHEM`, `CARMA`, `GEOSCHEMCHEM`, `MATRIX`, `MAM`, `MAM.data`, `TR`, `DNA`, `RRG`) and raises an error otherwise (`GEOS_ChemGridComp.F90`, `GetProvider_`).

| Key (main config) | Default in code | What it selects |
|---|---|---|
| `RATS_PROVIDER:` | `PCHEM` | provider of the radiatively active tracers (RATs: `OX`, `O3`, `O3PPMV`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22`) and of `H2O_TEND` |
| `ANALYSIS_OX_PROVIDER:` | `PCHEM` | provider of `OX`, `O3`, `O3PPMV` (overrides `RATS_PROVIDER` for these three); `OX_TEND` comes from the `OX` provider |
| `<SPECIES>_PROVIDER:` for `OX`, `O3`, `O3PPMV`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22` | `none` (no override) | per-species override |
| `CO2_PROVIDER:` | `None` | `None`, `GOCART` or `RRG`; exports `CO2` from that child |
| `RATS_DIAGNOSTICS:` | empty | extra species names exported from `RRG` for the long-wave diagnostics |
| `AERO_PROVIDER:` | `GOCART2G` | child exporting `AERO` and `AERO_DP`; `none` selects the parent's own empty `AERO` state |

(`GEOS_ChemGridComp.F90`, SetServices.) Also exported: `CO2SC001` from `GOCART` when `ENABLE_GOCART` and the registry has `doing_CO2`.

- Radiation receives `OX`, `O3`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22`, optionally `CO2` (when the provider is `GOCART` or `RRG` and `CO2: -2`) and the `AERO` state; Moist receives `AERO`; Surface receives `AERO_DP` and, when ocean biogeochemistry or `ATM_CO2 == 4` is on, `CO2SC001` (`../GEOS_PhysicsGridComp.F90`, SetServices). Physics exports `OX` and `O3PPMV` from here and reads `H2O_TEND`, `OX_TEND` during Run.
- `AERO` is an ESMF state item, not a field list. The provider fills it with optics and activation methods. Radiation reads the attribute `implements_aerosol_optics_method` on it (`@GEOSradiation_GridComp/GEOSirrad_GridComp/GEOS_IrradGridComp.F90`). GOCART-2G, CARMA, GMI, MAM and DNA set it true (`@GMI/GMIchem_GridCompMod.F90`, the routine that builds the `AERO` state); with no provider, Init sets it false and also sets `implements_aerosol_activation_properties_method` false (`GEOS_ChemGridComp.F90`, Init).
- Aerosol optics method labels differ between providers. GOCART2G (top-level state), CARMA and GMI attach the method under the label `run_aerosol_optics` (`@GOCART/ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`; `@CARMA/CARMAchem_GridCompMod.F90`; `@GMI/GMIchem_GridCompMod.F90`; all `ESMF_MethodAdd`). MAM, DNA, the GEOS-Chem aero coupler and the GOCART2G instance children (DU2G, SS2G, CA2G, NI2G, SU2G, on their own states) attach one under `aerosol_optics` (`@MAM/MAM_GridCompMod.F90`; `DNA_GridComp/DNA_GridCompMod.F90`; `GEOSCHEMchem_GridComp/@geos-chem/Interfaces/GEOS/geos_aerocoupler.F90`; `@GOCART/ESMF/GOCART2G_GridComp/DU2G_GridComp/DU2G_GridCompMod.F90` and its siblings; all `ESMF_MethodAdd`). Radiation executes `run_aerosol_optics` (`@GEOSradiation_GridComp/GEOSirrad_GridComp/GEOS_IrradGridComp.F90`), and so does the GMI photolysis code, which reads the attribute and executes the method on the `AERO` state (`@GMI/GMI_GridComp/GmiPhotolysis_GridCompClassMod.F90`). Observation 8 of the [Radiation document](../../@GEOSradiation_GridComp/docs/ARCHITECTURE.md) records the same mismatch; see Observation 10 here.
- The `AERO_ACI` export exists only when the preprocessor symbol `ENABLE_AERO_ACI` is defined; a search of the fixture's top `CMakeLists.txt`, `@cmake`, `presets` and the `GEOSphysics_GridComp` tree found no definition (`GEOS_ChemGridComp.F90`, SetServices).
- The Gcm writes `AERO_PROVIDER` back into the shared config as `GOCART.data` when the key is missing (`src/Components/@GEOSgcm_GridComp/GEOS_GcmGridComp.F90`, SetServices). *Inferred:* Chemistry's own default `GOCART2G` then applies only when a driver other than the Gcm is used. The template and `gcm_setup` set `AERO_PROVIDER` explicitly (`AGCM.rc.tmpl`; `gcm_setup`, default `GOCART2G`).

### 3.5 Connections between siblings
Made in SetServices with `MAPL_AddConnectivity`, always guarded by the `ENABLE_*` flags of both ends (`GEOS_ChemGridComp.F90`, SetServices). Full name lists are in that routine. The ones that matter:
- `CHEMENV` to children: `AIRDENS`, `DELP` to nearly every child; `AIRDENS_DRYP` to `GOCART`, `STRATCHEM`, `TR`; `CN_PRCP`, `NCN_PRCP` to several aerosol, chemistry and tracer children; `LFR` to `STRATCHEM`, `GMICHEM`, `GEOSCHEMCHEM`; `TPREC`, `LIGHT_NO_PROD` to `GMICHEM`; `BYNCY` to `GEOSCHEMCHEM` and `HEMCO`. `PCHEM` and `DNA` get none.
- Ozone: `PCHEM` `OX` to `GOCART` `O3` (CFC-12 photolysis), `PCHEM` `O3` to `RRG` and `ACHEM`, `PCHEM` `TO3` and `O3` (as `PCHEM_O3`) to `GEOSCHEMCHEM`. Other sources of the same `O3` import are coded per pair: `STRATCHEM` `O3CHEM` to `GOCART` and `ACHEM`, `GMICHEM` `OX` to `GOCART` and `ACHEM`, and `GEOSCHEMCHEM` `TRC_O3` to `GOCART` only when `PCHEM` is off.
- Oxidants and photolysis: `GMICHEM` (`OH`, `H2O2`, `NO3`, `HNO3`, `CH4`, `OX`) to `GOCART`, `GOCART2G`, `ACHEM`, `MAM`; `STRATCHEM` (`O3P`, `OH`, `O3CHEM`, `OCS_JRATE`) to `ACHEM`.
- Aerosol precursors and surface area: `ACHEM` to `MAM`, `MATRIX`, `CARMA`, `GOCART2G` (gas and secondary-organic-aerosol species); `GOCART` or `GOCART2G` to `STRATCHEM`, `GMICHEM`, `CARMA`, `GAAS` (aerosol state or sulfate); `CARMA` to `GMICHEM`/`STRATCHEM` (surface area) when `AERO_PROVIDER` is `CARMA`.
- Emissions: `HEMCO` to `GOCART2G` (`OC_ISOPRENE`, `OC_MTPA`, `OC_MTPO`, `OC_LIMO`) and, if GMI uses MEGAN through HEMCO, `GMI_ISOPRENE` to `GMICHEM`.
- A connection between siblings whose specs are equal is a direct connection: the importer uses the exporter's field (no copy, no lag). A coupler is created only where the specs differ (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_WireComponent`, "this a direct connection" versus "coupler is needed"). Chemistry's Run1 and Run2 do not call `MAPL_GenericRunCouplers` or `MAPL_GenericRunChildren` (`GEOS_ChemGridComp.F90`). *Inferred:* a coupler between two Chemistry children, if one were created, would therefore never run. Physics's calls `MAPL_GenericRunCouplers(STATE, CHEM, ...)` after each Chemistry phase run the couplers from Chemistry's exports to other Physics children (`../GEOS_PhysicsGridComp.F90`, Run).

### 3.6 Tracers, increments and imports
- **Imports**: the parent declares none itself (`MAPL_AddImportSpec` does not occur in `GEOS_ChemGridComp.F90`). Child imports that no sibling connection satisfies are added to the parent's import state by `MAPL_GenericSetServices` (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_WireComponent`: "Imports that are not internally satisfied have their specs put in the GC's import spec"). Physics connects those names from Moist, Turbulence, Surface and Radiation (`../GEOS_PhysicsGridComp.F90`, SetServices, "Chemistry Imports"). Which names are real therefore depends on which children are enabled. With `SYNCTQ` equal to 1 Physics terminates `T`, `TH` and fills them by hand at stages 1 and 3.
- **Friendly tracers**: Init gathers fields whose `FriendlyTo...` attribute names `ANALYSIS`, `TURBULENCE`, `DYNAMICS`, `MOIST`, `STRATCHEM`, `GMICHEM`, `GEOSCHEMCHEM` or `TR` from all children into the internal bundle `CHEM_TRACERS`, with `AddGCPrefix=.true.` so each field is renamed `<child>::<field>` (`GEOS_ChemGridComp.F90`, Init; `src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GridCompGetFriendlies0`). The renamed fields reference the same arrays (default `DoCopy=.false.`, `src/Shared/@MAPL/base/Base/Base_Base_implementation.F90`, `MAPL_FieldCreateRename`).
- **Reaching Physics**: Physics calls `MAPL_GridCompGetFriendlies` on Chemistry without `AddGCPrefix`, and MAPL recurses into children only when that option is true (`../GEOS_PhysicsGridComp.F90`, Initialize, `MAPL_GridCompGetFriendlies(GCS(CHEM), "TURBULENCE", ...)`; `src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GridCompGetFriendlies0`, comment "allow recursing ONLY when addGCPrefix"). *Inferred:* the children's tracers therefore reach Physics through Chemistry's internal bundle `CHEM_TRACERS`, whose fields keep the `FriendlyTo...` attributes of the originals. Physics builds `TR`, `TRADV`, `TRANA` from them ([Physics document](../../ARCHITECTURE.md), section 4.3).
- **Increments**: the parent exports bundles `CHEMTRI` and `CHEMTRIr1` (declared `MAPL_BundleItem`). The species in them are the rows of the table `CHEMTRI_increments::` in the main config; a name absent from `CHEM_TRACERS` stops Init (`Shared/Chem_Shared/Bundle_IncrementMod.F90`, `Initialize_IncBundle_init`). Run1 stores the tracer values before phase 1 and writes `(after - before)/DT` into `CHEMTRIr1`. Run2 does the same around phase 2 and adds the phase-1 value (type `CHMincR2`, `Compute_IncBundle`). Physics builds its 2D mass-weighted `MCHEMTRI` separately from Moist tracers (`../MBundle_IncrementMod.F90`).

### 3.7 State and restarts
The parent keeps a private Fortran state with the `enable_*` flags, `strict_child_timing` and the provider indices `RATS_PROVIDER`, `AERO_PROVIDER`, stored with `ESMF_UserCompSetInternalState` under `GEOSchem_GridComp_State`. Its MAPL internal state holds only the bundle `CHEM_TRACERS`. Restart files belong to the children; the template lists `<CHILD>_INTERNAL_RESTART_FILE` and `_IMPORT_RESTART_FILE` keys for `PCHEM`, `H2O`, `GOCART`, `GOCART.data`, `DU`, `SS`, `CA.*`, `NI`, `SU`, `CARMA`, `MAM`, `ACHEM`, `TR`, `HEMCO` (with the placeholder prefix `@MODELATM`), and keeps those for `STRATCHEM`, `GMICHEM`, `GEOSCHEMCHEM`, `OH`, `RRG` commented out (`src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`). The names given to `MAPL_AddChild` are the prefixes that child resource keys accept (`src/Shared/@MAPL/base/MAPL_Resource.F90`, `get_labels_with_prefix`) and the component names used in HISTORY collections (`GCMPROG.rc.tmpl` uses `'CHEMENV'`).

## 4. `CHEMENV` (`GEOS_ChemEnvGridComp.F90`)
A component with no children that computes the environment the siblings share. SetServices registers `Initialize_`, `Run1`, `Run2`.
- Imports (declared in SetServices): `PLE`, `TH`, `Q`, `T`, `PRECTOT`, `CN_PRCP`, convective mass fluxes and cloud fields, surface state, `ZLE`, `CAPE`, `INHB`, `BYNCY`, `ZLCL`, `ZLFC`, `PPBL`, `TROPP` and others. `PPBL`, `BYNCY`, `CAPE`, `INHB` have `RESTART = MAPL_RestartSkip`. The FIT lightning option adds `MCOR`, `RATIO_LOCAL`, `MIDLAT_ADJ`, fed by `ChemEnv_ExtData.yaml`.
- Run1: dry and moist air density `AIRDENS`, `AIRDENS_DRYP` (`Airdens`), layer thickness `DELP`, precipitation `TPREC`, `CN_PRCP`, `NCN_PRCP`, and lightning `LFR`, `BYNCY`, `CAPE`, `LIGHT_NO_PROD` (`Lightning_mod.F90`, `getLightning`; `LIGHT_NO_PROD` only if something requests it).
- Run2: samples fields at 10 a.m. and 2 p.m. local time into `OVP10_*`, `OVP14_*` exports (`OVP_apply_mask`, `src/Shared/@GMAO_Shared/GEOS_Shared/OVP.F90`). `OVP_init` reads the label `CHEM_DT:` (default `RUN_DT:`).
- Configuration (`ChemEnv.rc`, read by `Lightning_mod.F90`, `read_flash_source`, `read_lightning_config`): `flashSource:` one of `MOIST`, `FIT`, `HEMCO`, `LOPEZ` (code default `LOPEZ`), per-source scale factors `<SOURCE>_flashFactor_resvec_<CTM|FREE|REPLAY>:` by resolution, `useImportedCape:`, `minDeepCloudTop:`, `lightNOampFactor:`, `numberNOperFlash:`. The simulation type is `CTM` when `SimType:` is `CTM`, `REPLAY` when `REPLAY_MODE:` is set, else `FREE` (`GEOS_ChemEnvGridComp.F90`, `Initialize_`).

## 5. Dependencies and build
- `CMakeLists.txt` adds subdirectories `Shared GOCART TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM` (found under `name`, `@name` or `name@`) and lists `GEOSpchem_GridComp GEOSCHEMchem_GridComp H2O_GridComp DNA_GridComp HEMCO_GridComp RRG` as sub-components (a missing one is replaced by a generated stub). The library links `MAPL Chem_Shared Chem_Base GOCART_GridComp GOCART2G_GridComp TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM RRG ESMF::ESMF`. It installs `GEOS_ChemGridComp.rc`, `ChemEnv.rc`, `ChemEnv_ExtData.rc`, `ChemEnv_ExtData.yaml` into `etc`.
- The parent source also uses `Chem_Mod`, `Chem_UtilMod`, `Bundle_IncrementMod` and, in the fixture, `GMAO_Shared` and `MAPL` (`GEOS_ChemGridComp.F90`, `use` lines).
- Physics and the Agcm depend on this repository: Physics uses `Bundle_IncrementMod`, the Agcm uses `Chem_GroupMod` (`../../GEOS_AgcmGridComp.F90`, `Pack_Chem_Groups`). A checkout without this repository does not build Physics.
- `CMakeLists.txt` calls macros from the fixture's `@cmake` (`esma_set_this`, `esma_add_library`). *Inferred:* the repository builds only inside a fixture that provides them. `esma_mepo_style` finds `@HEMCO` and `@geos-chem`.
- The registry file `Chem_Registry.rc` is read from the run directory at SetServices (`Shared/Chem_Base/Chem_RegistryMod.F90`, `Chem_RegistryCreate`); the shipped copy turns `doing_PC` and the CO, CO2, CFC groups on and the aerosol groups off (`Shared/Chem_Base/Chem_Registry.rc`). The parent uses `doing_CO2`, `doing_OCS`, `doing_SU` (and `i_SU`..`j_SU`), `doing_BC`, `doing_DU`, `doing_OC`, `doing_SS`.

## 6. Parent and child repositories
- Parent repository: `GEOSgcm_GridComp` ([overview](../../../../docs/ARCHITECTURE.md)); component parent: [Physics](../../ARCHITECTURE.md). Fixture: [GEOSgcm](../../../../../../../docs/ARCHITECTURE.md), section 2.5.
- Configuration templates and setup scripts that edit the `.rc` files live in `GEOSgcm_App` (`src/Applications/@GEOSgcm_App`; its own document is written in that repository).
- Children with documents: [GOCART](../@GOCART/docs/ARCHITECTURE.md), [ACHEM](../@ACHEM/docs/ARCHITECTURE.md), [CARMA](../@CARMA/docs/ARCHITECTURE.md), [GAAS](../@GAAS/docs/ARCHITECTURE.md), [GMI](../@GMI/docs/ARCHITECTURE.md), [MAM](../@MAM/docs/ARCHITECTURE.md), [MATRIX](../@MATRIX/docs/ARCHITECTURE.md), [QuickChem](../@QuickChem/docs/ARCHITECTURE.md), [RRG](../@RRG/docs/ARCHITECTURE.md), [StratChem](../@StratChem/docs/ARCHITECTURE.md), [TR](../@TR/docs/ARCHITECTURE.md).

## 7. Known traps
1. **`providerName` is reused.** After the `AERO_PROVIDER` block it holds the aerosol provider's name; the later tests `TRIM(providerName) == "GOCART2G"`, `"GOCART.data"`, `"GOCART"`, `"CARMA"` select the GMI/CARMA/GOCART couplings by that name (`GEOS_ChemGridComp.F90`, SetServices).
2. **SetServices writes to the shared config.** With `ENABLE_CARMA`, it sets `SULFURIC_ACID_SOURCE:` to `tendency`, or `full_field` if GMI is also on; `CARMA` and `GMI` read the key later with other defaults (`tendency`, `GMIvalue`) (`GEOS_ChemGridComp.F90`, SetServices; `@CARMA/CARMAchem_GridCompMod.F90`; `@GMI/GMIchem_GridCompMod.F90`).
3. **The parent reads child resource files.** With `ENABLE_HEMCO` it opens `HEMCO_GridComp.rc`; with `ENABLE_HEMCO` and `ENABLE_GMICHEM` it also opens `GMI_GridComp.rc`, both from the run directory, and stops (status 98 or 99) if the HEMCO instance list and the GMI MEGAN flags disagree (`GEOS_ChemGridComp.F90`, SetServices).
4. **The `ENABLE_*` set changes the import list.** Physics connections name fields by string; adding or removing a child changes which of them have a destination (3.6).
5. **Run order is set by the parent.** Physics runs Chemistry phase 1 after Surface phase 1 and before Turbulence phase 1, and phase 2 after Turbulence phase 2 and before Radiation (`../GEOS_PhysicsGridComp.F90`, Run). A new child with one run phase runs only in Chemistry's Run2.
6. **First-call skip.** When `TH` is all zero, Run1 returns before any child runs, but Run2 still runs phase 2 of multi-phase children (`GEOS_ChemGridComp.F90`, Run1, Run2).
7. **Friendly names are exact strings.** `CHEM_TRACERS` collects only the eight names in 3.6; a new friendly name needs a change in Init and in Physics.
8. **`CHEM_TRACERS` entries are references.** An increment is correct only because the same arrays are read before and after the child runs.
9. **PCHEM as `ANALYSIS_OX_PROVIDER` needs `OX_FRIENDLIES` to include `ANALYSIS`**, otherwise it stops (`GEOSpchem_GridComp/GEOS_PChemGridComp.F90`, SetServices).

## 8. Observations (noted, not fixed)
1. `ENABLE_H2O` and `GetProvider_` name an `H2O` child whose directory does not exist here (section 2.3).
2. `RATS_DIAGNOSTICS:` entries are exported with `CHILD_ID = RRG` without checking `ENABLE_RRG`; `RRG` is -1 when the child is not added (`GEOS_ChemGridComp.F90`, SetServices).
3. `AERO_PROVIDER` has defaults of `GOCART2G` in Chemistry, `GOCART.data` in the Gcm and `GOCART2G` in `gcm_setup`. `IDFUPD.rc.tmpl` sets `AERO_PROVIDER: PCHEM`, a name `GetProvider_` accepts, while `GEOS_PChemGridComp.F90` declares no `AERO` export (`src/Applications/@GEOSgcm_App/IDFUPD.rc.tmpl`). The comment in `AGCM.rc.tmpl` beside `AERO_PROVIDER:` lists the options `GOCART2G, MAM, none`, while `GetProvider_` accepts every enabled child name (3.4).
4. `chemReg` created in SetServices is not destroyed there (`Chem_RegistryDestroy` is not called in the file).
5. Run2's `DEBUG` block lists `speciesName` in a different order from SetServices (`O3PPMV` last versus third); only printed names are affected.
6. The Init comment calls overriding child attributes when `AERO_PROVIDER` is `none` "an architecture violation" (`GEOS_ChemGridComp.F90`, Init).
7. `OVP_init` reads `CHEM_DT:`; the template defines `GOCART2G_DT`, `HEMCO_DT`, ... from `@CHEM_DT`, not a `CHEM_DT:` key (`AGCM.rc.tmpl`; grep of `src/Applications`).
8. `DNA_GridComp.rc` has `number_of_bands: 19` while `DNA_GridCompMod.F90` loops over 18 bands and does not read that key. `DNA_ExtData.yaml` lists absolute template paths under `/home/` of an individual user.
9. Questions for maintainers: What is the role of the `H2O` child today? Is `QuickChem` meant to be added as a Chemistry child? What happens if two Chemistry children are connected with unequal specs, given that Chemistry does not run its internal couplers? Which of the three `AERO_PROVIDER` defaults applies in which driver?
10. Aerosol optics method labels: MAM, DNA and the GEOS-Chem aero coupler attach `aerosol_optics` on the `AERO` state, while Radiation executes `run_aerosol_optics` (section 3.4). Which label is intended for each provider, and how does a provider with only `aerosol_optics` reach Radiation?
