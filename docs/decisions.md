# GEOSchem_GridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `GEOSchem_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSchem_GridComp` @ `d6ad829` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[ARCHITECTURE.md](ARCHITECTURE.md).

### Children chosen by run-time switches, not by the build
- Pattern: the Chemistry parent reads one `ENABLE_<name>` key per child from a private resource file and calls `MAPL_AddChild` only for the enabled ones. All children are compiled in; every key defaults to `.FALSE.` in code (`GEOS_ChemGridComp.F90`, SetServices; `GEOS_ChemGridComp.rc`, `ENABLE_*`). [Verified]
- Likely reason: one executable serves several configurations (parameterized, GMI, StratChem, GEOS-Chem, data-driven aerosol), selected by the setup scripts. [Inferred]
- Question for maintainers: the same keys are also read by `construct_extdata_yaml_list.py` to pick ExtData files; why are the two readers separate?

### Providers chosen by name, resolved to child indices
- Pattern: `RATS_PROVIDER`, `ANALYSIS_OX_PROVIDER`, per-species `<SPECIES>_PROVIDER`, `CO2_PROVIDER` and `AERO_PROVIDER` name a child; `GetProvider_` maps the name to a child index and the parent re-exports that child's field or state with `MAPL_AddExportSpec(..., CHILD_ID=...)` (`GEOS_ChemGridComp.F90`, SetServices, `GetProvider_`). Physics then connects the parent's exports to Radiation, Moist and Surface by name (`../GEOS_PhysicsGridComp.F90`, SetServices). [Verified]
- Likely reason: Physics, Radiation and Moist stay unchanged when the chemistry scheme changes. [Inferred]
- Question for maintainers: the comment above the child index variables says a new component must be added to `GetProvider_`; how is the list kept in step with the children?

### Two run phases, with single-phase children folded into phase 2
- Pattern: the parent registers `Run1` and `Run2`. Run1 calls phase 1 of children with more than one run phase; Run2 calls phase 2 of those and phase 1 of the single-phase children (`GEOS_ChemGridComp.F90`, Run1, Run2; `MAPL_Get(..., NumRunPhases=...)`). Physics runs Run1 before Turbulence and Run2 after it (`../GEOS_PhysicsGridComp.F90`, Run). [Verified]
- Likely reason: emissions and deposition in the first phase see the state after Surface; transport-like and chemical updates in the second phase see the state after Turbulence. [Inferred]
- Question for maintainers: a source comment dates the phase logic to 2014 (`GEOS_ChemGridComp.F90`, Run1); is the split between emissions and chemistry the same in every child?

### Run alarm shared by the whole component
- Pattern: Run1 and Run2 act only while the component's MAPL run alarm rings. Run1 leaves the alarm on, Run2 turns it off (`GEOS_ChemGridComp.F90`, Run1, Run2). Children may test their own alarm as well. [Verified]
- Likely reason: Chemistry can run at a longer interval than the heartbeat, and phase 2 must still see the ringing alarm. [Inferred]

### Shared environment computed by a child, not by the parent
- Pattern: `CHEMENV` is always added. It computes `AIRDENS`, `AIRDENS_DRYP`, `DELP`, precipitation, lightning flash rate and the overpass diagnostics, and the parent connects these to the other children (`GEOS_ChemGridComp.F90`, SetServices; `GEOS_ChemEnvGridComp.F90`, Run1, Run2). [Verified]
- Likely reason: the sibling children obtain the same fields and one component owns the lightning options. [Inferred]

### Friendly tracers collected into one internal bundle
- Pattern: Init gathers every child's fields marked `FriendlyTo` `ANALYSIS`, `TURBULENCE`, `DYNAMICS`, `MOIST`, `STRATCHEM`, `GMICHEM`, `GEOSCHEMCHEM` or `TR` into the internal bundle `CHEM_TRACERS`, with the child's name as a prefix (`GEOS_ChemGridComp.F90`, Init; `src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GridCompGetFriendlies0`). Increments are computed by differencing this bundle around each phase (`Shared/Chem_Shared/Bundle_IncrementMod.F90`). [Verified]
- Likely reason: Physics can list all chemistry tracers without knowing the child set. [Inferred]
- Question for maintainers: the fixed list of eight names in Init is not derived from the children; is that intended?

### Sibling coupling written as connections in the parent, guarded by enable flags
- Pattern: all chemistry-to-chemistry exchanges are `MAPL_AddConnectivity` calls in the parent's SetServices, each under the `enable_` flags of its two ends and sometimes the provider's name (`GEOS_ChemGridComp.F90`, SetServices). [Verified]
- Likely reason: children do not know about each other; the parent is the single place where the configuration is wired. [Inferred]

### Upstream codes reached through thin interface components
- Pattern: HEMCO and GEOS-Chem are compiled from nested upstream repositories by `CMakeLists.txt` files here that glob their sources; `HEMCO_GridComp` wraps HEMCO, and `GEOSCHEMchem_GridComp` copies the GEOS-Chem interface file `Interfaces/GCHP/Chem_GridCompMod.F90` into the build tree (`Shared/HEMCO/CMakeLists.txt`; `GEOSCHEMchem_GridComp/CMakeLists.txt`). [Verified]
- Likely reason: the upstream projects are developed elsewhere and used unchanged in several models. [Inferred]

### Child repositories behind a missing-directory stub
- Pattern: `CMakeLists.txt` lists sub-components, and a missing directory becomes a generated stub, which is how `H2O_GridComp`, absent from the repository, still builds (`CMakeLists.txt`, `alldirs`; fixture `@cmake/esma_support/esma_add_library.cmake`, `esma_create_stub_component`). [Verified]
- Likely reason: sparse or reduced checkouts can build without every child. [Inferred]
