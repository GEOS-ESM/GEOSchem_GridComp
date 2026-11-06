# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSchem_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSchem_GridComp` @ `3c925f2` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; key `GEOSchem_GridComp`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** only `GOCART` is built; the Chemistry
> gridded component is present but not compiled. On `release/MAPL-v3`,
> `CMakeLists.txt` calls `esma_add_subdirectories(GOCART)` ("Only GOCART is needed
> on release/MAPL-v3 (for GOCART2G_GridComp)"; "Restore subdirs as components are
> ported to MAPL3."); the line adding `Shared GOCART TR GMI StratChem MAM MATRIX
> CARMA GAAS ACHEM` is commented out. The `esma_add_library` block that builds
> `GEOS_ChemGridComp.F90` and `GEOS_ChemEnvGridComp.F90` is commented out too
> ("GEOS_ChemGridComp.F90 USE-s non-ported components"), with the in-repository
> children `GEOSpchem_GridComp`, `GEOSCHEMchem_GridComp`, `H2O_GridComp`,
> `DNA_GridComp`, `HEMCO_GridComp` and `RRG` ("not yet ported to MAPL3"). Not
> built; not yet ported to MAPL3: the Chemistry parent and `CHEMENV`, the
> in-repository children, and `Shared` (`Chem_Base`, `Chem_Shared`, `HEMCO`). The
> code of these is the MAPL2 code, kept for a later port; statements about
> `GEOS_ChemGridComp.F90` describe the file as it is. Of the nested child
> repositories only `@GOCART` is added to the build here; the others are listed in
> the commented-out line.

## What this is

The Chemistry and aerosol component of the GEOS atmosphere. The Physics
composite adds it as child `CHEMISTRY` and runs it in two phases
(`../GEOS_PhysicsGridComp.F90`, SetServices, Run). The repository holds the
Chemistry parent (`GEOS_ChemGridComp.F90`), the always-present environment
child `CHEMENV` (`GEOS_ChemEnvGridComp.F90`), four in-repository children
(`GEOSpchem_GridComp`, `HEMCO_GridComp`, `DNA_GridComp`,
`GEOSCHEMchem_GridComp`), shared libraries (`Shared/`) and build glue. Most
chemistry and aerosol children are separate nested repositories
(`@GOCART`, `@ACHEM`, `@GMI`, `@StratChem`, `@CARMA`, `@MAM`, `@MATRIX`,
`@GAAS`, `@TR`, `@RRG`, `@QuickChem`, and the upstream codes `@HEMCO` and
`@geos-chem`): edit them only in their own repository.

## Documents

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Module map, the Chemistry parent (children and `ENABLE_*` keys, run phases, providers of radiatively active gases and aerosols, sibling connections, tracer and increment bundles), `CHEMENV`, build, traps, observations. Read first. |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| `@GOCART/docs/ARCHITECTURE.md` | The GOCART repository (aerosol model) |
| `@ACHEM/docs/ARCHITECTURE.md` | The ACHEM repository (gas and aqueous phase chemistry) |
| `@CARMA/docs/ARCHITECTURE.md` | The CARMA repository (aerosol and cloud microphysics) |
| `@GAAS/docs/ARCHITECTURE.md` | The GAAS repository (aerosol analysis system) |
| `@GMI/docs/ARCHITECTURE.md` | The GMI repository (combined troposphere/stratosphere chemistry) |
| `@MAM/docs/ARCHITECTURE.md` | The MAM repository (modal aerosol microphysics) |
| `@MATRIX/docs/ARCHITECTURE.md` | The MATRIX repository (aerosol microphysics) |
| `@QuickChem/docs/ARCHITECTURE.md` | The QuickChem repository (quick chemistry (not wired in)) |
| `@RRG/docs/ARCHITECTURE.md` | The RRG repository (carbon gases production and loss) |
| `@StratChem/docs/ARCHITECTURE.md` | The StratChem repository (stratospheric chemistry) |
| `@TR/docs/ARCHITECTURE.md` | The TR repository (passive tracers) |
| `../ARCHITECTURE.md` | The Physics composite that runs Chemistry |
| `../../../docs/ARCHITECTURE.md` | The `GEOSgcm_GridComp` repository |

## Rules

- Nested child repositories are gitignored here (`.gitignore`); `git status`
  in this repository does not show their changes. Use
  `git -C <dir> rev-parse --show-toplevel` before editing or committing.
- A pull request must update `CHANGELOG.md` unless a skip label is set
  (`.github/workflows/changelog-enforcer.yml`). Default code owners:
  chemistry gatekeepers (`.github/CODEOWNERS`).
- Note bugs and surprises under "Observations" in `docs/ARCHITECTURE.md`; do not
  fix them unless asked. Keep the documents current when a child, an
  `ENABLE_*` key, a provider key or a sibling connection changes. Cite files
  and routines, not line numbers.
- Adding a child to the Chemistry parent touches several places: the `use`
  line, the `enable_` flag and `ENABLE_` key, the `MAPL_AddChild` line,
  `GetProvider_`, the sibling connections, `CMakeLists.txt`,
  `GEOS_ChemGridComp.rc` and the `GEOSgcm_App` child lists.

## Build and test

- Build status on this line: see the note at the top. Only `GOCART` is built
  (`CMakeLists.txt`); `GEOS_ChemGridComp.F90` is not compiled, so the Chemistry
  component is not part of any executable built from this line. The block that
  builds it, when restored, calls the fixture's `@cmake` macros
  (`esma_set_this`, `esma_add_library`) and so builds only inside the GEOSgcm
  fixture. See the fixture `AGENTS.md` (`../../../../../../AGENTS.md` in the fixture).
- No tests are registered in this repository's own `CMakeLists.txt` files.
  CI builds the fixture with several compilers
  (`.github/workflows/workflow.yml`) and runs short GCM runs
  (`.circleci/config.yml`).
