# Critique Report: Crystal MCP Server (Non-MLFF)

## Scope and Constraints
- Reviewed source under `src/`, `docs/`, `tests/`, and top-level config/scripts.
- Excluded generated or third-party artifacts (`dist/`, `node_modules/`, `htmlcov/`, `__pycache__/`).
- MLFF-specific issues intentionally omitted per request.

## Critical Breakages (Tool Calls Fail or Crash)

1) Transformation tools send the wrong payload keys to Python.
- Evidence: TypeScript uses `structure` and `matrix` (`src/types/tools.ts:198-205`), but Python expects `structure_dict` and `scaling_matrix` (`src/python/structure_tools.py:166-174`).
- Impact: `structure_tools.py` raises `TypeError: unexpected keyword argument` on every call to supercell, slab, defect, alloy, heterostructure, adsorbate, strain.
- Fix: Normalize keys in TS before dispatch (`structure_dict`, `scaling_matrix`) and update schemas/docs to match.

2) Python entrypoints pass unfiltered kwargs into functions.
- Evidence: `structure_tools.py` and `space_group_scanner.py` call functions with `**params` directly (`src/python/structure_tools.py:872-888`, `src/python/space_group_scanner.py:206-228`).
- Impact: Extra keys like `operation`, `structure`, `matrix`, `fix_atoms`, `strain_type` trigger runtime exceptions.
- Fix: Filter allowed keys before calling, or accept `**kwargs` and ignore unknowns.

3) `generate_slab` schema includes `fix_atoms`, but Python does not accept it.
- Evidence: TS schema allows `fix_atoms` (`src/types/tools.ts:241-244`); Python signature lacks it (`src/python/structure_tools.py:309-320`).
- Impact: Any request with `fix_atoms` crashes the backend.
- Fix: Add `fix_atoms` support or remove from schema.

4) `apply_strain` can be called without `strain_tensor` but Python requires it.
- Evidence: TS makes `strain_tensor` optional (`src/types/tools.ts:639-647`); Python requires it (`src/python/structure_tools.py:548-558`).
- Impact: Valid TS inputs (strain_type/strain_value only) crash Python.
- Fix: Make `strain_tensor` required in schema or implement `strain_type` handling in Python.

5) `space_group_scanner.py` passes unsupported kwargs to `generate_from_spacegroup`.
- Evidence: `generate_from_spacegroup` accepts `factor`, not `volume_factor` or `max_attempts` (`src/python/generators/bulk/spacegroups.py:97-117`). `space_group_scanner.py` passes both (`src/python/space_group_scanner.py:72-112`).
- Impact: Space-group scan fails immediately with `TypeError`.
- Fix: Map `volume_factor` -> `factor`, remove `max_attempts` or implement in generator.

6) `formatStructureOutput` assumes fields that many tools do not provide.
- Evidence: It requires `space_group` and `metadata.density` (`src/utils/formatting.ts:16-47`), but outputs from `structure_tools.py` and `nanostructure_generator.py` omit them.
- Impact: Even successful Python calls crash while formatting.
- Fix: Make formatting resilient to missing fields or standardize outputs to include space-group/metadata everywhere.

7) `generate_space_group_scan` formatter expects schema not produced by backend.
- Evidence: Formatter expects `result.structure.space_group` (`src/utils/formatting.ts:115-118`) but backend structures from `generate_from_spacegroup` lack `space_group` in the structure dict (`src/python/generators/bulk/spacegroups.py:275-287`).
- Impact: Space-group scan responses throw when formatted.
- Fix: Include `space_group` in returned structure or adjust formatter to use top-level `spacegroup_number`/`spacegroup_symbol`.

8) Advanced tools are implemented but not exposed to clients.
- Evidence: Server handles `generate_prototype`, `generate_twisted_bilayer`, `generate_high_entropy_alloy`, `generate_2d_material`, `generate_mof`, `generate_cage` (`src/server.ts:38-75`), but `TOOL_DEFINITIONS` does not include them (`src/types/tools.ts:652-887`).
- Impact: Tools are undiscoverable via ListTools in standard MCP clients.
- Fix: Add schemas and metadata for these tools in `TOOL_DEFINITIONS`.

## High-Severity Issues (Incorrect or Misleading Results)

9) `generate_slab` interprets `thickness` as Angstroms, not layers.
- Evidence: TS describes `thickness` as number of layers (`src/types/tools.ts:229-233`), but Python uses it as `min_slab_size` in Angstroms (`src/python/structure_tools.py:333-446`).
- Impact: Slabs are the wrong thickness and fixed-layer counts are incorrect.
- Fix: Convert layers to Angstroms using interplanar spacing or rename the parameter.

10) `generate_molecular_crystal` accepts `lattice_params` but passes raw dict to PyXtal.
- Evidence: Python calls `crystal.from_random(..., lattice=lattice_params)` where `lattice_params` is a dict (`src/python/crystal_generator.py:546-585`).
- Impact: Likely runtime error or ignored lattice input.
- Fix: Convert to PyXtal lattice object or reject partial lattice params.

11) `space_group` accepts symbols in TS but Python expects ints.
- Evidence: TS schemas allow string symbols (`src/types/tools.ts:27-37`, `src/types/tools.ts:488-491`); Python passes `space_group` to `Group()` as an int (`src/python/crystal_generator.py:393`).
- Impact: Symbol inputs crash.
- Fix: Convert symbols to numbers or use symbol-aware APIs.

12) `comprehensive_generate` can crash when `list_available` is passed.
- Evidence: Python strips only `operation` and `category`, leaving `list_available` in params, and re-raises `TypeError` on unexpected kwargs (`src/python/comprehensive_structures.py:52-107`).
- Impact: Some requests crash with non-JSON output.
- Fix: Strip `list_available` in Python or TS before dispatch.

13) `comprehensive_generate` lacks schema for `operation_info`.
- Evidence: Python supports `operation_info` with `operation_name` (`src/python/comprehensive_structures.py:32-44`), but TS schema has no `operation_name` (`src/types/tools.ts:97-155`).
- Impact: Feature exists but is unreachable through validated MCP calls.
- Fix: Add `operation_name` to schema and docs.

14) `add_adsorbate` claims string molecule input is supported, but Python rejects it.
- Evidence: TS schema allows string (`src/types/tools.ts:620-631`); Python expects dict and errors on strings (`src/python/structure_tools.py:768-824`).
- Impact: Valid TS inputs fail in Python.
- Fix: Accept string names in Python (build molecule) or disallow in schema.

15) `make_supercell` schema claims `[nx, ny, nz]` but does not allow it.
- Evidence: Schema only permits `z.array(z.array(z.number()))` or string (`src/types/tools.ts:198-205`).
- Impact: Legit vector scaling is rejected at validation.
- Fix: Accept `z.array(z.number()).length(3)` in the schema.

## Medium-Severity Issues (Robustness, Consistency, UX)

16) Output schema inconsistency: `atoms` vs `sites` and missing `space_group`/`metadata`.
- Evidence: `structure_tools.py` produces `atoms` only; generator registry often uses both `sites` and `atoms` (`src/python/generators/bulk/spacegroups.py:80-87`).
- Impact: Clients cannot rely on a stable structure shape.
- Fix: Standardize a single output schema and enforce it across generators.

17) CIF export may be invalid if lattice parameters are missing.
- Evidence: `generateCIF` only writes cell lengths if `lattice.a` exists (`src/utils/structure-formats.ts:18-40`).
- Impact: Some structures with only `matrix` produce incomplete CIFs.
- Fix: Derive `a/b/c/alpha/beta/gamma` from `matrix` when missing.

18) Visualization HTML fails to load.
- Evidence: HTML uses `$()` but no jQuery is included (`src/python/visualization.py:53-74`).
- Impact: Browser errors, blank viewer.
- Fix: Replace with vanilla JS or include jQuery.

19) Visualization dependencies are not in `requirements.txt`.
- Evidence: `visualization.py` imports `matplotlib` but `requirements.txt` omits it.
- Impact: `generate_visualization` fails at import time.
- Fix: Add `matplotlib` to requirements or gate PNG output behind optional dependency.

20) Cross-platform issues in the Python bridge.
- Evidence: Uses `/tmp` for temp files and `python3` hardcoded (`src/utils/python-bridge.ts:366-370`, `src/utils/python-bridge.ts:417-434`).
- Impact: Fails on Windows by default.
- Fix: Use `os.tmpdir()` and allow configurable python executable.

21) `checkPythonAvailable` does not validate required packages.
- Evidence: Only runs `python3 --version` (`src/utils/python-bridge.ts:417-434`).
- Impact: Startup passes even when core dependencies are missing; failures show up later.
- Fix: Add a light dependency check (import pyxtal/spglib/ase) at startup.

22) `generate_crystal` ignores partial lattice params.
- Evidence: Only uses `lattice_params` if all six keys exist (`src/python/crystal_generator.py:357-383`).
- Impact: Users supply partial params that are silently ignored.
- Fix: Either reject partial inputs or fill missing values.

23) `space_group_scanner` ignores `num_atoms`.
- Evidence: Parameter exists but is never used (`src/python/space_group_scanner.py:20-52`, `src/python/space_group_scanner.py:206-228`).
- Impact: User intent is silently ignored.
- Fix: Map `num_atoms` into composition scaling or reject unsupported input.

24) `comprehensive_generate` skips input validation.
- Evidence: No `safeParse` in `src/tools/generation/comprehensive-generator.ts`.
- Impact: Bad inputs propagate to Python, increasing crash risk.
- Fix: Add zod validation like other tools.

## Documentation and Tests

25) README examples do not match current schemas.
- Evidence: Example uses `scaling_matrix` while schema expects `matrix` (README vs `src/types/tools.ts:198-205`).
- Impact: Copy-paste examples fail validation.
- Fix: Align docs with actual schema field names.

26) Tests are out of sync with actual outputs.
- Evidence: Jest expectations refer to old output headers (see `tests/test_mcp_features.test.ts` vs `src/utils/formatting.ts:16-47`).
- Impact: Tests fail even when tool logic is corrected.
- Fix: Update expected strings to match current formatter outputs.

## Notes
- MLFF-specific issues are intentionally excluded from this report.
- Several backend files still declare "DEPRECATED" but are actively used, which is a maintainability and signaling problem.

