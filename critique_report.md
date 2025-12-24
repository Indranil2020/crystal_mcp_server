# Critique Report: Crystal MCP Server

## Scope and Method
- Reviewed source under `src/`, `docs/`, `tests/`, and top-level config/scripts.
- Excluded generated or third-party artifacts (`dist/`, `node_modules/`, `htmlcov/`, `__pycache__/`) from code inspection.
- Focused on MCP server correctness, protocol integration, Python bridge reliability, and backend/tool contract alignment.

## Critical Issues (Breaks Core Functionality)

1) Tool handlers read Python results from the wrong path.
- Evidence: `executePythonWithJSON` returns parsed JSON at `result.data`, but handlers dereference `result.data.data` (example: `src/tools/generation/generate-crystal.ts:51`, `src/tools/optimization/calculate-energy-mlff.ts:35`).
- Impact: `pythonResult` becomes `undefined`, and most tools crash when checking `pythonResult.success`.
- Fix: Replace all `result.data.data` reads with `result.data`. Add unit tests for each tool to verify a success response path.

2) Python entrypoints pass unfiltered kwargs, causing `TypeError` for unexpected keys.
- Evidence: `structure_tools.py` calls `make_supercell(**params)` etc while `params` includes `operation` and mismatched keys (`src/python/structure_tools.py:875`); MLFF and validators do the same (`src/python/mlff_calculator.py:567`, `src/python/validators.py:599`).
- Impact: Most Python tools crash on first call because of unexpected keyword arguments.
- Fix: Filter params before calling (`{k: v for k, v in params.items() if k in allowed_keys}`), or change Python signatures to accept `**kwargs` and ignore unknowns.

3) Transformation tool payload keys do not match Python expectations.
- Evidence: TypeScript sends `structure` and `matrix` (`src/tools/transformation/supercell.ts:23-33`), but Python expects `structure_dict` and `scaling_matrix` (`src/python/structure_tools.py:166-174`).
- Impact: Every transformation tool is broken even if the `operation` key is filtered.
- Fix: Normalize keys in TypeScript before dispatch (`structure_dict` + `scaling_matrix`), and align schemas to match Python.

4) `generate_space_group_scan` references a missing backend script.
- Evidence: `src/tools/generation/space-group-scan.ts:32` calls `space_group_scanner.py`, which does not exist in `src/python/`.
- Impact: The tool is dead on arrival.
- Fix: Implement `space_group_scanner.py` or route to an existing generator in `comprehensive_structures.py`.

5) Ground state search input/output contracts are completely mismatched.
- Evidence: TS expects `optimization_settings` and formats output via `formatGroundStateSearchOutput` (`src/tools/optimization/ground-state-search.ts:21-58`, `src/utils/formatting.ts:130-210`), but Python returns `best_structure` and does not accept those inputs (`src/python/ground_state_searcher.py:18-122`).
- Impact: The tool either fails or returns data that the formatter cannot consume.
- Fix: Align TS schema and formatter to the Python output, or update Python to return the richer schema used by TS.

6) Formatting code assumes fields that do not exist in most structures.
- Evidence: `formatStructureOutput` requires `space_group` and `metadata.density` (`src/utils/formatting.ts:16-47`). Outputs from structure tools, nanostructures, and MLFF lack these fields (e.g., `src/python/structure_tools.py:109-152`, `src/python/nanostructure_generator.py:215-260`, `src/python/mlff_calculator.py:232-270`).
- Impact: Many successful tool calls crash during formatting.
- Fix: Make formatting resilient to missing fields or enforce a consistent structure schema across all Python outputs.

7) MLFF energy tool uses the wrong operation name.
- Evidence: TS sends `operation: "calculate_energy"` (`src/tools/optimization/calculate-energy-mlff.ts:23`) but Python only accepts `"energy"` (`src/python/mlff_calculator.py:565-570`).
- Impact: Tool always returns "Unknown operation" even if all dependencies are installed.
- Fix: Change TS to send `"energy"` or update Python to accept `"calculate_energy"`.

8) Validation tool sends extra `structure` key, breaking Python.
- Evidence: TS passes `{ structure_dict, ...parsed.data }` (`src/tools/analysis/validation.ts:26-33`); Python uses `validate_structure(**params)` (`src/python/validators.py:599`).
- Impact: Python raises `TypeError: unexpected keyword argument 'structure'`.
- Fix: Only send `structure_dict` (and relevant check params) from TS.

## High-Severity Issues (Crashes or Severe Misbehavior)

9) JSON parsing in the Python bridge can crash the Node process.
- Evidence: `JSON.parse` is called without a try/catch (`src/utils/python-bridge.ts:150-171`).
- Impact: Any stray stdout line (warnings, logs) can terminate the MCP server.
- Fix: Wrap `JSON.parse` in try/catch and return a structured error.

10) Python timeouts do not terminate the subprocess.
- Evidence: The timeout flag is set but the Python process is never killed (`src/utils/python-bridge.ts:224-270`).
- Impact: Long-running jobs leak processes and tie up CPU/memory even after timeout.
- Fix: Call `shell.terminate()` or `childProcess.kill()` when timing out.

11) File I/O helpers claim defensive behavior but throw synchronously.
- Evidence: `fs.accessSync`, `fs.writeFileSync`, and `JSON.parse` are used without guards (`src/utils/file-io.ts:54-71`, `src/utils/file-io.ts:198-214`).
- Impact: A single permission error or malformed file can crash the server.
- Fix: Wrap these in try/catch and return `Result` failures instead.

12) `apply_strain` schema and Python implementation disagree on tensor format.
- Evidence: TS expects flattened length-9 array (`src/types/tools.ts:639-647`), Python expects a 3x3 matrix (`src/python/structure_tools.py:548-585`).
- Impact: Valid TS inputs crash or produce wrong deformation.
- Fix: Accept flattened arrays in Python or reshape in TS before dispatch.

13) `validate_structure` checks do not match schema defaults.
- Evidence: TS default checks include `symmetry` and `overlaps` (`src/types/tools.ts:316-327`), but Python only understands `distances`, `lattice`, and `stoichiometry` (`src/python/validators.py:247-306`).
- Impact: Requested checks are ignored silently, producing incomplete validation.
- Fix: Align check names or implement the missing checks in Python.

14) Backends accept string space-group symbols but then treat them as ints.
- Evidence: TS schemas allow string symbols (`src/types/tools.ts:27-37`, `src/types/tools.ts:488-491`); Python calls `Group(space_group)` expecting an int (`src/python/crystal_generator.py:393`).
- Impact: Symbol inputs crash at runtime.
- Fix: Convert symbols to numbers or use a library that accepts symbols.

15) `generate_molecular_crystal` ignores `min_distance` and mishandles `lattice_params`.
- Evidence: `min_distance` is never used and `lattice_params` is passed directly to `pyxtal` without conversion (`src/python/crystal_generator.py:540-585`).
- Impact: Inputs appear accepted but have no effect, eroding trust in tool output.
- Fix: Apply `min_distance` as a constraint or validation gate; convert lattice params to a PyXtal lattice object as done in `generate_crystal`.

16) Export tool advertises formats it cannot generate.
- Evidence: Schema includes `contcar`, `extxyz`, `pdb`, `pwscf`, `castep`, `lammps` (`src/types/tools.ts:402-419`), but implementation only outputs `cif`, `poscar`, `xyz`, `json` (`src/tools/export/export-structure.ts:28-53`).
- Impact: Requests for advertised formats silently drop outputs.
- Fix: Implement missing formats or remove them from the schema.

17) Visualization HTML depends on jQuery but does not load it.
- Evidence: HTML calls `$()` without including jQuery (`src/python/visualization.py:53-74`).
- Impact: Generated HTML errors on load and renders nothing.
- Fix: Use plain JS (`document.addEventListener`) or include jQuery explicitly.

## Medium-Severity Issues (Incorrect Results or Hidden Drift)

18) Symmetry analysis uses incomplete element-to-Z mappings.
- Evidence: Unknown elements default to Z=1 (H) in `dict_to_spglib_cell` (`src/python/symmetry_analyzer.py:71-92`) and reverse mapping is limited (`src/python/symmetry_analyzer.py:284-304`).
- Impact: Symmetry results can be wrong for most transition/metals and heavier elements.
- Fix: Use a full periodic table mapping (e.g., from pymatgen).

19) Density and packing calculations use truncated element data.
- Evidence: Atomic masses and radii only cover elements up to Zn (`src/python/validators.py:517-575`).
- Impact: Density and packing fraction are incorrect for many materials.
- Fix: Use a complete dataset or rely on pymatgen for these properties.

20) MLFF loaders assume packages exist and hard-code CPU.
- Evidence: Imports are unconditional and device is fixed to CPU (`src/python/mlff_calculator.py:52-120`).
- Impact: Missing packages crash the tool; GPU acceleration is unavailable even when supported.
- Fix: Guard imports, surface dependency errors, and allow device selection.

21) Schema allows `GPMin` optimizer but implementation does not support it.
- Evidence: `GPMin` is valid in TS (`src/types/tools.ts:262-269`) but not implemented in Python (`src/python/mlff_calculator.py:404-420`).
- Impact: Requests validated by TS fail in Python.
- Fix: Remove from schema or implement the optimizer.

22) MLFF stress output shape does not match TypeScript type.
- Evidence: Python returns Voigt 6-vector (`src/python/mlff_calculator.py:537-540`), but TS expects 3x3 matrix (`src/types/crystal.ts:148-159`).
- Impact: Clients expecting a matrix will misinterpret stress values.
- Fix: Convert to a 3x3 matrix or update TS types to match.

23) Slab thickness is interpreted as Angstroms, not layers.
- Evidence: TS describes `thickness` as number of layers (`src/types/tools.ts:229-233`), Python passes it to `SlabGenerator(min_slab_size=...)` in Angstroms (`src/python/structure_tools.py:333-446`).
- Impact: Users get much thinner or thicker slabs than requested.
- Fix: Convert layers to Angstroms using interplanar spacing, or rename the parameter.

24) Many schemas accept file paths for `structure`, but backends only accept dicts.
- Evidence: TS allows `string` in `structure` across multiple schemas (`src/types/tools.ts:198-205` and similar), but Python `dict_to_structure` only accepts dicts (`src/python/structure_tools.py:64-100`).
- Impact: Passing a file path consistently fails with "invalid structure".
- Fix: Implement file loading in TS or Python, or remove `string` from schemas.

25) Tool metadata annotations are defined but discarded.
- Evidence: `TOOL_DEFINITIONS` includes `annotations` (`src/types/tools.ts:437-442`) but `ListTools` ignores them (`src/server.ts:68-73`).
- Impact: Clients lose read-only/destructive hints for safe tool selection.
- Fix: Include `annotations` in the `ListTools` response.

26) Version reporting is inconsistent.
- Evidence: Server reports `2.0.0` (`src/server.ts:55-59`), but package.json is `1.0.0`.
- Impact: Clients and logs cannot trust the version.
- Fix: Use a single source of truth (package.json).

27) Python temp file path and Python executable are not cross-platform.
- Evidence: Uses `/tmp` and `python3` by default (`src/utils/python-bridge.ts:366-370`, `src/utils/python-bridge.ts:417-434`).
- Impact: Windows users fail out of the box.
- Fix: Use `os.tmpdir()` and allow `python` on Windows via config/env.

## Documentation and Testing Gaps

28) README examples do not match current schemas and implementations.
- Example: `make_supercell` shows `scaling_matrix` but schema expects `matrix` (`README.md:64-80`, `src/types/tools.ts:198-205`).
- Example: `generate_space_group_scan` is documented but has no backend script.
- Fix: Update docs to match actual tool names and parameter keys.

29) Tests are misaligned with actual output formats and tool behavior.
- Evidence: Tests expect strings that no longer exist in formatter output (`tests/test_mcp_features.test.ts:11-40`, `src/utils/formatting.ts:16-47`).
- Impact: CI will fail even when tools work.
- Fix: Update tests after correcting formatter and tool outputs.

30) Deprecated Python modules are still the primary backend for many tools.
- Evidence: `crystal_generator.py`, `nanostructure_generator.py`, `molecule_generator.py`, `ground_state_searcher.py` are marked DEPRECATED but used by TS.
- Impact: Architectural drift and duplicated logic between legacy and modular generators.
- Fix: Migrate TS tools to the modular generator registry or remove deprecation.

## Maintainability and Hygiene

31) Unused modules and configs add maintenance noise.
- Evidence: `src/utils/validation.ts` and `src/config/*.json` are unused.
- Impact: Harder audits and higher cognitive load.
- Fix: Remove or integrate these artifacts.

32) Duplicate export logic exists in multiple places.
- Evidence: CIF/POSCAR/XYZ generation exists in both `generate-crystal.ts` and `export-structure.ts`.
- Impact: Inconsistent output and duplicated bugs.
- Fix: Consolidate into shared helpers.

33) Nonstandard units and labels reduce scientific clarity.
- Evidence: Volume displayed with "u" glyphs instead of "A^3" (`src/utils/formatting.ts:31-44`).
- Impact: Confusing for scientific users and inconsistent with most tooling.
- Fix: Use explicit ASCII units (e.g., "A^3").

34) Repo hygiene issues likely pollute version control.
- Evidence: `.pytest_cache` and `__pycache__` exist in `src/python/`.
- Impact: Noisy diffs and larger repository size.
- Fix: Add to `.gitignore` and remove from repo.

