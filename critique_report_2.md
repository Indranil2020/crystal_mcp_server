CRITICAL FAILURES
- Zeolite generation builds random Si/Al/O/cation positions with no framework connectivity or collision checks, so the output is not a zeolite structure and is not scientifically valid for DFT/MD. src/python/generators/bulk/zeolites.py:135
- generate_slab_supercell discards the input structure geometry and composition and builds a synthetic slab from the first element only, producing wrong slabs for any compound. src/python/generators/bulk/supercells.py:263
- generate_double_perovskite accepts ordering="random" but returns a simplified 7-atom cell that violates the A2BB'X6 stoichiometry without warning. src/python/generators/bulk/perovskites.py:152
- LAMMPS export always writes an orthogonal box and omits triclinic tilt factors, so non-orthogonal cells are exported with the wrong simulation cell. src/python/generators/output_formats/converters.py:328

HIGH-RISK DESIGN ISSUES
- generate_from_spacegroup applies requested a/b/c even when they violate the crystal system, while top-level spacegroup metadata remains the original PyXtal group, so reported symmetry can be inconsistent with the returned structure. src/python/generators/bulk/spacegroups.py:335
- generate_crystal silently perturbs stoichiometry when num_atoms does not match the input ratios (first species is adjusted), with no warning to the caller. src/python/crystal_generator.py:520
- Workflow editing uses species/positions/cell in Cartesian Angstroms with no adapter to the atoms/lattice fractional schema used by other tools, so edited structures are easily misinterpreted downstream. src/python/generators/workflow/iterative_editing.py:168
- 2D vacuum semantics are inconsistent (c = vacuum in TMD/Xene generators) and 2D utilities update lattice/coords without recomputing cartesian, which can yield insufficient vacuum and stale cartesian coordinates. src/python/generators/two_d/tmds.py:109 src/python/generators/two_d/xenes.py:150 src/python/generators/two_d/base.py:149
- Adsorbate placement and heterostructure stacking assume a global z axis and reconstruct the lattice from lengths/angles, which misorients structures for tilted or non-orthogonal cells. src/python/structure_tools.py:880 src/python/structure_tools.py:1026
- Cell setting transformations use incomplete element maps and emit only fractional coords, so heavy elements become Z=0/X and downstream exporters/visualizers can use wrong coordinates. src/python/generators/bulk/cell_settings.py:97

MEDIUM / LOW-RISK IMPROVEMENTS
- VASP export writes fractional coords even when direct=False, so Cartesian POSCAR output is wrong. src/python/generators/output_formats/converters.py:156
- PDB export assumes orthogonal cells and uses a/b/c scaling for coordinates, which is wrong for non-orthogonal lattices. src/python/generators/output_formats/converters.py:548
- TS POSCAR/XYZ exporters depend on lattice.matrix/cartesian; if missing, POSCAR lacks lattice vectors and XYZ falls back to fractional coords as Cartesian. src/utils/structure-formats.ts:121 src/utils/structure-formats.ts:172
- export_structure ignores include_metadata and dft_software options, so requested metadata/format behavior is silently dropped. src/types/tools.ts:424 src/tools/export/export-structure.ts:24
- StructureInputSchema accepts string|any but no file loading is implemented; passing a file path succeeds schema validation and then fails in Python. src/types/tools.ts:34
- Some generator parameters are accepted but unused (perovskite tilting, molecule formula), which misleads users about applied physics. src/python/generators/bulk/perovskites.py:47 src/python/molecule_generator.py:21

ACTIONABLE FIXES
- Mark toy generators (zeolite, slab_supercell, perovskite ordering=random) as non-scientific or require an explicit allow_toy=true flag, and add validation or remove them from production-facing tools. src/python/generators/bulk/zeolites.py:135 src/python/generators/bulk/supercells.py:245 src/python/generators/bulk/perovskites.py:152
- Enforce stoichiometry consistency in generate_crystal and warn or fail when num_atoms forces a composition change; add an explicit allow_stoichiometry_adjustment flag. src/python/crystal_generator.py:520
- Standardize the structure schema across tools (atoms/lattice with fractional + cartesian) and add adapters for workflow editing and 2D utilities; always recompute cartesian after lattice or coord edits. src/python/generators/workflow/iterative_editing.py:168 src/python/generators/two_d/base.py:149
- Fix exporters: add LAMMPS triclinic tilt factors, convert coordinates correctly for VASP Cartesian mode, and compute cartesian for XYZ/PDB from the lattice matrix when missing. src/python/generators/output_formats/converters.py:103 src/utils/structure-formats.ts:165
- Normalize vacuum semantics across all 2D generators (use c = thickness + vacuum), and add tests that verify minimum vacuum separation in the generated structures. src/python/generators/two_d/tmds.py:109 src/python/generators/two_d/xenes.py:150
- After lattice overrides in spacegroup generation, recompute symmetry and surface it in the top-level response (or reject incompatible lattice parameters). src/python/generators/bulk/spacegroups.py:335

SCIENTIFIC RISK ASSESSMENT
- Zeolite, slab_supercell, and perovskite ordering=random outputs can invalidate DFT/TDDFT/MD studies due to wrong composition or non-physical structures. src/python/generators/bulk/zeolites.py:135 src/python/generators/bulk/supercells.py:263 src/python/generators/bulk/perovskites.py:152
- Non-orthogonal cells exported to LAMMPS or POSCAR/XYZ without proper cell handling can corrupt forces, stresses, and energies. src/python/generators/output_formats/converters.py:328 src/utils/structure-formats.ts:165
- 2D structures may have insufficient vacuum or stale cartesian coordinates, leading to periodic image interactions or wrong geometries in downstream simulations. src/python/generators/two_d/tmds.py:109 src/python/generators/two_d/base.py:149
- Workflow-edited structures can be misread as fractional when they are Cartesian, propagating incorrect geometries into subsequent transformations and exports. src/python/generators/workflow/iterative_editing.py:168
- Adsorbate/heterostructure tools can misplace atoms on tilted or non-orthogonal surfaces, invalidating surface energies and adsorption geometries. src/python/structure_tools.py:880 src/python/structure_tools.py:1026

OPEN QUESTIONS / ASSUMPTIONS
- Should toy/placeholder generators be removed from the public tool registry or gated behind an explicit allow_toy flag?
- Do you want spacegroup generation to reject incompatible lattice parameters, or to return a downgraded symmetry with a warning?
- Should workflow editing adopt the core atoms/lattice schema, or should explicit adapters be introduced for both directions?
- Should vacuum be defined as spacing between periodic images for all 2D generators (c = thickness + vacuum) to ensure consistent DFT usage?
