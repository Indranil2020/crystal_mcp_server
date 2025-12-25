# Critique Fix Task List

This checklist enumerates every unresolved or partially resolved item from:
- `docs/expert_critique_report.md`
- `critique_report.md`
- `docs/critique_status_report.md`

Status will be updated as fixes land.

## Protocol and Server
- [ ] Expose tool annotations in ListTools response
- [ ] Remove unsafe schema casts in ListTools
- [ ] Use package.json version in server and startup logs

## Type Schemas and Validation
- [ ] Add Matrix3x3 schema and use it for scaling/strain inputs
- [ ] Convert nanostructure schema to discriminated union with type-specific params
- [ ] Enforce strain tensor shape + require strain_tensor or (strain_type + strain_value)
- [ ] Validate min_distance keys with element-pair regex
- [ ] Add shared base schemas for structure/space-group/composition inputs

## Crystal Generator and Crystallography Data
- [ ] Preserve Wyckoff letter/multiplicity/site symmetry per atom
- [ ] Validate user-provided Wyckoff positions against space group
- [ ] Compute and populate reciprocal lattice matrix
- [ ] Normalize space-group symbols to numbers (string support)
- [ ] Improve element validation using full Element table (remove hardcoded list)
- [ ] Decide behavior for partial lattice_params (reject or complete)

## Structure Outputs and Consistency
- [ ] Standardize atoms/sites output (provide consistent alias)
- [ ] Provide space_group info consistently in structure_tools outputs
- [ ] Fix CIF export when only matrix is present

## Scientific Accuracy Improvements
- [ ] Improve heterostructure stacking (lattice matching or clearer model)
- [ ] Improve adsorbate positioning (rotation sampling / contact checks)
- [ ] Improve interstitial placement (site search or better heuristics)
- [ ] Add strain-size warnings for large strains
- [ ] Document or refine packing fraction radii selection
- [ ] Document defect formation energy sources or allow overrides

## Export and Formatting
- [ ] Align export_structure schema with supported formats
- [ ] Error on unsupported export formats instead of silent skip
- [ ] Remove non-ASCII unit symbols in outputs

## Validation and Error Handling
- [ ] Expand validate_structure to handle schema checks (symmetry/overlaps/cell/charge)
- [ ] Add Python response shape validation in python-bridge
- [ ] Use MODEL_NOT_AVAILABLE error path where applicable

## Tooling and MLFF
- [ ] Implement GPMin optimizer support or remove from schema
- [ ] Decide on constrain_symmetry support (implement or remove)

## Documentation
- [ ] Fix outdated API reference parameter names
- [ ] Add Wyckoff primer and symmetry tolerance explanation
- [ ] Document strain tensor conventions
- [ ] Add crystallographic references and unit conventions

## Tests
- [ ] Fix tests using old parameter names
- [ ] Add tests for advanced tools
- [ ] Add additional edge-case tests (invalid inputs)

## Misc
- [ ] Fix bulk module import failure for missing exports
- [ ] Cross-platform python path default (avoid hardcoded python3)
- [ ] Add dependency sanity check in checkPythonAvailable
- [ ] Fix molecule generator behavior for unknown molecules
