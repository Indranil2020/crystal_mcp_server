#!/usr/bin/env python3
"""
Molecule Editor - MCP Entry Point

Enables iterative molecular modification through natural language or structured operations.
Integrates with existing build_molecule and build_molecular_cluster tools.

Input: JSON file path via argv[1] containing:
  - molecule: {smiles: str, annotations?: object}
  - operations: str (NL) or array of structured operations
  - validate?: bool (default: true)
  - optimize_geometry?: bool (default: true)
  - return_3d?: bool (default: true)

Output: JSON to stdout containing:
  - success: bool
  - structure: {lattice, atoms, metadata}
  - annotations: {stereocenters, rings, functional_groups, atom_labels, edit_history}
  - edit_result: {operations_applied, validation}
"""

import sys
import json
import os
import importlib.util
from typing import Dict, Any, List, Optional, Union
import numpy as np

# Ensure src/python is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Check module availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None
EDITING_AVAILABLE = importlib.util.find_spec("generators.molecule.editing") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem  # noqa: F401 - used for type checking

if EDITING_AVAILABLE:
    from generators.molecule.editing import (
        AnnotatedMolecule,
        MolecularEditor,
        EditOperation,
        EditOperationType,
        parse_nl_instruction,
        ValidationPipeline,
        ValidationResult,
    )


def create_error_response(code: str, message: str, details: Optional[Dict] = None) -> Dict[str, Any]:
    """Create a standardized error response."""
    response = {
        "success": False,
        "error": {
            "code": code,
            "message": message,
        }
    }
    if details:
        response["error"]["details"] = details
    return response


def parse_input_molecule(mol_data: Dict[str, Any]) -> Optional["AnnotatedMolecule"]:
    """
    Parse input molecule from request data.

    Args:
        mol_data: Dictionary with 'smiles' and optional 'annotations'

    Returns:
        AnnotatedMolecule or None if parsing fails
    """
    if not EDITING_AVAILABLE:
        return None

    smiles = mol_data.get("smiles")
    if not smiles:
        return None

    annotations = mol_data.get("annotations")

    # If annotations provided, reconstruct from dict
    if annotations:
        full_data = {"smiles": smiles, **annotations}
        return AnnotatedMolecule.from_dict(full_data)

    # Otherwise create from SMILES with auto-annotation
    return AnnotatedMolecule.from_smiles(smiles)


def parse_operations(
    operations: Union[str, List[Dict]],
    molecule: "AnnotatedMolecule"
) -> List["EditOperation"]:
    """
    Parse operations from either natural language or structured format.

    Args:
        operations: Either NL string or list of operation dicts
        molecule: The annotated molecule (for context in NL parsing)

    Returns:
        List of EditOperation objects
    """
    if not EDITING_AVAILABLE:
        return []

    # Natural language input
    if isinstance(operations, str):
        return parse_nl_instruction(operations, molecule)

    # Structured operations
    parsed_ops = []
    for op_dict in operations:
        if isinstance(op_dict, dict):
            op = EditOperation.from_dict(op_dict)
            if op is not None:
                parsed_ops.append(op)

    return parsed_ops


def molecule_to_structure(
    mol: "AnnotatedMolecule",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Convert AnnotatedMolecule to MCP structure format.

    Args:
        mol: The annotated molecule
        vacuum: Vacuum padding in Angstroms

    Returns:
        Structure dictionary compatible with MCP response format
    """
    coords = mol.coordinates
    atoms = []

    if mol._mol is not None and RDKIT_AVAILABLE:
        for atom in mol._mol.GetAtoms():
            atoms.append(atom.GetSymbol())
    else:
        # Fallback: parse from formula
        atoms = ["X"] * len(coords) if coords else []

    n_atoms = len(atoms)

    if n_atoms == 0 or not coords:
        return {}

    coords_array = np.array(coords)
    min_coords = coords_array.min(axis=0)
    max_coords = coords_array.max(axis=0)
    box_size = max_coords - min_coords + 2 * vacuum

    # Center molecule in box
    center_offset = box_size / 2 - (max_coords + min_coords) / 2
    centered_coords = coords_array + center_offset

    # Create cell matrix
    cell = np.diag(box_size)

    # Create sites
    sites = []
    for atom, coord in zip(atoms, centered_coords):
        frac_coord = coord / box_size
        sites.append({
            "element": atom,
            "coords": frac_coord.tolist(),
            "cartesian": coord.tolist(),
            "species": [{"element": atom, "occupation": 1.0}]
        })

    return {
        "lattice": {
            "matrix": cell.tolist(),
            "volume": float(np.prod(box_size)),
            "a": float(box_size[0]),
            "b": float(box_size[1]),
            "c": float(box_size[2]),
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0
        },
        "atoms": sites,
        "sites": sites,
        "space_group": {
            "number": 1,
            "symbol": "P1",
            "crystal_system": "triclinic"
        },
        "metadata": {
            "formula": mol.formula,
            "natoms": n_atoms,
            "smiles": mol.smiles,
            "canonical_smiles": mol.canonical_smiles,
            "molecular_weight": mol.molecular_weight,
            "pbc": [False, False, False],
            "source": "edit",
        }
    }


def run_validation(mol: "AnnotatedMolecule") -> "ValidationResult":
    """
    Run validation pipeline on edited molecule.

    Args:
        mol: The annotated molecule to validate

    Returns:
        ValidationResult
    """
    if not EDITING_AVAILABLE:
        return None

    pipeline = ValidationPipeline()
    return pipeline.validate(mol._mol)


def edit_molecule(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Main editing function.

    Args:
        params: Request parameters

    Returns:
        Response dictionary
    """
    # Check dependencies
    if not RDKIT_AVAILABLE:
        return create_error_response(
            "DEPENDENCY_MISSING",
            "RDKit is required for molecular editing",
            {"missing": "rdkit"}
        )

    if not EDITING_AVAILABLE:
        return create_error_response(
            "DEPENDENCY_MISSING",
            "Editing module not available",
            {"missing": "generators.molecule.editing"}
        )

    # Extract parameters
    mol_data = params.get("molecule")
    operations = params.get("operations")
    should_validate = params.get("validate", True)
    optimize_geometry = params.get("optimize_geometry", True)
    return_3d = params.get("return_3d", True)
    vacuum = float(params.get("vacuum", 10.0))

    # Validate required inputs
    if not mol_data:
        return create_error_response(
            "INVALID_INPUT",
            "Missing 'molecule' parameter",
            {"required": ["molecule", "operations"]}
        )

    if not operations:
        return create_error_response(
            "INVALID_INPUT",
            "Missing 'operations' parameter",
            {"required": ["molecule", "operations"]}
        )

    # Parse input molecule
    molecule = parse_input_molecule(mol_data)
    if molecule is None:
        return create_error_response(
            "INVALID_MOLECULE",
            f"Failed to parse molecule from SMILES: {mol_data.get('smiles')}",
            {"smiles": mol_data.get("smiles")}
        )

    # Parse operations
    parsed_ops = parse_operations(operations, molecule)

    if not parsed_ops:
        # Provide helpful feedback for NL parsing failures
        if isinstance(operations, str):
            return create_error_response(
                "PARSE_FAILED",
                f"Could not parse natural language instruction: '{operations}'",
                {
                    "instruction": operations,
                    "suggestions": [
                        "Use simpler verbs like 'hydroxylate', 'methylate', 'epimerize'",
                        "Specify target clearly: 'at C-7', 'on the benzene ring'",
                        "Use structured operations for complex edits",
                    ],
                    "supported_operations": [op.value for op in EditOperationType],
                }
            )
        else:
            return create_error_response(
                "INVALID_OPERATIONS",
                "No valid operations found in input",
                {"input": operations}
            )

    # Create editor and apply operations
    editor = MolecularEditor(molecule)
    results = editor.apply_operations(parsed_ops)

    # Check if any operations succeeded
    successful_ops = [r for r in results if r.success]
    failed_ops = [r for r in results if not r.success]

    if not successful_ops:
        return create_error_response(
            "ALL_OPERATIONS_FAILED",
            "All edit operations failed",
            {
                "operations": [r.to_dict() for r in results],
                "errors": [r.message for r in failed_ops],
            }
        )

    # Finalize molecule
    edited_mol = editor.finalize(
        optimize=optimize_geometry,
        add_hydrogens=True,
        generate_3d=return_3d
    )

    if edited_mol is None:
        return create_error_response(
            "FINALIZE_FAILED",
            "Failed to finalize edited molecule",
            {"state": editor.get_state().__dict__}
        )

    # Run validation if requested
    validation_result = None
    if should_validate:
        validation_result = run_validation(edited_mol)

    # Build response
    structure = molecule_to_structure(edited_mol, vacuum) if return_3d else {}

    annotations = {
        "stereocenters": [s.to_dict() for s in edited_mol.stereocenters],
        "rings": [r.to_dict() for r in edited_mol.rings],
        "functional_groups": [f.to_dict() for f in edited_mol.functional_groups],
        "atom_labels": edited_mol.atom_labels,
        "edit_history": edited_mol.edit_history,
    }

    edit_result = {
        "operations_applied": [r.to_dict() for r in results],
        "successful_count": len(successful_ops),
        "failed_count": len(failed_ops),
    }

    if validation_result:
        edit_result["validation"] = validation_result.to_dict()

    response = {
        "success": True,
        "structure": structure,
        "annotations": annotations,
        "edit_result": edit_result,
    }

    # Add warnings if some operations failed
    if failed_ops:
        response["warnings"] = [
            f"Operation failed: {r.message}" for r in failed_ops
        ]

    return response


def main():
    """Main entry point."""
    # Check for input file
    if len(sys.argv) < 2:
        print(json.dumps(create_error_response(
            "USAGE_ERROR",
            "Usage: python molecule_editor.py <input.json>"
        )))
        sys.exit(1)

    input_file = sys.argv[1]

    # Check file exists
    if not os.path.exists(input_file):
        print(json.dumps(create_error_response(
            "FILE_NOT_FOUND",
            f"Input file not found: {input_file}"
        )))
        sys.exit(1)

    # Read input file
    with open(input_file, 'r') as f:
        params = json.load(f)

    # Execute edit
    result = edit_molecule(params)

    # Output result
    print(json.dumps(result))


if __name__ == "__main__":
    main()
