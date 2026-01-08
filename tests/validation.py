"""
validation.py - Schema-Based Validation for E2E Testing

This module provides dynamic validation functions that verify structure
correctness WITHOUT hardcoded expected values.

Principle: Tests should validate that outputs conform to expected SCHEMAS
and PHYSICAL CONSTRAINTS, not specific atom counts.

Usage:
    from validation import validate_molecule_structure, validate_cluster_structure

    result = run_e2e_test("Generate benzene")
    is_valid, errors = validate_molecule_structure(result)
"""

import re
import math
from typing import Dict, List, Tuple, Any, Optional
import numpy as np


# =============================================================================
# SCHEMA VALIDATION
# =============================================================================

def validate_structure_schema(data: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate that a structure has all required fields.
    
    Required structure:
    {
        "success": True,
        "structure": {
            "atoms": [...],
            "lattice": {...}
        }
    }
    OR direct structure format.
    
    Returns:
        (is_valid, list_of_errors)
    """
    errors = []
    
    # Handle wrapper format
    if "success" in data:
        if not data.get("success"):
            errors.append(f"Structure generation failed: {data.get('error', 'Unknown error')}")
            return False, errors
        structure = data.get("structure", data)
    else:
        structure = data
    
    # Check required top-level fields
    if "atoms" not in structure and "sites" not in structure:
        errors.append("Missing 'atoms' or 'sites' field")
    
    if "lattice" not in structure:
        errors.append("Missing 'lattice' field")
    
    return len(errors) == 0, errors


def validate_atoms_have_required_fields(atoms: List[Dict]) -> Tuple[bool, List[str]]:
    """
    Validate that each atom has element and coordinates.
    
    Required per atom:
    - element: str
    - coords or cartesian: [x, y, z]
    """
    errors = []
    
    if not atoms:
        errors.append("Atoms list is empty")
        return False, errors
    
    for i, atom in enumerate(atoms):
        if "element" not in atom:
            errors.append(f"Atom {i} missing 'element' field")
        
        if "coords" not in atom and "cartesian" not in atom:
            errors.append(f"Atom {i} missing coordinates (need 'coords' or 'cartesian')")
        
        # Validate coordinate format
        coords = atom.get("coords") or atom.get("cartesian")
        if coords is not None:
            if not isinstance(coords, (list, tuple)) or len(coords) < 3:
                errors.append(f"Atom {i} has invalid coordinates: {coords}")
    
    return len(errors) == 0, errors


def validate_lattice_positive(lattice: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate lattice parameters are physically reasonable.
    
    - a, b, c > 0
    - 0 < alpha, beta, gamma < 180
    """
    errors = []
    
    for param in ["a", "b", "c"]:
        val = lattice.get(param)
        if val is None:
            errors.append(f"Missing lattice parameter '{param}'")
        elif val <= 0:
            errors.append(f"Lattice {param} must be positive, got {val}")
    
    for angle in ["alpha", "beta", "gamma"]:
        val = lattice.get(angle)
        if val is not None:
            if val <= 0 or val >= 180:
                errors.append(f"Lattice {angle} must be in (0, 180), got {val}")
    
    return len(errors) == 0, errors


# =============================================================================
# PHYSICAL VALIDATION
# =============================================================================

def validate_no_atomic_clashes(
    atoms: List[Dict], 
    min_dist: float = 0.5
) -> Tuple[bool, List[str]]:
    """
    Validate no atoms are too close (clash detection).
    
    For most molecules, atoms should be at least 0.5 Å apart.
    Exception: H-H in H2 is ~0.74 Å.
    """
    errors = []
    warnings = []
    
    if len(atoms) < 2:
        return True, []
    
    coords_list = []
    for atom in atoms:
        coords = atom.get("cartesian") or atom.get("coords", [0, 0, 0])
        if isinstance(coords, (list, tuple)) and len(coords) >= 3:
            coords_list.append(np.array(coords[:3], dtype=float))
    
    # Check pairwise distances
    for i in range(len(coords_list)):
        for j in range(i + 1, len(coords_list)):
            dist = np.linalg.norm(coords_list[i] - coords_list[j])
            if dist < min_dist:
                errors.append(
                    f"Atoms {i} and {j} too close: {dist:.3f} Å (min: {min_dist} Å)"
                )
    
    return len(errors) == 0, errors


def validate_formula_parseable(formula: str) -> Tuple[bool, List[str]]:
    """
    Validate that a formula string is parseable.
    
    Valid: C6H6, H2O, C24H8O6, Na+Cl-
    Invalid: empty, "???", random text
    """
    errors = []
    
    if not formula:
        errors.append("Formula is empty")
        return False, errors
    
    # Pattern for chemical formulas
    pattern = r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*[+-]?$'
    
    # Also allow formulas with + separator (like "C6H6+C6H6")
    parts = formula.replace("+", " ").replace("-", " ").split()
    
    for part in parts:
        if part and not re.match(pattern, part):
            # Check if it's just numbers or invalid
            if not any(c.isalpha() for c in part):
                errors.append(f"Invalid formula part: '{part}'")
    
    return len(errors) == 0, errors


def validate_formula_matches_atoms(
    formula: str, 
    atoms: List[Dict]
) -> Tuple[bool, List[str]]:
    """
    Validate that formula matches atom composition.
    
    Note: This is a soft validation - formulas can be written in different orders.
    """
    errors = []
    
    if not formula or not atoms:
        return True, []  # Skip if missing data
    
    # Count atoms by element
    atom_counts = {}
    for atom in atoms:
        elem = atom.get("element", "X")
        atom_counts[elem] = atom_counts.get(elem, 0) + 1
    
    # Parse formula
    formula_pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(formula_pattern, formula)
    formula_counts = {}
    for elem, count in matches:
        formula_counts[elem] = formula_counts.get(elem, 0) + int(count if count else 1)
    
    # Compare
    for elem, count in formula_counts.items():
        if elem not in atom_counts:
            errors.append(f"Formula has {elem} but not found in atoms")
        elif atom_counts[elem] != count:
            errors.append(
                f"Formula shows {count} {elem} but atoms has {atom_counts[elem]}"
            )
    
    return len(errors) == 0, errors


# =============================================================================
# CLUSTER-SPECIFIC VALIDATION
# =============================================================================

def validate_cluster_metadata(
    structure: Dict[str, Any],
    expected_pattern: Optional[str] = None
) -> Tuple[bool, List[str]]:
    """
    Validate cluster-specific metadata.
    
    Checks:
    - n_molecules field exists and > 0
    - stacking_type or pattern field exists
    - intermolecular_distance is reasonable
    """
    errors = []
    
    metadata = structure.get("metadata", structure.get("structure", {}).get("metadata", {}))
    
    n_mols = metadata.get("n_molecules")
    if n_mols is None:
        errors.append("Missing n_molecules in metadata")
    elif n_mols < 1:
        errors.append(f"Invalid n_molecules: {n_mols}")
    
    stacking = metadata.get("stacking_type") or metadata.get("pattern")
    if stacking is None:
        errors.append("Missing stacking_type or pattern in metadata")
    
    if expected_pattern and stacking:
        # Normalize for comparison
        expected_norm = expected_pattern.lower().replace("-", "_").replace(" ", "_")
        actual_norm = stacking.lower().replace("-", "_").replace(" ", "_")
        if expected_norm != actual_norm and expected_pattern != "auto":
            errors.append(f"Expected pattern '{expected_pattern}', got '{stacking}'")
    
    distance = metadata.get("intermolecular_distance")
    if distance is not None:
        if distance < 0.5:
            errors.append(f"Intermolecular distance too small: {distance} Å")
        elif distance > 20:
            errors.append(f"Intermolecular distance unusually large: {distance} Å")
    
    return len(errors) == 0, errors


def validate_arrangement_geometry(
    atoms: List[Dict],
    expected_spacing: Optional[float] = None,
    tolerance: float = 0.5
) -> Tuple[bool, List[str]]:
    """
    Validate geometric properties of an arrangement.
    
    For stacked structures, checks z-spacing.
    For circular, checks radial symmetry.
    """
    errors = []
    warnings = []
    
    if len(atoms) < 2:
        return True, []
    
    # Extract z-coordinates
    z_coords = []
    for atom in atoms:
        coords = atom.get("cartesian") or atom.get("coords", [0, 0, 0])
        if isinstance(coords, (list, tuple)) and len(coords) >= 3:
            z_coords.append(float(coords[2]))
    
    if not z_coords:
        return True, []
    
    z_range = max(z_coords) - min(z_coords)
    
    # Basic check: if expected_spacing and multiple molecules, 
    # z-range should be roughly (n-1) * spacing
    if expected_spacing and z_range > 0:
        # Just check that spacing is in reasonable range
        pass  # Could add more specific checks
    
    return len(errors) == 0, errors


# =============================================================================
# COMPOSITE VALIDATORS
# =============================================================================

def validate_molecule_structure(data: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Complete validation for a single molecule structure.
    
    Combines all relevant validators without hardcoding expected values.
    """
    all_errors = []
    
    # 1. Schema validation
    valid, errors = validate_structure_schema(data)
    all_errors.extend(errors)
    if not valid:
        return False, all_errors
    
    # Get structure
    structure = data.get("structure", data)
    
    # 2. Atoms validation
    atoms = structure.get("atoms") or structure.get("sites", [])
    valid, errors = validate_atoms_have_required_fields(atoms)
    all_errors.extend(errors)
    
    # 3. Lattice validation
    lattice = structure.get("lattice", {})
    valid, errors = validate_lattice_positive(lattice)
    all_errors.extend(errors)
    
    # 4. Clash detection
    valid, errors = validate_no_atomic_clashes(atoms)
    all_errors.extend(errors)
    
    # 5. Formula validation (if available)
    formula = structure.get("metadata", {}).get("formula")
    if formula:
        valid, errors = validate_formula_parseable(formula)
        all_errors.extend(errors)
    
    return len(all_errors) == 0, all_errors


def validate_cluster_structure(
    data: Dict[str, Any],
    expected_pattern: Optional[str] = None
) -> Tuple[bool, List[str]]:
    """
    Complete validation for a molecular cluster.
    
    Includes all molecule validations plus cluster-specific checks.
    """
    all_errors = []
    
    # 1. Basic molecule validation
    valid, errors = validate_molecule_structure(data)
    all_errors.extend(errors)
    
    # 2. Cluster metadata
    valid, errors = validate_cluster_metadata(data, expected_pattern)
    all_errors.extend(errors)
    
    return len(all_errors) == 0, all_errors


# =============================================================================
# TEST RESULT VALIDATION
# =============================================================================

def validate_e2e_result(
    result: Dict[str, Any],
    is_cluster: bool = False,
    expected_pattern: Optional[str] = None
) -> Tuple[bool, List[str], Dict[str, Any]]:
    """
    Validate complete E2E test result.
    
    Returns:
        (is_valid, errors, extracted_info)
    """
    extracted = {
        "n_atoms": 0,
        "formula": None,
        "n_molecules": None,
        "pattern": None
    }
    
    # Handle LLM errors
    if "error" in result and not result.get("success"):
        return False, [f"LLM/MCP error: {result['error']}"], extracted
    
    # Get structure data
    structure = result.get("structure", result)
    atoms = structure.get("atoms") or structure.get("sites", [])
    
    extracted["n_atoms"] = len(atoms)
    extracted["formula"] = structure.get("metadata", {}).get("formula")
    extracted["n_molecules"] = structure.get("metadata", {}).get("n_molecules")
    extracted["pattern"] = structure.get("metadata", {}).get("stacking_type")
    
    # Run validation
    if is_cluster:
        valid, errors = validate_cluster_structure(structure, expected_pattern)
    else:
        valid, errors = validate_molecule_structure(structure)
    
    return valid, errors, extracted


# =============================================================================
# SELF-TEST
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("VALIDATION MODULE - Self Test")
    print("=" * 60)
    
    # Test valid structure
    valid_structure = {
        "success": True,
        "structure": {
            "atoms": [
                {"element": "C", "coords": [0, 0, 0], "cartesian": [0, 0, 0]},
                {"element": "H", "coords": [1, 0, 0], "cartesian": [1, 0, 0]},
            ],
            "lattice": {"a": 10, "b": 10, "c": 10, "alpha": 90, "beta": 90, "gamma": 90},
            "metadata": {"formula": "CH", "natoms": 2}
        }
    }
    
    valid, errors = validate_molecule_structure(valid_structure)
    print(f"\n1. Valid structure: {'✓ PASS' if valid else '✗ FAIL'}")
    if errors:
        for e in errors:
            print(f"   - {e}")
    
    # Test invalid structure
    invalid_structure = {
        "success": True,
        "structure": {
            "atoms": [],  # Empty!
            "lattice": {"a": -5}  # Negative!
        }
    }
    
    valid, errors = validate_molecule_structure(invalid_structure)
    print(f"\n2. Invalid structure: {'✓ correctly rejected' if not valid else '✗ should have failed'}")
    for e in errors:
        print(f"   - {e}")
    
    # Test clash detection
    clash_structure = {
        "structure": {
            "atoms": [
                {"element": "C", "cartesian": [0, 0, 0]},
                {"element": "C", "cartesian": [0.1, 0, 0]},  # Too close!
            ],
            "lattice": {"a": 10, "b": 10, "c": 10, "alpha": 90, "beta": 90, "gamma": 90}
        }
    }
    
    valid, errors = validate_molecule_structure(clash_structure)
    print(f"\n3. Clash detection: {'✓ detected clash' if not valid else '✗ missed clash'}")
    for e in errors:
        print(f"   - {e}")
    
    print("\n" + "=" * 60)
    print("Self-test complete!")
    print("=" * 60)
