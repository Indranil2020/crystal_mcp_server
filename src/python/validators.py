#!/usr/bin/env python3
"""
Validators - Structure Quality Checks

This module provides comprehensive validation utilities for crystal structures.
All functions follow defensive programming without try/except blocks.
"""

from typing import Dict, List, Optional, Any, Tuple
import json
import sys
import numpy as np
import spglib
from pymatgen.core import Structure, Element


class ValidationResult:
    """Result of validation check."""
    
    def __init__(self):
        self.valid = True
        self.errors = []
        self.warnings = []
        self.metrics = {}
    
    def add_error(self, message: str, severity: str = "error", details: Any = None):
        """Add error message."""
        self.errors.append({
            "type": "validation",
            "severity": severity,
            "message": message,
            "details": details or {}
        })
        if severity == "error":
            self.valid = False
    
    def add_warning(self, message: str, details: Any = None):
        """Add warning message."""
        self.warnings.append({
            "message": message,
            "details": details or {}
        })
    
    def add_metric(self, name: str, value: Any):
        """Add validation metric."""
        self.metrics[name] = value
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "valid": self.valid,
            "errors": self.errors,
            "warnings": self.warnings,
            "metrics": self.metrics
        }


def validate_lattice_parameters(
    a: float, b: float, c: float,
    alpha: float, beta: float, gamma: float
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate lattice parameters.

    Args:
        a, b, c: Lattice lengths in Angstroms
        alpha, beta, gamma: Lattice angles in degrees

    Returns:
        Tuple of (is_valid, errors, warnings)
    """
    errors = []
    warnings = []

    # Check lengths are positive
    if a <= 0:
        errors.append(f"Lattice parameter a must be positive, got {a}")
    if b <= 0:
        errors.append(f"Lattice parameter b must be positive, got {b}")
    if c <= 0:
        errors.append(f"Lattice parameter c must be positive, got {c}")

    # Check angles are in valid range
    if alpha <= 0 or alpha >= 180:
        errors.append(f"Angle alpha must be in range (0, 180), got {alpha}")
    if beta <= 0 or beta >= 180:
        errors.append(f"Angle beta must be in range (0, 180), got {beta}")
    if gamma <= 0 or gamma >= 180:
        errors.append(f"Angle gamma must be in range (0, 180), got {gamma}")

    # Check for unusual values
    if a < 2.0 or b < 2.0 or c < 2.0:
        warnings.append("Very small lattice parameters (< 2.0 Angstroms)")
    if a > 50.0 or b > 50.0 or c > 50.0:
        warnings.append("Very large lattice parameters (> 50.0 Angstroms)")

    if alpha < 30 or alpha > 150 or beta < 30 or beta > 150 or gamma < 30 or gamma > 150:
        warnings.append("Unusual lattice angles (< 30 deg or > 150 deg)")

    is_valid = len(errors) == 0
    return is_valid, errors, warnings


def validate_miller_indices(miller: List[int]) -> Tuple[bool, Optional[str]]:
    """
    Validate Miller indices for a crystallographic plane.

    Args:
        miller: List of 3 integers [h, k, l]

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(miller, (list, tuple)):
        return False, "Miller indices must be a list or tuple"

    if len(miller) != 3:
        return False, f"Miller indices must have 3 values, got {len(miller)}"

    if not all(isinstance(x, int) for x in miller):
        return False, "Miller indices must be integers"

    if all(x == 0 for x in miller):
        return False, "Miller indices cannot all be zero"

    return True, None


def validate_scaling_matrix(scaling: Any) -> Tuple[bool, Optional[str]]:
    """
    Validate a scaling matrix or vector for supercell creation.

    Args:
        scaling: Either [nx, ny, nz] vector or 3x3 matrix

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(scaling, (list, tuple, np.ndarray)):
        return False, "Scaling must be a list, tuple, or array"

    if isinstance(scaling, np.ndarray):
        scaling = scaling.tolist()

    # Check if it's a vector [nx, ny, nz]
    if len(scaling) == 3 and all(isinstance(x, (int, float)) for x in scaling):
        if any(x <= 0 for x in scaling):
            return False, "Scaling factors must be positive"
        if any(x != int(x) for x in scaling):
            return False, "Scaling factors must be integers"
        return True, None

    # Check if it's a 3x3 matrix
    if len(scaling) == 3:
        if all(isinstance(row, (list, tuple)) and len(row) == 3 for row in scaling):
            # Validate it's a proper transformation matrix
            matrix = np.array(scaling)
            det = np.linalg.det(matrix)
            if abs(det) < 1e-6:
                return False, "Scaling matrix is singular (determinant ~ 0)"
            return True, None

    return False, "Scaling must be [nx, ny, nz] or 3x3 matrix"


def check_atoms_too_close(
    structure_dict: Dict[str, Any],
    threshold: float = 1.0
) -> Dict[str, Any]:
    """
    Check if any atoms are too close to each other.

    Args:
        structure_dict: Structure to check
        threshold: Minimum distance threshold in Angstroms

    Returns:
        Dictionary with check results
    """
    if "atoms" not in structure_dict:
        return {
            "atoms_too_close": False,
            "min_distance": None,
            "close_pairs": [],
            "error": "No atoms in structure"
        }

    atoms = structure_dict["atoms"]
    if len(atoms) < 2:
        return {
            "atoms_too_close": False,
            "min_distance": None,
            "close_pairs": []
        }

    close_pairs = []
    min_distance = float('inf')

    # Check all pairs
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            atom1 = atoms[i]
            atom2 = atoms[j]

            if "cartesian" in atom1 and "cartesian" in atom2:
                pos1 = np.array(atom1["cartesian"])
                pos2 = np.array(atom2["cartesian"])
                distance = np.linalg.norm(pos1 - pos2)

                min_distance = min(min_distance, distance)

                if distance < threshold:
                    close_pairs.append({
                        "atom1_index": i,
                        "atom2_index": j,
                        "atom1_element": atom1.get("element", "?"),
                        "atom2_element": atom2.get("element", "?"),
                        "distance": float(distance)
                    })

    return {
        "atoms_too_close": len(close_pairs) > 0,
        "min_distance": float(min_distance) if min_distance != float('inf') else None,
        "close_pairs": close_pairs,
        "threshold": threshold
    }


def validate_structure(
    structure_dict: Dict[str, Any],
    checks: Optional[List[str]] = None,
    min_distance: Optional[float] = None,
    max_coordination: Optional[int] = None,
    expected_space_group: Optional[int] = None
) -> Dict[str, Any]:
    """
    Comprehensive structure validation.
    
    Args:
        structure_dict: Structure to validate
        checks: List of checks to perform
        min_distance: Minimum allowed interatomic distance
        max_coordination: Maximum coordination number
        expected_space_group: Expected space group number
    
    Returns:
        Dictionary with validation results
    """
    result = ValidationResult()
    
    # Default checks
    if checks is None:
        checks = ["distances", "symmetry", "overlaps"]
    
    # Validate structure format
    if not structure_dict or not isinstance(structure_dict, dict):
        result.add_error("Structure must be a dictionary")
        return {"success": True, **result.to_dict()}
    
    if "lattice" not in structure_dict:
        result.add_error("Structure missing lattice information")
        return {"success": True, **result.to_dict()}
    
    if "atoms" not in structure_dict and "sites" not in structure_dict:
        result.add_error("Structure missing atoms information")
        return {"success": True, **result.to_dict()}
    
    lattice = structure_dict["lattice"]
    atoms = structure_dict.get("atoms") or structure_dict.get("sites")
    
    # Check distances
    if "distances" in checks:
        distance_result = check_interatomic_distances(atoms, lattice, min_distance)
        if not distance_result["valid"]:
            for error in distance_result["errors"]:
                result.add_error(error["message"], error["severity"], error.get("details"))
        for warning in distance_result.get("warnings", []):
            result.add_warning(warning)
        if "min_distance" in distance_result.get("metrics", {}):
            result.add_metric("min_distance", distance_result["metrics"]["min_distance"])
        if "avg_distance" in distance_result.get("metrics", {}):
            result.add_metric("avg_distance", distance_result["metrics"]["avg_distance"])

    # Check overlaps (reuse distance logic with overlap threshold)
    if "overlaps" in checks and "distances" not in checks:
        overlap_result = check_interatomic_distances(atoms, lattice, min_distance or 0.5)
        if not overlap_result["valid"]:
            for error in overlap_result["errors"]:
                result.add_error(error["message"], error["severity"], error.get("details"))
        for warning in overlap_result.get("warnings", []):
            result.add_warning(warning)
    
    # Check lattice parameters
    if "lattice" in checks or "cell_parameters" in checks:
        lattice_result = check_lattice_parameters(lattice)
        if not lattice_result["valid"]:
            for error in lattice_result["errors"]:
                result.add_error(error["message"], error["severity"], error.get("details"))
        for warning in lattice_result.get("warnings", []):
            result.add_warning(warning)
    
    # Check stoichiometry
    if "stoichiometry" in checks:
        stoich_result = check_stoichiometry(atoms)
        for warning in stoich_result.get("warnings", []):
            result.add_warning(warning)

    # Check symmetry
    if "symmetry" in checks:
        symmetry_result = check_symmetry(atoms, lattice, expected_space_group)
        if not symmetry_result["valid"]:
            for error in symmetry_result["errors"]:
                result.add_error(error["message"], error["severity"], error.get("details"))
        for warning in symmetry_result.get("warnings", []):
            result.add_warning(warning)
        for metric, value in symmetry_result.get("metrics", {}).items():
            result.add_metric(metric, value)

    # Check charge balance
    if "charge_balance" in checks:
        charge_result = check_charge_balance(atoms)
        if not charge_result["valid"]:
            for error in charge_result["errors"]:
                result.add_error(error["message"], error["severity"], error.get("details"))
        for warning in charge_result.get("warnings", []):
            result.add_warning(warning)
        for metric, value in charge_result.get("metrics", {}).items():
            result.add_metric(metric, value)
    
    # Calculate density
    density = calculate_density(atoms, lattice)
    if density is not None:
        result.add_metric("density", density)
        
        if density < 0.1:
            result.add_warning(f"Very low density: {density:.3f} g/cm^3")
        elif density > 30:
            result.add_warning(f"Very high density: {density:.3f} g/cm^3")
    
    # Calculate packing fraction
    packing = calculate_packing_fraction(atoms, lattice)
    if packing is not None:
        result.add_metric("packing_fraction", packing)
        
        if packing > 0.8:
            result.add_warning(f"High packing fraction: {packing:.3f}")
    
    return {"success": True, **result.to_dict()}


def check_interatomic_distances(
    atoms: List[Dict],
    lattice: Dict,
    min_distance: Optional[float] = None
) -> Dict[str, Any]:
    """
    Check interatomic distances.
    
    Args:
        atoms: List of atoms
        lattice: Lattice parameters
        min_distance: Minimum allowed distance
    
    Returns:
        Validation result dictionary
    """
    result = ValidationResult()
    
    if not atoms or len(atoms) < 2:
        return result.to_dict()
    
    # Get lattice matrix
    if "matrix" not in lattice:
        result.add_error("Lattice missing matrix")
        return result.to_dict()
    
    lattice_matrix = np.array(lattice["matrix"])
    
    # Calculate all pairwise distances
    distances = []
    min_dist = float('inf')
    min_pair = None
    
    for i in range(len(atoms)):
        if "coords" not in atoms[i] or "element" not in atoms[i]:
            continue
        
        frac_i = np.array(atoms[i]["coords"])
        
        for j in range(i + 1, len(atoms)):
            if "coords" not in atoms[j] or "element" not in atoms[j]:
                continue
            
            frac_j = np.array(atoms[j]["coords"])
            
            # Calculate minimum image distance
            delta = frac_j - frac_i
            delta = delta - np.round(delta)  # Minimum image convention
            cart_delta = np.dot(delta, lattice_matrix)
            dist = float(np.linalg.norm(cart_delta))
            
            distances.append(dist)
            
            if dist < min_dist:
                min_dist = dist
                min_pair = (atoms[i]["element"], atoms[j]["element"], i, j)
    
    if len(distances) == 0:
        return result.to_dict()
    
    # Check minimum distance
    result.add_metric("min_distance", min_dist)
    result.add_metric("avg_distance", float(np.mean(distances)))
    
    if min_distance is not None and min_dist < min_distance:
        elem1, elem2, idx1, idx2 = min_pair
        result.add_error(
            f"Atoms {idx1} ({elem1}) and {idx2} ({elem2}) too close: {min_dist:.3f} Angstroms < {min_distance:.3f} Angstroms",
            "error",
            {"min_distance": min_dist, "threshold": min_distance}
        )
    elif min_dist < 0.5:
        elem1, elem2, idx1, idx2 = min_pair
        result.add_error(
            f"Atoms {idx1} ({elem1}) and {idx2} ({elem2}) overlapping: {min_dist:.3f} Angstroms",
            "error",
            {"min_distance": min_dist}
        )
    elif min_dist < 1.0:
        result.add_warning(f"Very short interatomic distance: {min_dist:.3f} Angstroms")
    
    return result.to_dict()


def check_lattice_parameters(lattice: Dict) -> Dict[str, Any]:
    """
    Check lattice parameters are reasonable.
    
    Args:
        lattice: Lattice dictionary
    
    Returns:
        Validation result dictionary
    """
    result = ValidationResult()
    
    # Check required fields
    required = ["a", "b", "c", "alpha", "beta", "gamma"]
    for field in required:
        if field not in lattice:
            result.add_error(f"Lattice missing {field}")
            return result.to_dict()
    
    a, b, c = lattice["a"], lattice["b"], lattice["c"]
    alpha, beta, gamma = lattice["alpha"], lattice["beta"], lattice["gamma"]
    
    # Check lengths
    if a <= 0 or b <= 0 or c <= 0:
        result.add_error(
            f"Lattice parameters must be positive: a={a}, b={b}, c={c}",
            "error"
        )
    
    if a < 1.0 or b < 1.0 or c < 1.0:
        result.add_warning(f"Very small lattice parameters: a={a:.3f}, b={b:.3f}, c={c:.3f} Angstroms")
    
    if a > 100 or b > 100 or c > 100:
        result.add_warning(f"Very large lattice parameters: a={a:.3f}, b={b:.3f}, c={c:.3f} Angstroms")
    
    # Check angles
    if not (0 < alpha < 180 and 0 < beta < 180 and 0 < gamma < 180):
        result.add_error(
            f"Lattice angles must be in (0, 180): alpha={alpha}, beta={beta}, gamma={gamma}",
            "error"
        )
    
    if alpha < 20 or beta < 20 or gamma < 20:
        result.add_warning(f"Very acute angles: alpha={alpha:.1f} deg, beta={beta:.1f} deg, gamma={gamma:.1f} deg")
    
    if alpha > 160 or beta > 160 or gamma > 160:
        result.add_warning(f"Very obtuse angles: alpha={alpha:.1f} deg, beta={beta:.1f} deg, gamma={gamma:.1f} deg")
    
    # Check volume
    if "volume" in lattice:
        volume = lattice["volume"]
        if volume <= 0:
            result.add_error(f"Volume must be positive: {volume}", "error")
        elif volume < 10:
            result.add_warning(f"Very small unit cell volume: {volume:.3f} A^3")
        elif volume > 100000:
            result.add_warning(f"Very large unit cell volume: {volume:.3f} A^3")
    
    return result.to_dict()


def check_symmetry(
    atoms: List[Dict],
    lattice: Dict,
    expected_space_group: Optional[int] = None
) -> Dict[str, Any]:
    """
    Check symmetry using spglib and compare with expected space group if provided.
    """
    result = ValidationResult()

    if "matrix" not in lattice:
        result.add_error("Lattice missing matrix")
        return result.to_dict()

    if not atoms:
        result.add_error("No atoms provided for symmetry check")
        return result.to_dict()

    coords = []
    numbers = []
    for atom in atoms:
        if "coords" not in atom or "element" not in atom:
            continue
        coords.append(atom["coords"])
        element = atom["element"]
        if Element.is_valid_symbol(element):
            numbers.append(Element(element).Z)
        else:
            numbers.append(0)

    if not coords:
        result.add_error("No valid coordinates for symmetry check")
        return result.to_dict()

    cell = (np.array(lattice["matrix"]), np.array(coords), np.array(numbers))
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-3)
    if dataset is None:
        result.add_warning("Symmetry detection failed (spglib returned no dataset)")
        return result.to_dict()

    sg_number = dataset.get("number")
    sg_symbol = dataset.get("international", "")
    result.add_metric("space_group_number", int(sg_number) if sg_number is not None else None)
    result.add_metric("space_group_symbol", sg_symbol)

    if expected_space_group is not None and sg_number is not None:
        if int(sg_number) != int(expected_space_group):
            result.add_error(
                f"Expected space group {expected_space_group}, got {sg_number}",
                "error",
                {"expected": expected_space_group, "actual": int(sg_number)}
            )

    return result.to_dict()


def check_charge_balance(atoms: List[Dict]) -> Dict[str, Any]:
    """
    Check charge balance if per-atom charges are provided.
    """
    result = ValidationResult()
    charges = []

    for atom in atoms:
        if "charge" in atom:
            charges.append(float(atom["charge"]))
        elif "oxidation_state" in atom:
            charges.append(float(atom["oxidation_state"]))

    if not charges:
        result.add_warning("Charge balance check skipped (no per-atom charges provided)")
        return result.to_dict()

    total_charge = float(np.sum(charges))
    result.add_metric("total_charge", total_charge)

    if abs(total_charge) > 1e-3:
        result.add_error(
            f"Structure is not charge balanced (total charge {total_charge:.3f})",
            "error",
            {"total_charge": total_charge}
        )

    return result.to_dict()


def check_stoichiometry(atoms: List[Dict]) -> Dict[str, Any]:
    """
    Check stoichiometry for unusual compositions.
    
    Args:
        atoms: List of atoms
    
    Returns:
        Validation result dictionary
    """
    result = ValidationResult()
    
    if not atoms:
        result.add_warning("No atoms in structure")
        return result.to_dict()
    
    # Count elements
    element_counts = {}
    for atom in atoms:
        if "element" not in atom:
            continue
        element = atom["element"]
        element_counts[element] = element_counts.get(element, 0) + 1
    
    n_elements = len(element_counts)
    total_atoms = sum(element_counts.values())
    
    if n_elements == 0:
        result.add_warning("No valid elements found")
    elif n_elements > 10:
        result.add_warning(f"Very complex composition: {n_elements} different elements")
    
    if total_atoms > 1000:
        result.add_warning(f"Very large structure: {total_atoms} atoms")
    
    return result.to_dict()


def calculate_density(atoms: List[Dict], lattice: Dict) -> Optional[float]:
    """
    Calculate crystal density in g/cm^3.
    
    Args:
        atoms: List of atoms
        lattice: Lattice parameters
    
    Returns:
        Density in g/cm^3 or None
    """
    if not atoms or "volume" not in lattice:
        return None

    # Complete atomic masses (in atomic mass units)
    atomic_masses = {
        'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
        'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
        'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
        'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
        'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.380,
        'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
        'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
        'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.91,
        'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71,
        'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
        'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
        'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93,
        'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
        'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
        'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
        'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209.0, 'At': 210.0,
        'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.04,
        'Pa': 231.04, 'U': 238.03, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0
    }

    total_mass = 0.0
    for atom in atoms:
        if "element" not in atom:
            continue
        element = atom["element"]
        mass = atomic_masses.get(element, 1.0)  # Default to 1.0 for unknown
        total_mass += mass
    
    volume_cm3 = lattice["volume"] * 1e-24  # A^3 to cm^3
    avogadro = 6.022e23
    
    density = (total_mass / avogadro) / volume_cm3
    return float(density)


def calculate_packing_fraction(atoms: List[Dict], lattice: Dict) -> Optional[float]:
    """
    Estimate packing fraction.
    
    Args:
        atoms: List of atoms
        lattice: Lattice parameters
    
    Returns:
        Packing fraction or None
    """
    if not atoms or "volume" not in lattice:
        return None
    
    def atomic_radius(element: str) -> Optional[float]:
        if not Element.is_valid_symbol(element):
            return None
        el = Element(element)
        radius = el.atomic_radius or el.van_der_waals_radius or el.covalent_radius
        return float(radius) if radius else None

    total_volume = 0.0
    for atom in atoms:
        if "element" not in atom:
            continue
        element = atom["element"]
        radius = atomic_radius(element) or 1.0  # Default to 1.0 Angstrom for unknown
        atom_volume = (4.0 / 3.0) * np.pi * radius ** 3
        total_volume += atom_volume
    
    packing_fraction = total_volume / lattice["volume"]
    return float(packing_fraction)


def main():
    """
    Main entry point for command-line execution.
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python validators.py <input_json_file>"
            }
        }))
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input file
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    # Validate structure
    result = validate_structure(**params)
    
    # Output result
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
