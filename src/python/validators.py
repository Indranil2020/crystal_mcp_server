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
from pymatgen.core import Structure


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
        warnings.append("Very small lattice parameters (< 2.0 Å)")
    if a > 50.0 or b > 50.0 or c > 50.0:
        warnings.append("Very large lattice parameters (> 50.0 Å)")

    if alpha < 30 or alpha > 150 or beta < 30 or beta > 150 or gamma < 30 or gamma > 150:
        warnings.append("Unusual lattice angles (< 30° or > 150°)")

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
            try:
                matrix = np.array(scaling)
                det = np.linalg.det(matrix)
                if abs(det) < 1e-6:
                    return False, "Scaling matrix is singular (determinant ~ 0)"
                return True, None
            except Exception:
                return False, "Invalid scaling matrix"

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
        checks = ["distances", "lattice", "stoichiometry"]
    
    # Validate structure format
    if not structure_dict or not isinstance(structure_dict, dict):
        result.add_error("Structure must be a dictionary")
        return {"success": True, **result.to_dict()}
    
    if "lattice" not in structure_dict:
        result.add_error("Structure missing lattice information")
        return {"success": True, **result.to_dict()}
    
    if "atoms" not in structure_dict:
        result.add_error("Structure missing atoms information")
        return {"success": True, **result.to_dict()}
    
    lattice = structure_dict["lattice"]
    atoms = structure_dict["atoms"]
    
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
    
    # Check lattice parameters
    if "lattice" in checks:
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
    
    # Calculate density
    density = calculate_density(atoms, lattice)
    if density is not None:
        result.add_metric("density", density)
        
        if density < 0.1:
            result.add_warning(f"Very low density: {density:.3f} g/cm³")
        elif density > 30:
            result.add_warning(f"Very high density: {density:.3f} g/cm³")
    
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
            f"Atoms {idx1} ({elem1}) and {idx2} ({elem2}) too close: {min_dist:.3f} Å < {min_distance:.3f} Å",
            "error",
            {"min_distance": min_dist, "threshold": min_distance}
        )
    elif min_dist < 0.5:
        elem1, elem2, idx1, idx2 = min_pair
        result.add_error(
            f"Atoms {idx1} ({elem1}) and {idx2} ({elem2}) overlapping: {min_dist:.3f} Å",
            "error",
            {"min_distance": min_dist}
        )
    elif min_dist < 1.0:
        result.add_warning(f"Very short interatomic distance: {min_dist:.3f} Å")
    
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
        result.add_warning(f"Very small lattice parameters: a={a:.3f}, b={b:.3f}, c={c:.3f} Å")
    
    if a > 100 or b > 100 or c > 100:
        result.add_warning(f"Very large lattice parameters: a={a:.3f}, b={b:.3f}, c={c:.3f} Å")
    
    # Check angles
    if not (0 < alpha < 180 and 0 < beta < 180 and 0 < gamma < 180):
        result.add_error(
            f"Lattice angles must be in (0, 180): α={alpha}, β={beta}, γ={gamma}",
            "error"
        )
    
    if alpha < 20 or beta < 20 or gamma < 20:
        result.add_warning(f"Very acute angles: α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")
    
    if alpha > 160 or beta > 160 or gamma > 160:
        result.add_warning(f"Very obtuse angles: α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")
    
    # Check volume
    if "volume" in lattice:
        volume = lattice["volume"]
        if volume <= 0:
            result.add_error(f"Volume must be positive: {volume}", "error")
        elif volume < 10:
            result.add_warning(f"Very small unit cell volume: {volume:.3f} ų")
        elif volume > 100000:
            result.add_warning(f"Very large unit cell volume: {volume:.3f} ų")
    
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
    Calculate crystal density in g/cm³.
    
    Args:
        atoms: List of atoms
        lattice: Lattice parameters
    
    Returns:
        Density in g/cm³ or None
    """
    if not atoms or "volume" not in lattice:
        return None
    
    # Atomic masses (simplified)
    atomic_masses = {
        'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
        'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
        'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
        'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845,
        'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.380
    }
    
    total_mass = 0.0
    for atom in atoms:
        if "element" not in atom:
            continue
        element = atom["element"]
        mass = atomic_masses.get(element, 1.0)
        total_mass += mass
    
    volume_cm3 = lattice["volume"] * 1e-24  # ų to cm³
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
    
    # Covalent radii (simplified, in Angstroms)
    covalent_radii = {
        'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
        'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
        'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
        'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
        'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32,
        'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22
    }
    
    total_volume = 0.0
    for atom in atoms:
        if "element" not in atom:
            continue
        element = atom["element"]
        radius = covalent_radii.get(element, 1.0)
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
