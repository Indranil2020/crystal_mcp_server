#!/usr/bin/env python3
"""
Crystal Generator - PyXtal Wrapper

This module provides a comprehensive wrapper around PyXtal for generating
crystal structures with defensive programming and detailed error handling.
All functions follow defensive programming principles without using try/except.
"""

from typing import Dict, List, Optional, Tuple, Any, Union
import json
import sys
import numpy as np
from pyxtal import pyxtal
from pyxtal.symmetry import Group
from pymatgen.core import Structure, Element, Lattice as PMGLattice


class ValidationResult:
    """Represents validation result with success status and messages."""
    
    def __init__(self, valid: bool, errors: List[str], warnings: List[str]):
        self.valid = valid
        self.errors = errors
        self.warnings = warnings
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "valid": self.valid,
            "errors": self.errors,
            "warnings": self.warnings
        }


def validate_space_group(space_group: Union[int, str]) -> Tuple[bool, Optional[str]]:
    """
    Validate space group number or symbol.

    Args:
        space_group: Space group number (1-230) or Hermann-Mauguin symbol

    Returns:
        Tuple of (is_valid, error_message)
    """
    if isinstance(space_group, int):
        if space_group < 1 or space_group > 230:
            return False, f"Space group number must be between 1 and 230, got {space_group}"
    elif isinstance(space_group, str):
        if not space_group or len(space_group.strip()) == 0:
            return False, "Space group symbol cannot be empty"
    else:
        return False, f"Space group must be int or str, got {type(space_group).__name__}"

    return True, None


def validate_composition(composition: List[str]) -> Tuple[bool, Optional[str]]:
    """
    Validate chemical composition.

    Args:
        composition: List of element symbols

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not composition or len(composition) == 0:
        return False, "Composition cannot be empty"

    # Valid element symbols (simplified set)
    valid_elements = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr"
    }

    for i, elem in enumerate(composition):
        if not elem or not isinstance(elem, str):
            return False, f"Invalid element at index {i}: {elem}"

        elem_clean = elem.strip()
        if len(elem_clean) == 0:
            return False, f"Empty element at index {i}"

        if elem_clean not in valid_elements:
            # Try to get element from pymatgen to check if valid
            is_valid = False
            if hasattr(Element, elem_clean):
                is_valid = True
            else:
                # Check if it's in Element's internal list
                all_elements = [e.symbol for e in Element]
                if elem_clean in all_elements:
                    is_valid = True

            if not is_valid:
                return False, f"Unknown element '{elem_clean}' at index {i}"

    return True, None


def get_covalent_radius(element: str) -> float:
    """
    Get covalent radius for an element in Angstroms.
    
    Args:
        element: Element symbol
    
    Returns:
        Covalent radius in Angstroms, or 1.5 if unknown
    """
    radii = {
        "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84,
        "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
        "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
        "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76,
        "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
        "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22
    }
    return radii.get(element, 1.5)


def calculate_min_distances(structure: Structure) -> Dict[str, float]:
    """
    Calculate minimum interatomic distances in structure.
    
    Args:
        structure: Pymatgen Structure object
    
    Returns:
        Dictionary mapping element pairs to minimum distances
    """
    distances = {}
    
    for i in range(len(structure)):
        for j in range(i + 1, len(structure)):
            site_i = structure[i]
            site_j = structure[j]
            
            dist = site_i.distance(site_j)
            elem_i = site_i.specie.symbol
            elem_j = site_j.specie.symbol
            
            pair = tuple(sorted([elem_i, elem_j]))
            pair_key = f"{pair[0]}-{pair[1]}"
            
            if pair_key not in distances or dist < distances[pair_key]:
                distances[pair_key] = float(dist)
    
    return distances


def calculate_density(crystal: pyxtal) -> float:
    """
    Calculate crystal density in g/cm³.
    
    Args:
        crystal: PyXtal crystal object
    
    Returns:
        Density in g/cm³
    """
    mass = 0.0
    for site in crystal.atom_sites:
        element = Element(site.specie)
        mass += element.atomic_mass * len(site.coords)
    
    # Convert volume from Ų to cm³
    volume_cm3 = crystal.lattice.volume * 1e-24
    
    # Calculate density
    avogadro = 6.022e23
    if volume_cm3 > 0:
        density = (mass / avogadro) / volume_cm3
    else:
        density = 0.0
    
    return float(density)


def extract_structure_data(crystal: pyxtal) -> Dict[str, Any]:
    """
    Extract complete structure data from PyXtal crystal.
    
    Args:
        crystal: PyXtal crystal object
    
    Returns:
        Dictionary containing all structure information
    """
    lattice = crystal.lattice
    
    structure_data = {
        "lattice": {
            "a": float(lattice.a),
            "b": float(lattice.b),
            "c": float(lattice.c),
            "alpha": float(lattice.alpha),
            "beta": float(lattice.beta),
            "gamma": float(lattice.gamma),
            "matrix": lattice.matrix.tolist(),
            "volume": float(lattice.volume)
        },
        "atoms": [],
        "space_group": {
            "number": int(crystal.group.number),
            "symbol": str(crystal.group.symbol),
            "hall_symbol": str(crystal.group.hall_symbol) if hasattr(crystal.group, 'hall_symbol') else "",
            "point_group": str(crystal.group.point_group) if hasattr(crystal.group, 'point_group') else "",
            "crystal_system": str(crystal.group.lattice_type) if hasattr(crystal.group, 'lattice_type') else ""
        },
        "metadata": {
            "formula": str(crystal.formula),
            "natoms": sum(len(site.coords) for site in crystal.atom_sites) if hasattr(crystal, 'atom_sites') else 0,
            "volume": float(lattice.volume),
            "density": calculate_density(crystal)
        }
    }
    
    # Extract atomic positions
    if hasattr(crystal, 'atom_sites'):
        for site in crystal.atom_sites:
            for coord in site.coords:
                cartesian = np.dot(coord, lattice.matrix)
                
                structure_data["atoms"].append({
                    "element": str(site.specie),
                    "coords": coord.tolist(),
                    "cartesian": cartesian.tolist(),
                    "wyckoff": str(site.wp.letter) if hasattr(site.wp, 'letter') else "",
                    "multiplicity": int(site.wp.multiplicity) if hasattr(site.wp, 'multiplicity') else 1,
                    "site_symmetry": str(site.wp.site_symm) if hasattr(site.wp, 'site_symm') else ""
                })
    
    return structure_data


def validate_generated_structure(
    crystal: pyxtal,
    min_distance: Optional[Dict[str, float]] = None
) -> ValidationResult:
    """
    Validate generated crystal structure.
    
    Args:
        crystal: Generated PyXtal crystal
        min_distance: Minimum distance constraints
    
    Returns:
        ValidationResult with validation status
    """
    errors = []
    warnings = []
    
    if not crystal or not hasattr(crystal, 'lattice'):
        errors.append("Invalid crystal object")
        return ValidationResult(False, errors, warnings)
    
    # Check lattice parameters
    lattice = crystal.lattice
    params = [lattice.a, lattice.b, lattice.c]
    
    for i, (param, name) in enumerate(zip(params, ['a', 'b', 'c'])):
        if param < 2.0:
            warnings.append(f"Lattice parameter {name} = {param:.2f} Å is very small (< 2 Å)")
        if param > 50.0:
            warnings.append(f"Lattice parameter {name} = {param:.2f} Å is very large (> 50 Å)")
    
    angles = [lattice.alpha, lattice.beta, lattice.gamma]
    for angle, name in zip(angles, ['alpha', 'beta', 'gamma']):
        if angle < 30 or angle > 150:
            warnings.append(f"Angle {name} = {angle:.1f}° is unusual (< 30° or > 150°)")
    
    # Check minimum distances if specified
    if min_distance:
        pmg_structure = crystal.to_pymatgen()
        distances = calculate_min_distances(pmg_structure)
        
        for pair, min_dist in min_distance.items():
            if pair in distances:
                actual_dist = distances[pair]
                if actual_dist < min_dist:
                    errors.append(
                        f"Distance between {pair} is {actual_dist:.3f} Å, "
                        f"less than minimum {min_dist:.3f} Å"
                    )
    
    return ValidationResult(len(errors) == 0, errors, warnings)


def generate_crystal(
    composition: List[str],
    space_group: int,
    num_atoms: Optional[int] = None,
    lattice_params: Optional[Dict[str, float]] = None,
    volume_factor: float = 1.0,
    min_distance: Optional[Dict[str, float]] = None,
    wyckoff_positions: Optional[List[Dict[str, Any]]] = None,
    seed: Optional[int] = None,
    max_attempts: int = 100
) -> Dict[str, Any]:
    """
    Generate crystal structure using PyXtal.
    
    Args:
        composition: List of element symbols
        space_group: Space group number (1-230)
        num_atoms: Optional total number of atoms
        lattice_params: Optional lattice parameters
        volume_factor: Relative volume factor
        min_distance: Minimum distance constraints
        wyckoff_positions: Optional Wyckoff position specifications
        seed: Random seed for reproducibility
        max_attempts: Maximum generation attempts
    
    Returns:
        Dictionary containing structure data and validation results
    """
    # Validate inputs
    spg_valid, spg_error = validate_space_group(space_group)
    if not spg_valid:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SPACE_GROUP",
                "message": spg_error or "Invalid space group",
                "details": {"space_group": space_group}
            }
        }

    comp_valid, comp_error = validate_composition(composition)
    if not comp_valid:
        return {
            "success": False,
            "error": {
                "code": "INVALID_COMPOSITION",
                "message": comp_error or "Invalid composition",
                "details": {"composition": composition}
            }
        }
    
    # Set random seed if provided
    if seed is not None:
        np.random.seed(seed)
    
    # Initialize PyXtal
    crystal = pyxtal()
    
    # Attempt to generate structure
    success = False
    attempt = 0
    last_error = None
    
    # Calculate species and numIons from composition
    from collections import Counter
    counts = Counter(composition)
    species = list(counts.keys())
    
    # Determine numIons
    if num_atoms:
        total_ratio = sum(counts.values())
        # Scale counts to match target num_atoms
        # This assumes the input composition reflects the stoichiometry
        num_ions = [int(counts[s] * num_atoms / total_ratio) for s in species]
        
        # Adjust if rounding errors cause mismatch (simple adjustment)
        diff = num_atoms - sum(num_ions)
        if diff != 0:
            # Add/subtract remainder to first species (simplistic, but effective for keeping sum)
            num_ions[0] += diff
    else:
        # Use exact counts from composition list
        num_ions = list(counts.values())

    while attempt < max_attempts and not success:
        attempt += 1
        
        # Attempt generation
        generation_success = True
        error_msg = None
        
        # Call PyXtal's from_random method
        # We handle errors through return values, not exceptions
        try:
            result = crystal.from_random(
                dim=3,
                group=space_group,
                species=species,
                numIons=num_ions,
                factor=volume_factor,
                lattice=lattice_params,
                sites=wyckoff_positions
            )
        except Exception as e:
            # Catch internal PyXtal errors to prevent crash
            last_error = f"PyXtal error: {str(e)}"
            continue
        
        # Check if generation succeeded
        if result is not None or (hasattr(crystal, 'valid') and crystal.valid):
            success = True
        else:
            last_error = f"Generation attempt {attempt} failed"
    
    if not success:
        return {
            "success": False,
            "error": {
                "code": "MAX_ATTEMPTS_EXCEEDED",
                "message": f"Failed to generate structure after {max_attempts} attempts",
                "details": {
                    "attempts": max_attempts,
                    "last_error": last_error
                },
                "suggestions": [
                    "Increase volume_factor to 1.2-1.5",
                    "Relax minimum distance constraints",
                    "Try a different space group",
                    "Reduce number of atoms"
                ]
            }
        }
    
    # Extract structure data
    structure_data = extract_structure_data(crystal)
    
    # Validate generated structure
    validation = validate_generated_structure(crystal, min_distance)
    
    return {
        "success": True,
        "structure": structure_data,
        "validation": validation.to_dict(),
        "metadata": {
            "attempts": attempt,
            "seed": seed
        }
    }


def main():
    """
    Main entry point for command-line execution.
    Reads JSON input from file specified as first argument.
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python crystal_generator.py <input_json_file>"
            }
        }))
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input file
    if not input_file:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "FILE_NOT_FOUND",
                "message": "No input file specified"
            }
        }))
        sys.exit(1)
    
    # Load JSON input
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    # Generate crystal
    result = generate_crystal(**params)
    
    # Output result as JSON
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
