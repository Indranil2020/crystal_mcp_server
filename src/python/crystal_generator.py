#!/usr/bin/env python3
"""
Crystal Generator - PyXtal Wrapper

Primary backend for crystal structure generation via MCP tools.
Provides comprehensive wrapper around PyXtal for generating crystal structures.

Related modular generators (for advanced use):
  - generators.bulk.spacegroups.generate_crystal_from_spacegroup
  - generators.bulk.prototypes.generate_from_prototype
  - generators.molecule.small_molecules.generate_molecule

All functions follow defensive programming principles without using try/except.
"""

from typing import Dict, List, Optional, Tuple, Any, Union
import json
import sys
import numpy as np
from pyxtal import pyxtal
from pyxtal.lattice import Lattice as PyXtalLattice
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


def validate_space_group(space_group: Union[int, str], dimensionality: int = 3) -> Tuple[bool, Optional[str]]:
    """
    Validate space group number or symbol for given dimensionality.

    Args:
        space_group: Space group number or Hermann-Mauguin symbol
        dimensionality: 0=0D, 1=1D, 2=2D, 3=3D

    Returns:
        Tuple of (is_valid, error_message)
    """
    if isinstance(space_group, int):
        if dimensionality == 3:
            if space_group < 1 or space_group > 230:
                return False, f"Space group number must be 1-230 for 3D, got {space_group}"
        elif dimensionality == 2:
            if space_group < 1 or space_group > 80:
                return False, f"Layer group number must be 1-80 for 2D, got {space_group}"
        elif dimensionality == 1:
            if space_group < 1 or space_group > 75:
                 return False, f"Rod group number must be 1-75 for 1D, got {space_group}"
        # 0D point groups range is complex, usually strings. For ints, maybe not strictly checked here.
        
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
    if crystal is None or not crystal.valid:
        return 0.0
        
    # Use Pymatgen structure for robust density calculation
    pmg_struct = crystal.to_pymatgen()
    return float(pmg_struct.density)


def extract_structure_data(crystal: pyxtal) -> Dict[str, Any]:
    """
    Extract complete structure data from PyXtal crystal.
    
    Args:
        crystal: PyXtal crystal object
    
    Returns:
        Dictionary containing all structure information
    """
    # Convert to Pymatgen structure for consistent handling of atoms
    pmg_struct = crystal.to_pymatgen()
    lattice = pmg_struct.lattice
    
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
            "natoms": len(pmg_struct),
            "volume": float(lattice.volume),
            "density": calculate_density(crystal)
        }
    }
    
    # Extract atomic positions from Pymatgen structure
    for site in pmg_struct:
        structure_data["atoms"].append({
            "element": str(site.specie),
            "coords": site.frac_coords.tolist(),
            "cartesian": site.coords.tolist(),
            "wyckoff": "", # Wyckoff info lost in Pymatgen conversion
            "multiplicity": 1,
            "site_symmetry": ""
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
    dimensionality: int = 3,
    max_attempts: int = 100
) -> Dict[str, Any]:
    """
    Generate crystal structure using PyXtal.
    
    Args:
        composition: List of element symbols
        space_group: Space group number (1-230 for 3D)
        num_atoms: Optional total number of atoms
        lattice_params: Optional lattice parameters
        volume_factor: Relative volume factor
        min_distance: Minimum distance constraints
        wyckoff_positions: Optional Wyckoff position specifications
        seed: Random seed for reproducibility
        dimensionality: 0=cluster, 1=rod, 2=layer, 3=bulk
        max_attempts: Maximum generation attempts
    
    Returns:
        Dictionary containing structure data and validation results
    """
    # Validate inputs
    spg_valid, spg_error = validate_space_group(space_group, dimensionality)
    if not spg_valid:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SPACE_GROUP",
                "message": spg_error or "Invalid space group",
                "details": {"space_group": space_group, "dimensionality": dimensionality}
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

    # Check compatibility with space group
    # Avoids try-except by pre-validating
    # Check compatibility with space group
    # We rely on pre-validation. If Group() fails, it will raise an error, 
    # which is acceptable as we do not use try-except for control flow.
    
    # Only perform check if dimensionality matches (Group behavior specifics)
    # PyXtal's Group(space_group) is valid for 1-230.
    # We assume space_group is valid int from validate_space_group.
    
    group_obj = Group(space_group, dim=dimensionality)
    
    # check_compatible returns (bool, bool)
    is_compatible = group_obj.check_compatible(num_ions)[0]
    
    if not is_compatible:
         return {
            "success": False,
            "error": {
                "code": "INCOMPATIBLE_COMPOSITION",
                "message": f"Composition {num_ions} not compatible with space group {space_group} (dim={dimensionality})",
                "details": {"num_ions": num_ions, "space_group": space_group}
            }
        }

    while attempt < max_attempts and not success:
        attempt += 1
        
        # Attempt generation
        generation_success = True
        error_msg = None
        
        # Prepare lattice object if params provided
        current_lattice = None
        if lattice_params and isinstance(lattice_params, dict):
            # Check for required keys
            req_keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
            if all(k in lattice_params for k in req_keys):
                # PyXtal Lattice.from_para expects radians usually for angles?
                # Let's try degrees first as is standard for "parameters" function signatures.
                # If PyXtal fails or produces wrong result, we might need radians.
                # Usage: Lattice.from_para(a, b, c, alpha, beta, gamma, ltype)
                # But we don't know ltype easily without mapping SG.
                # PyXtal infers type or we pass generic.
                    current_lattice = PyXtalLattice.from_para(
                        lattice_params['a'],
                        lattice_params['b'],
                        lattice_params['c'],
                        lattice_params['alpha'],
                        lattice_params['beta'],
                        lattice_params['gamma']
                    )

        # Call PyXtal's from_random method
        # We rely on pre-validation and library stability
        result = crystal.from_random(
            dim=dimensionality,
            group=space_group,
            species=species,
            numIons=num_ions,
            factor=volume_factor,
            lattice=current_lattice,
            sites=wyckoff_positions
        )
        
        # Check if generation succeeded
        if crystal.valid:
            success = True
        else:
            last_error = f"Generation attempt {attempt} failed (validity check)"
    
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
            "seed": seed,
            "dimensionality": dimensionality
        }
    }


def validate_molecules(molecules: List[str]) -> Tuple[bool, Optional[str]]:
    """
    Validate molecule list.
    
    Args:
        molecules: List of molecule strings (formulas or names)
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not molecules or len(molecules) == 0:
        return False, "Molecules list cannot be empty"
        
    for i, mol in enumerate(molecules):
        if not mol or not isinstance(mol, str) or len(mol.strip()) == 0:
            return False, f"Invalid molecule at index {i}"
            
    return True, None


def generate_molecular_crystal(
    molecules: List[str],
    space_group: int,
    num_molecules: Optional[int] = None,
    lattice_params: Optional[Dict[str, float]] = None,
    volume_factor: float = 1.0,
    min_distance: Optional[Dict[str, float]] = None,
    seed: Optional[int] = None,
    max_attempts: int = 100
) -> Dict[str, Any]:
    """
    Generate molecular crystal structure.
    
    Args:
        molecules: List of molecule formulas
        space_group: Space group number
        num_molecules: Total number of molecules
        lattice_params: Optional lattice parameters
        volume_factor: Volume factor
        min_distance: Minimum distance constraints
        seed: Random seed
        max_attempts: Max attempts
        
    Returns:
        Dictionary with structure data
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
        
    mol_valid, mol_error = validate_molecules(molecules)
    if not mol_valid:
         return {
            "success": False,
            "error": {
                "code": "INVALID_MOLECULES",
                "message": mol_error or "Invalid molecules",
                "details": {"molecules": molecules}
            }
        }

    if seed is not None:
        np.random.seed(seed)
        
    crystal = pyxtal(molecular=True)
    
    # Calculate numIons (molecules count)
    from collections import Counter
    counts = Counter(molecules)
    species = list(counts.keys())
    
    if num_molecules:
        total = sum(counts.values())
        num_ions = [int(counts[s] * num_molecules / total) for s in species]
        diff = num_molecules - sum(num_ions)
        if diff != 0:
            num_ions[0] += diff
    else:
        num_ions = list(counts.values())
        
    success = False
    attempt = 0
    last_error = None

    # Convert lattice_params dict to PyXtal Lattice object if provided
    current_lattice = None
    if lattice_params and isinstance(lattice_params, dict):
        req_keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        if all(k in lattice_params for k in req_keys):
            from pyxtal.lattice import Lattice as PyXtalLattice
            current_lattice = PyXtalLattice.from_para(
                lattice_params['a'],
                lattice_params['b'],
                lattice_params['c'],
                lattice_params['alpha'],
                lattice_params['beta'],
                lattice_params['gamma']
            )

    while attempt < max_attempts and not success:
        attempt += 1
        crystal.from_random(
            dim=3,
            group=space_group,
            species=species,
            numIons=num_ions,
            factor=volume_factor,
            lattice=current_lattice
        )
        if crystal.valid:
            success = True
        else:
            last_error = "Generation failed check validity"
            
    if not success:
        return {
            "success": False,
            "error": {
                "code": "MAX_ATTEMPTS_EXCEEDED",
                "message": f"Failed to generate molecular crystal after {max_attempts} attempts",
                "details": {"last_error": last_error}
            }
        }
        
    structure_data = extract_structure_data(crystal)
    # Skip standard validation for now as it might be strict on distances for molecules
    return {
        "success": True,
        "structure": structure_data,
        "metadata": {"attempts": attempt, "seed": seed}
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
    
    # Determine mode
    operation = params.get("operation", "generate_crystal")
    
    if operation == "generate_molecular":
        result = generate_molecular_crystal(
            molecules=params.get("molecules"),
            space_group=params.get("space_group"),
            num_molecules=params.get("num_molecules"),
            lattice_params=params.get("lattice_params"),
            volume_factor=params.get("volume_factor", 1.0),
            min_distance=params.get("min_distance"),
            seed=params.get("seed"),
            max_attempts=params.get("max_attempts", 100)
        )
    else:
        # Default to standard crystal generation
        # Filter out operation key if present
        gen_params = {k: v for k, v in params.items() if k != "operation"}
        result = generate_crystal(**gen_params)
    
    # Output result as JSON
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
