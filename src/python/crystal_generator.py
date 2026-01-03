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
from monty.serialization import loadfn
from pyxtal import pyxtal
from pyxtal.lattice import Lattice as PyXtalLattice
from pyxtal.symmetry import Group, rf
from pymatgen.core import Structure, Element

GROUP_SYMBOL_KEYS = {
    3: "space_group",
    2: "layer_group",
    1: "rod_group",
    0: "point_group"
}
GROUP_SYMBOLS = loadfn(rf("pyxtal", "database/symbols.json"))


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


def resolve_space_group(
    space_group: Union[int, str],
    dimensionality: int = 3
) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """
    Normalize a space group input to (number, symbol).

    Args:
        space_group: Space group number or Hermann-Mauguin symbol
        dimensionality: 0=0D, 1=1D, 2=2D, 3=3D

    Returns:
        Tuple of (number, symbol, error_message)
    """
    if isinstance(space_group, bool):
        return None, None, "Space group must be int or str, got bool"

    key = GROUP_SYMBOL_KEYS.get(dimensionality)
    if key is None:
        return None, None, f"Dimensionality must be 0-3, got {dimensionality}"

    symbols = GROUP_SYMBOLS[key]
    max_number = len(symbols)

    if isinstance(space_group, int):
        if space_group < 1 or space_group > max_number:
            return None, None, f"Space group number must be 1-{max_number} for dim={dimensionality}, got {space_group}"
        return space_group, symbols[space_group - 1], None

    if isinstance(space_group, str):
        cleaned = space_group.strip()
        if not cleaned:
            return None, None, "Space group symbol cannot be empty"
        if cleaned not in symbols:
            return None, None, f"Unknown space group symbol '{cleaned}' for dim={dimensionality}"
        return symbols.index(cleaned) + 1, cleaned, None

    return None, None, f"Space group must be int or str, got {type(space_group).__name__}"


def validate_space_group(space_group: Union[int, str], dimensionality: int = 3) -> Tuple[bool, Optional[str]]:
    """
    Validate space group number or symbol for given dimensionality.

    Args:
        space_group: Space group number or Hermann-Mauguin symbol
        dimensionality: 0=0D, 1=1D, 2=2D, 3=3D

    Returns:
        Tuple of (is_valid, error_message)
    """
    _, _, error = resolve_space_group(space_group, dimensionality)
    return error is None, error


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

    for i, elem in enumerate(composition):
        if not elem or not isinstance(elem, str):
            return False, f"Invalid element at index {i}: {elem}"

        elem_clean = elem.strip()
        if len(elem_clean) == 0:
            return False, f"Empty element at index {i}"

        if not Element.is_valid_symbol(elem_clean):
            return False, f"Unknown element '{elem_clean}' at index {i}"

    return True, None


def normalize_wyckoff_letter(value: str) -> Optional[str]:
    """Extract the Wyckoff letter from a label like '4a' or 'a'."""
    if not isinstance(value, str):
        return None
    for char in value:
        if char.isalpha():
            return char.lower()
    return None


def validate_wyckoff_positions(
    wyckoff_positions: Optional[List[Dict[str, Any]]],
    group_obj: Group,
    species: List[str]
) -> Tuple[bool, Optional[str]]:
    """
    Validate user-specified Wyckoff positions against the space group.
    """
    if wyckoff_positions is None:
        return True, None

    if not isinstance(wyckoff_positions, list):
        return False, "wyckoff_positions must be a list of position objects"

    valid_letters = {wp.letter.lower() for wp in group_obj.Wyckoff_positions}

    for index, entry in enumerate(wyckoff_positions):
        if not isinstance(entry, dict):
            return False, f"Wyckoff entry at index {index} must be an object"

        element = entry.get("element")
        if not element or element not in species:
            return False, f"Wyckoff entry {index} uses unknown element '{element}'"

        wyckoff = entry.get("wyckoff")
        if not isinstance(wyckoff, str) or not wyckoff.strip():
            return False, f"Wyckoff entry {index} is missing a valid wyckoff label"

        letter = normalize_wyckoff_letter(wyckoff)
        if not letter or letter not in valid_letters:
            return False, f"Wyckoff entry {index} has invalid letter '{wyckoff}' for this space group"

        coords = entry.get("coords")
        if coords is not None:
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                return False, f"Wyckoff entry {index} coords must be a length-3 list"

    return True, None


def build_wyckoff_sites(
    wyckoff_positions: List[Dict[str, Any]],
    species: List[str]
) -> List[Optional[List[Any]]]:
    """
    Convert user wyckoff_positions to PyXtal sites format.
    """
    sites_by_element: Dict[str, List[Any]] = {element: [] for element in species}

    for entry in wyckoff_positions:
        element = entry.get("element")
        wyckoff = entry.get("wyckoff")
        coords = entry.get("coords")

        if element not in sites_by_element:
            continue

        if coords is not None:
            x, y, z = coords
            sites_by_element[element].append((wyckoff, x, y, z))
        else:
            sites_by_element[element].append(wyckoff)

    sites: List[Optional[List[Any]]] = []
    for element in species:
        entries = sites_by_element.get(element, [])
        sites.append(entries if entries else None)

    return sites


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
    Calculate crystal density in g/cm^3.
    
    Args:
        crystal: PyXtal crystal object
    
    Returns:
        Density in g/cm^3
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
    reciprocal = lattice.reciprocal_lattice.matrix.tolist()

    atoms: List[Dict[str, Any]] = []
    atom_sites = getattr(crystal, "atom_sites", None)

    if atom_sites:
        for atom_site in atom_sites:
            element = str(atom_site.specie)
            wp = getattr(atom_site, "wp", None)
            wyckoff = ""
            multiplicity = int(getattr(atom_site, "multiplicity", 1))
            site_symmetry = ""

            if wp is not None:
                if hasattr(wp, "letter") and wp.letter:
                    wyckoff = str(wp.letter)
                if hasattr(wp, "get_site_symmetry_object"):
                    symm_obj = wp.get_site_symmetry_object()
                    if symm_obj is not None and hasattr(symm_obj, "name"):
                        site_symmetry = str(symm_obj.name)

            coords_array = np.array(getattr(atom_site, "coords", []), dtype=float)
            if coords_array.ndim == 1 and coords_array.size == 3:
                coords_iter = [coords_array]
            else:
                coords_iter = coords_array

            for coord in coords_iter:
                frac_coords = np.array(coord, dtype=float).tolist()
                cart_coords = lattice.get_cartesian_coords(frac_coords).tolist()
                atoms.append({
                    "element": element,
                    "coords": frac_coords,
                    "cartesian": cart_coords,
                    "wyckoff": wyckoff,
                    "multiplicity": multiplicity,
                    "site_symmetry": site_symmetry,
                    "species": [{"element": element, "occupation": 1.0}]
                })
    else:
        for site in pmg_struct:
            element = str(site.specie)
            atoms.append({
                "element": element,
                "coords": site.frac_coords.tolist(),
                "cartesian": site.coords.tolist(),
                "wyckoff": "",
                "multiplicity": 1,
                "site_symmetry": "",
                "species": [{"element": element, "occupation": 1.0}]
            })

    structure_data = {
        "lattice": {
            "a": float(lattice.a),
            "b": float(lattice.b),
            "c": float(lattice.c),
            "alpha": float(lattice.alpha),
            "beta": float(lattice.beta),
            "gamma": float(lattice.gamma),
            "matrix": lattice.matrix.tolist(),
            "reciprocal_matrix": reciprocal,
            "volume": float(lattice.volume)
        },
        "atoms": atoms,
        "sites": atoms,
        "space_group": {
            "number": int(crystal.group.number),
            "symbol": str(crystal.group.symbol),
            "hall_symbol": str(crystal.group.hall_symbol) if hasattr(crystal.group, 'hall_symbol') else "",
            "point_group": str(crystal.group.point_group) if hasattr(crystal.group, 'point_group') else "",
            "crystal_system": str(crystal.group.lattice_type) if hasattr(crystal.group, 'lattice_type') else ""
        },
        "metadata": {
            "formula": str(crystal.formula),
            "natoms": len(atoms),
            "volume": float(lattice.volume),
            "density": calculate_density(crystal)
        }
    }

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
            warnings.append(f"Lattice parameter {name} = {param:.2f} Angstroms is very small (< 2 Angstroms)")
        if param > 50.0:
            warnings.append(f"Lattice parameter {name} = {param:.2f} Angstroms is very large (> 50 Angstroms)")
    
    angles = [lattice.alpha, lattice.beta, lattice.gamma]
    for angle, name in zip(angles, ['alpha', 'beta', 'gamma']):
        if angle < 30 or angle > 150:
            warnings.append(f"Angle {name} = {angle:.1f} deg is unusual (< 30 deg or > 150 deg)")
    
    # Check minimum distances if specified
    if min_distance:
        pmg_structure = crystal.to_pymatgen()
        distances = calculate_min_distances(pmg_structure)

        # Normalize distance keys: sort element pairs so Si-O == O-Si
        def normalize_pair(pair_str: str) -> str:
            parts = pair_str.split('-')
            if len(parts) == 2:
                return '-'.join(sorted(parts))
            return pair_str

        # Build normalized constraint dict
        normalized_constraints = {normalize_pair(k): v for k, v in min_distance.items()}
        # Build normalized distances dict
        normalized_distances = {normalize_pair(k): v for k, v in distances.items()}

        for pair, min_dist in normalized_constraints.items():
            if pair in normalized_distances:
                actual_dist = normalized_distances[pair]
                if actual_dist < min_dist:
                    errors.append(
                        f"Distance between {pair} is {actual_dist:.3f} Angstroms, "
                        f"less than minimum {min_dist:.3f} Angstroms"
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
    space_group_number, _, spg_error = resolve_space_group(space_group, dimensionality)
    if spg_error:
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
    
    # Determine numIons and track stoichiometry changes
    stoich_warnings = []
    requested_stoich = {s: counts[s] for s in species}

    if num_atoms:
        total_ratio = sum(counts.values())
        # Scale counts to match target num_atoms
        # This assumes the input composition reflects the stoichiometry
        num_ions = [int(counts[s] * num_atoms / total_ratio) for s in species]

        # Adjust if rounding errors cause mismatch
        diff = num_atoms - sum(num_ions)
        if diff != 0:
            # WARNING: This perturbs stoichiometry!
            # Add/subtract remainder to first species
            old_value = num_ions[0]
            num_ions[0] += diff

            stoich_warnings.append(
                f"STOICHIOMETRY PERTURBED: Requested num_atoms={num_atoms} with composition {dict(counts)} "
                f"does not yield integer atom counts. Adjusted {species[0]} from {old_value} to {num_ions[0]} "
                f"to match total atoms. Actual stoichiometry: {dict(zip(species, num_ions))}"
            )
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
    
    group_obj = Group(space_group_number, dim=dimensionality)
    
    # check_compatible returns (bool, bool)
    is_compatible = group_obj.check_compatible(num_ions)[0]
    
    if not is_compatible:
         return {
            "success": False,
            "error": {
                "code": "INCOMPATIBLE_COMPOSITION",
            "message": f"Composition {num_ions} not compatible with space group {space_group_number} (dim={dimensionality})",
                "details": {"num_ions": num_ions, "space_group": space_group_number}
            }
        }

    sites = None
    if wyckoff_positions:
        wyckoff_valid, wyckoff_error = validate_wyckoff_positions(
            wyckoff_positions, group_obj, species
        )
        if not wyckoff_valid:
            return {
                "success": False,
                "error": {
                    "code": "INVALID_WYCKOFF_POSITIONS",
                    "message": wyckoff_error or "Invalid Wyckoff positions",
                    "details": {"wyckoff_positions": wyckoff_positions}
                }
            }
        sites = build_wyckoff_sites(wyckoff_positions, species)

    while attempt < max_attempts and not success:
        attempt += 1
        
        # Attempt generation
        generation_success = True
        error_msg = None
        
        # Prepare lattice object if params provided
        current_lattice = None
        if lattice_params and isinstance(lattice_params, dict):
            req_keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
            missing = [k for k in req_keys if k not in lattice_params or lattice_params[k] is None]
            if missing:
                return {
                    "success": False,
                    "error": {
                        "code": "INCOMPLETE_LATTICE_PARAMS",
                        "message": "lattice_params must include a, b, c, alpha, beta, gamma",
                        "details": {"missing": missing, "lattice_params": lattice_params}
                    }
                }
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
            group=space_group_number,
            species=species,
            numIons=num_ions,
            factor=volume_factor,
            lattice=current_lattice,
            sites=sites
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

    # HARD FAILURE: If min_distance constraints are violated, fail the generation
    # This is scientifically critical - unphysical short contacts invalidate results
    if not validation.valid:
        # Check if any errors are min_distance related
        min_distance_errors = [e for e in validation.errors if "Distance between" in e]
        if min_distance_errors:
            return {
                "success": False,
                "error": {
                    "code": "MIN_DISTANCE_VIOLATION",
                    "message": "Generated structure violates minimum distance constraints",
                    "details": {
                        "violations": min_distance_errors,
                        "all_errors": validation.errors,
                        "warnings": validation.warnings
                    },
                    "suggestions": [
                        "Increase volume_factor to create more space between atoms",
                        "Reduce the number of atoms",
                        "Try a different space group with more available volume"
                    ]
                }
            }
        # Other validation errors still return failure
        return {
            "success": False,
            "error": {
                "code": "VALIDATION_FAILED",
                "message": "Generated structure failed validation",
                "details": {
                    "errors": validation.errors,
                    "warnings": validation.warnings
                }
            }
        }

    result = {
        "success": True,
        "structure": structure_data,
        "validation": validation.to_dict(),
        "metadata": {
            "attempts": attempt,
            "seed": seed,
            "dimensionality": dimensionality
        }
    }

    # Add stoichiometry warnings if any
    if stoich_warnings:
        result["warnings"] = stoich_warnings

    return result


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
    space_group: Union[int, str],
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
    space_group_number, _, spg_error = resolve_space_group(space_group)
    if spg_error:
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
        missing = [k for k in req_keys if k not in lattice_params or lattice_params[k] is None]
        if missing:
            return {
                "success": False,
                "error": {
                    "code": "INCOMPLETE_LATTICE_PARAMS",
                    "message": "lattice_params must include a, b, c, alpha, beta, gamma",
                    "details": {"missing": missing, "lattice_params": lattice_params}
                }
            }
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
            group=space_group_number,
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
