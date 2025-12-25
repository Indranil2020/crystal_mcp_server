#!/usr/bin/env python3
"""
Structure Tools - Transformation Utilities

This module provides structure transformation capabilities including
supercell generation, slab creation, defect introduction, and strain application.
All functions follow defensive programming without try/except blocks.
"""

from typing import Dict, List, Optional, Any, Tuple
import json
import sys
import numpy as np
import spglib
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator

CRYSTAL_SYSTEM_RANGES = [
    (1, 2, "triclinic"),
    (3, 15, "monoclinic"),
    (16, 74, "orthorhombic"),
    (75, 142, "tetragonal"),
    (143, 167, "trigonal"),
    (168, 194, "hexagonal"),
    (195, 230, "cubic")
]


def crystal_system_from_number(number: Optional[int]) -> str:
    """Map space group number to crystal system."""
    if number is None:
        return ""
    for start, end, name in CRYSTAL_SYSTEM_RANGES:
        if start <= number <= end:
            return name
    return ""


def get_space_group_info(structure: Structure, symprec: float = 1e-3) -> Dict[str, Any]:
    """Extract space group info using spglib dataset."""
    if structure is None or len(structure) == 0:
        return {"number": None, "symbol": "", "hall_symbol": "", "point_group": "", "crystal_system": ""}

    cell = (
        structure.lattice.matrix,
        structure.frac_coords,
        structure.atomic_numbers
    )
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
    if dataset is None:
        return {"number": None, "symbol": "", "hall_symbol": "", "point_group": "", "crystal_system": ""}

    number = dataset.get("number")
    number_int = int(number) if number is not None else None
    symbol = dataset.get("international", "") or ""
    hall_symbol = dataset.get("hall", "") or ""
    point_group = dataset.get("pointgroup", "") or ""

    return {
        "number": number_int,
        "symbol": symbol,
        "hall_symbol": hall_symbol,
        "point_group": point_group,
        "crystal_system": crystal_system_from_number(number_int)
    }


def apply_deformation(structure: Structure, deformation: np.ndarray) -> Structure:
    """Apply a lattice deformation matrix to a structure."""
    new_matrix = np.dot(deformation, structure.lattice.matrix)
    new_lattice = Lattice(new_matrix)
    return Structure(new_lattice, structure.species, structure.frac_coords, coords_are_cartesian=False)


def validate_supercell_matrix(scaling: Any) -> Tuple[bool, Optional[str]]:
    """
    Validate a supercell scaling matrix or vector.

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
    if len(scaling) == 3:
        if all(isinstance(x, (int, float)) for x in scaling):
            if any(x <= 0 for x in scaling):
                return False, "Scaling factors must be positive"
            if any(x != int(x) for x in scaling):
                return False, "Scaling factors must be integers"
            return True, None

        # Could be a 3x3 matrix
        if all(isinstance(row, (list, tuple)) and len(row) == 3 for row in scaling):
            matrix = np.array(scaling)
            if matrix.shape != (3, 3):
                return False, "Invalid matrix shape"
            det = np.linalg.det(matrix)
            if abs(det) < 1e-6:
                return False, "Scaling matrix is singular (determinant ~ 0)"
            return True, None

    return False, "Scaling must be [nx, ny, nz] or 3x3 matrix"


def calculate_supercell_size(scaling: Any) -> int:
    """
    Calculate the multiplier for a supercell scaling matrix.

    Args:
        scaling: Either [nx, ny, nz] vector or 3x3 matrix

    Returns:
        Integer multiplier for number of atoms
    """
    if not isinstance(scaling, (list, tuple, np.ndarray)):
        return 1

    if isinstance(scaling, np.ndarray):
        scaling = scaling.tolist()

    # If it's a vector [nx, ny, nz], multiply them
    if len(scaling) == 3 and all(isinstance(x, (int, float)) for x in scaling):
        return int(scaling[0] * scaling[1] * scaling[2])

    # If it's a 3x3 matrix, use determinant
    if len(scaling) == 3:
        if all(isinstance(row, (list, tuple)) and len(row) == 3 for row in scaling):
            matrix = np.array(scaling)
            det = abs(np.linalg.det(matrix))
            return int(round(det))

    return 1


def dict_to_structure(structure_dict: Dict[str, Any]) -> Optional[Structure]:
    """
    Convert structure dictionary to Pymatgen Structure.
    
    Args:
        structure_dict: Structure dictionary
    
    Returns:
        Pymatgen Structure or None
    """
    if not structure_dict or not isinstance(structure_dict, dict):
        return None
    
    if "lattice" not in structure_dict or ("atoms" not in structure_dict and "sites" not in structure_dict):
        return None
    
    lattice_data = structure_dict["lattice"]
    if "matrix" not in lattice_data:
        return None
    
    lattice_matrix = lattice_data["matrix"]
    lattice = Lattice(lattice_matrix)
    
    atoms_data = structure_dict.get("atoms") or structure_dict.get("sites")
    if not atoms_data:
        return None
    
    species = []
    frac_coords = []
    
    for atom in atoms_data:
        element = atom.get("element")
        if element is None and "species" in atom:
            species_list = atom.get("species")
            if isinstance(species_list, list) and species_list:
                element = species_list[0].get("element")

        atom_coords = atom.get("coords") or atom.get("abc")
        if atom_coords is None and "cartesian" in atom:
            atom_coords = lattice.get_fractional_coords(atom["cartesian"])

        if element is None or atom_coords is None:
            return None
        species.append(element)
        frac_coords.append(atom_coords)
    
    structure = Structure(lattice, species, frac_coords, coords_are_cartesian=False)
    return structure


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """
    Convert Pymatgen Structure to dictionary.
    
    Args:
        structure: Pymatgen Structure
    
    Returns:
        Structure dictionary
    """
    lattice = structure.lattice
    
    reciprocal = lattice.reciprocal_lattice.matrix.tolist()
    space_group = get_space_group_info(structure)

    sites = []
    for site in structure:
        element = str(site.specie)
        sites.append({
            "element": element,
            "coords": site.frac_coords.tolist(),
            "cartesian": site.coords.tolist(),
            "species": [{"element": element, "occupation": 1.0}]
        })

    result = {
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
        "atoms": sites,
        "sites": sites,
        "space_group": space_group,
        "metadata": {
            "natoms": len(structure),
            "formula": structure.composition.formula
        }
    }

    return result


def make_supercell(
    structure_dict: Dict[str, Any],
    scaling_matrix: List[List[int]],
    wrap_atoms: bool = True,
    preserve_symmetry: bool = False,
    min_distance_check: bool = True
) -> Dict[str, Any]:
    """
    Create supercell from structure.
    
    Args:
        structure_dict: Original structure
        scaling_matrix: Scaling matrix (3x3) or [nx, ny, nz]
        wrap_atoms: Wrap atoms back into cell
        preserve_symmetry: Try to preserve space group
        min_distance_check: Check for overlapping atoms
    
    Returns:
        Dictionary with supercell structure
    """
    # Convert to Pymatgen Structure
    structure = dict_to_structure(structure_dict)
    if structure is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure to Pymatgen Structure",
                "details": {}
            }
        }
    
    # Handle presets
    if isinstance(scaling_matrix, str):
        presets = {
            "sqrt3": [[2, 1, 0], [-1, 1, 0], [0, 0, 1]], # Det=3, sqrt(3)xsqrt(3)R30 for hexagonal
            "root2": [[1, 1, 0], [1, -1, 0], [0, 0, 1]], # Det=2, sqrt(2)xsqrt(2)R45 for square
            "2x2x2": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
            "3x3x3": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        }
        if scaling_matrix.lower() in presets:
            scaling_matrix = presets[scaling_matrix.lower()]
        else:
             return {
                "success": False, 
                "error": {
                    "code": "INVALID_SCALING_MATRIX",
                    "message": f"Unknown matrix preset: {scaling_matrix}",
                    "details": {"available": list(presets.keys())}
                }
            }

    # Validate scaling matrix
    if isinstance(scaling_matrix, list):
        if len(scaling_matrix) == 3 and all(isinstance(x, (int, float)) for x in scaling_matrix):
            # Convert [nx, ny, nz] to 3x3 diagonal matrix
            scaling_matrix = [
                [int(scaling_matrix[0]), 0, 0],
                [0, int(scaling_matrix[1]), 0],
                [0, 0, int(scaling_matrix[2])]
            ]
        elif len(scaling_matrix) != 3 or any(len(row) != 3 for row in scaling_matrix):
            return {
                "success": False,
                "error": {
                    "code": "INVALID_SCALING_MATRIX",
                    "message": "Scaling matrix must be 3x3 or [nx, ny, nz]",
                    "details": {"scaling_matrix": scaling_matrix}
                }
            }
    else:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SCALING_MATRIX",
                "message": "Scaling matrix must be a list",
                "details": {"type": type(scaling_matrix).__name__}
            }
        }
    
    # Calculate volume multiplier
    scaling_array = np.array(scaling_matrix)
    volume_multiplier = float(np.abs(np.linalg.det(scaling_array)))
    
    if volume_multiplier < 1e-6:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SCALING_MATRIX",
                "message": "Scaling matrix has zero determinant",
                "details": {"determinant": volume_multiplier}
            }
        }
    
    # Create supercell
    supercell = structure.copy()
    supercell.make_supercell(scaling_matrix)
    
    # Wrap atoms if requested
    if wrap_atoms:
        for i, site in enumerate(supercell):
            frac_coords = site.frac_coords
            wrapped_coords = frac_coords - np.floor(frac_coords)
            supercell.translate_sites(i, wrapped_coords - frac_coords)
    
    # Check minimum distances if requested
    if min_distance_check:
        min_dist = float('inf')
        for i in range(len(supercell)):
            for j in range(i + 1, len(supercell)):
                dist = supercell.get_distance(i, j)
                if dist < min_dist:
                    min_dist = dist
        
        if min_dist < 0.5:  # Less than 0.5 Angstrom
            return {
                "success": False,
                "error": {
                    "code": "ATOMS_TOO_CLOSE",
                    "message": f"Minimum distance {min_dist:.3f} Angstroms is too small",
                    "details": {"min_distance": min_dist},
                    "suggestion": "Check scaling matrix or original structure"
                }
            }
    
    # Convert to dictionary
    supercell_dict = structure_to_dict(supercell)
    original_info = get_space_group_info(structure)
    new_info = supercell_dict.get("space_group", {})
    preserved_symmetry = (
        original_info.get("number") is not None
        and original_info.get("number") == new_info.get("number")
    )

    return {
        "success": True,
        "supercell": supercell_dict,
        "transformation_matrix": scaling_array.tolist(),
        "original_space_group": original_info.get("number"),
        "new_space_group": new_info.get("number"),
        "volume_multiplier": volume_multiplier,
        "preserved_symmetry": preserved_symmetry
    }


def generate_slab(
    structure_dict: Dict[str, Any],
    miller_indices: List[int],
    thickness: int,
    vacuum: float,
    center_slab: bool = True,
    fix_bottom_layers: int = 0,
    orthogonalize: bool = False,
    min_slab_size: float = 0.0,
    symmetric: bool = True
) -> Dict[str, Any]:
    """
    Generate surface slab from bulk structure.
    
    Args:
        structure_dict: Bulk structure
        miller_indices: Miller indices [h, k, l]
        thickness: Number of layers (minimum)
        vacuum: Vacuum thickness in Angstroms
        center_slab: Center slab in cell
        fix_bottom_layers: Number of bottom layers to fix
        orthogonalize: Create orthogonal cell
        min_slab_size: Minimum slab size in Angstroms
        symmetric: Make symmetric slab
    
    Returns:
        Dictionary with slab structure
    """
    # Convert to Pymatgen Structure
    structure = dict_to_structure(structure_dict)
    if structure is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure",
                "details": {}
            }
        }
    
    # Validate miller indices
    if not miller_indices or len(miller_indices) != 3:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MILLER_INDICES",
                "message": "Miller indices must be [h, k, l]",
                "details": {"miller_indices": miller_indices}
            }
        }
    
    if all(x == 0 for x in miller_indices):
        return {
            "success": False,
            "error": {
                "code": "INVALID_MILLER_INDICES",
                "message": "Miller indices cannot be [0, 0, 0]",
                "details": {}
            }
        }
    
    # Validate parameters
    if thickness < 1:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"Thickness must be >= 1 layer, got {thickness}",
                "details": {"thickness": thickness}
            }
        }

    if vacuum < 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"Vacuum must be >= 0 Angstroms, got {vacuum}",
                "details": {"vacuum": vacuum}
            }
        }

    # CRITICAL: thickness is in LAYERS, must convert to Angstroms
    # Calculate d-spacing for the given Miller indices
    lattice = structure.lattice
    h, k, l = miller_indices

    # For cubic/orthorhombic: simple formula
    # For general case: use reciprocal lattice
    recip_lattice = lattice.reciprocal_lattice
    # d-spacing = 1 / |G| where G is reciprocal lattice vector (h,k,l)
    G_vector = h * recip_lattice.matrix[0] + k * recip_lattice.matrix[1] + l * recip_lattice.matrix[2]
    d_spacing = 1.0 / np.linalg.norm(G_vector) if np.linalg.norm(G_vector) > 1e-10 else 3.0

    # Convert layers to Angstroms
    thickness_angstrom = thickness * d_spacing

    # Use the larger of calculated thickness or provided min_slab_size
    actual_min_slab_size = max(min_slab_size, thickness_angstrom) if min_slab_size is not None else thickness_angstrom

    # Generate slab using Pymatgen
    slabgen = SlabGenerator(
        structure,
        miller_indices,
        min_slab_size=actual_min_slab_size,
        min_vacuum_size=vacuum,
        center_slab=center_slab,
        primitive=False,
        max_normal_search=1
    )
    
    slabs = slabgen.get_slabs(symmetrize=symmetric)
    
    if not slabs or len(slabs) == 0:
        return {
            "success": False,
            "error": {
                "code": "SLAB_GENERATION_FAILED",
                "message": "Failed to generate slab",
                "details": {
                    "miller_indices": miller_indices,
                    "thickness": thickness
                }
            }
        }
    
    # Use first slab
    slab = slabs[0]
    
    # Orthogonalize if requested
    if orthogonalize:
        # Pymatgen's get_orthogonal_c_slab might raise errors if impossible
        # We allow it to propagate
        slab = slab.get_orthogonal_c_slab()
    
    # Calculate fixed atoms
    fixed_atoms = []
    if fix_bottom_layers > 0:
        # Sort atoms by z-coordinate
        z_coords = [site.coords[2] for site in slab]
        sorted_indices = np.argsort(z_coords)
        
        # Calculate atoms per layer (approximately)
        n_atoms_per_layer = len(slab) // max(thickness, 1)
        n_fixed = min(fix_bottom_layers * n_atoms_per_layer, len(slab) // 2)
        
        fixed_atoms = sorted_indices[:n_fixed].tolist()
    
    # Convert to dictionary
    slab_dict = structure_to_dict(slab)
    original_info = get_space_group_info(structure)
    
    # Calculate surface area
    a = slab.lattice.a
    b = slab.lattice.b
    gamma = slab.lattice.gamma * np.pi / 180
    surface_area = float(a * b * np.sin(gamma))
    
    return {
        "success": True,
        "slab": slab_dict,
        "surface_area": surface_area,
        "slab_thickness": float(slab.lattice.c - vacuum),
        "vacuum_thickness": float(vacuum),
        "termination": "unknown",  # Would require surface analysis
        "miller_indices": miller_indices,
        "fixed_atoms": fixed_atoms,
        "metadata": {
            "original_space_group": original_info.get("number"),
            "slab_space_group": slab_dict.get("space_group", {}).get("number"),
            "n_original_atoms": len(structure),
            "n_slab_atoms": len(slab)
        }
    }


def create_defect(
    structure_dict: Dict[str, Any],
    defect_type: str,
    defect_site: int,
    defect_species: Optional[str] = None,
    concentration: float = 1.0
) -> Dict[str, Any]:
    """
    Create point defect in structure.
    
    Args:
        structure_dict: Original structure
        defect_type: Type of defect (vacancy, substitution, interstitial)
        defect_site: Atom index for defect
        defect_species: Element for substitution
        concentration: Defect concentration (future use)
    
    Returns:
        Dictionary with defected structure
    """
    # Convert to Pymatgen Structure
    structure = dict_to_structure(structure_dict)
    if structure is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure",
                "details": {}
            }
        }
    
    # Validate defect type
    valid_types = ["vacancy", "substitution", "interstitial"]
    if defect_type not in valid_types:
        return {
            "success": False,
            "error": {
                "code": "INVALID_DEFECT_TYPE",
                "message": f"Defect type must be one of: {', '.join(valid_types)}",
                "details": {"defect_type": defect_type}
            }
        }
    
    # Validate defect site
    if defect_site < 0 or defect_site >= len(structure):
        return {
            "success": False,
            "error": {
                "code": "INVALID_DEFECT_SITE",
                "message": f"Defect site {defect_site} out of range [0, {len(structure)-1}]",
                "details": {"defect_site": defect_site, "natoms": len(structure)}
            }
        }
    
    # Create defected structure
    defected = structure.copy()
    
    if defect_type == "vacancy":
        defected.remove_sites([defect_site])
    
    elif defect_type == "substitution":
        if not defect_species:
            return {
                "success": False,
                "error": {
                    "code": "MISSING_DEFECT_SPECIES",
                    "message": "Substitution requires defect_species",
                    "details": {}
                }
            }
        defected.replace(defect_site, defect_species)
    
    elif defect_type == "interstitial":
        if not defect_species:
            return {
                "success": False,
                "error": {
                    "code": "MISSING_DEFECT_SPECIES",
                    "message": "Interstitial requires defect_species",
                    "details": {}
                }
            }
        # Add interstitial at a low-overlap position (grid search)
        grid = np.linspace(0.1, 0.9, 4)
        candidates = np.array([(x, y, z) for x in grid for y in grid for z in grid])
        distances = structure.lattice.get_all_distances(candidates, structure.frac_coords)
        min_distances = distances.min(axis=1)
        best_index = int(np.argmax(min_distances))
        new_coords = candidates[best_index]
        defected.append(defect_species, new_coords)
    
    # Convert to dictionary
    defected_dict = structure_to_dict(defected)
    
    return {
        "success": True,
        "defected_structure": defected_dict,
        "defect_type": defect_type,
        "defect_site": defect_site,
        "defect_species": defect_species
    }


def apply_strain(
    structure_dict: Dict[str, Any],
    strain_tensor: Optional[List] = None,
    strain_type: Optional[str] = None,
    strain_value: Optional[float] = None
) -> Dict[str, Any]:
    """
    Apply strain to structure.

    Args:
        structure_dict: Original structure
        strain_tensor: Strain tensor as 3x3 matrix or flattened 9-element array
        strain_type: Predefined strain type ("biaxial", "uniaxial", "hydrostatic")
        strain_value: Strain magnitude (e.g., 0.05 for 5% strain)

    Returns:
        Dictionary with strained structure
    """
    # Convert to Pymatgen Structure
    structure = dict_to_structure(structure_dict)
    if structure is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure",
                "details": {}
            }
        }

    # Build strain tensor from inputs
    if strain_tensor is not None:
        tensor = np.array(strain_tensor)
        # Handle flattened 9-element array
        if tensor.ndim == 1 and len(tensor) == 9:
            tensor = tensor.reshape(3, 3)
        elif tensor.ndim != 2 or tensor.shape != (3, 3):
            return {
                "success": False,
                "error": {
                    "code": "INVALID_STRAIN_TENSOR",
                    "message": "Strain tensor must be 3x3 matrix or 9-element array",
                    "details": {"shape": list(tensor.shape)}
                }
            }
    elif strain_type is not None and strain_value is not None:
        # Build tensor from strain_type and strain_value
        s = strain_value
        if strain_type == "hydrostatic":
            tensor = np.array([[s, 0, 0], [0, s, 0], [0, 0, s]])
        elif strain_type == "biaxial":
            tensor = np.array([[s, 0, 0], [0, s, 0], [0, 0, 0]])
        elif strain_type == "uniaxial":
            tensor = np.array([[s, 0, 0], [0, 0, 0], [0, 0, 0]])
        else:
            return {
                "success": False,
                "error": {
                    "code": "INVALID_STRAIN_TYPE",
                    "message": f"Unknown strain type: {strain_type}",
                    "details": {"valid_types": ["biaxial", "uniaxial", "hydrostatic"]}
                }
            }
    else:
        return {
            "success": False,
            "error": {
                "code": "MISSING_STRAIN_PARAMS",
                "message": "Must provide either strain_tensor or (strain_type and strain_value)",
                "details": {}
            }
        }

    # Validate strain tensor values
    warnings = []
    max_strain = float(np.max(np.abs(tensor)))

    # Check for unrealistic strain values (> 50% is typically unphysical)
    if max_strain > 0.5:
        return {
            "success": False,
            "error": {
                "code": "STRAIN_TOO_LARGE",
                "message": f"Strain magnitude {max_strain:.1%} exceeds 50% - this is typically unphysical",
                "details": {"max_strain": max_strain, "limit": 0.5}
            }
        }

    if max_strain > 0.1:
        warnings.append(f"Large strain magnitude detected: {max_strain:.1%} (> 10%)")
    elif max_strain > 0.05:
        warnings.append(f"Moderate strain magnitude: {max_strain:.1%} (> 5%)")

    # Check for asymmetry in tensor (shear strains)
    if not np.allclose(tensor, tensor.T, atol=1e-10):
        warnings.append("Strain tensor is asymmetric - includes rotational component")

    # Apply strain using full deformation matrix
    deformation = np.eye(3) + tensor
    strained = apply_deformation(structure, deformation)

    # Convert to dictionary
    strained_dict = structure_to_dict(strained)

    result = {
        "success": True,
        "strained_structure": strained_dict,
        "strain_tensor": tensor.tolist(),
        "deformation_matrix": deformation.tolist()
    }
    if warnings:
        result["warnings"] = warnings

    return result


def create_alloy(
    structure_dict: Dict[str, Any],
    substitutions: Dict[str, Dict[str, Any]],
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Create an alloy by random substitution.
    
    Args:
        structure_dict: Base structure
        substitutions: Dict mapping element -> {element, concentration}
                       e.g. {"Si": {"element": "Ge", "concentration": 0.5}}
        seed: Random seed
    
    Returns:
        Alloy structure
    """
    if seed is not None:
        np.random.seed(seed)
        
    structure = dict_to_structure(structure_dict)
    if structure is None:
        return {"success": False, "error": "Invalid structure"}
        
    alloy = structure.copy()
    
    # Process each substitution rule
    for target_elem, rule in substitutions.items():
        if "element" not in rule or "concentration" not in rule:
             continue
             
        new_elem = rule["element"]
        conc = float(rule["concentration"])
        
        # Find sites with target element
        # Only check symbol, ignore oxidation state for simplicity
        indices = [i for i, site in enumerate(alloy) if site.specie.symbol == target_elem]
        
        if not indices:
            continue
            
        n_replace = int(round(len(indices) * conc))
        if n_replace == 0:
            continue
            
        # Randomly select indices to replace
        replace_indices = np.random.choice(indices, size=n_replace, replace=False)
        
        for idx in replace_indices:
            alloy.replace(idx, new_elem)
            
    return {
        "success": True,
        "alloy_structure": structure_to_dict(alloy)
    }


def create_heterostructure(
    substrate_dict: Dict[str, Any],
    overlayer_dict: Dict[str, Any],
    interface_distance: float = 3.0,
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Create a simple vertical heterostructure (stacking).
    Does NOT perform advanced lattice matching analysis yet (assumes pre-matched).
    
    Args:
        substrate_dict: Substrate structure
        overlayer_dict: Overlayer structure
        interface_distance: Distance between layers
        vacuum: Vacuum padding
    
    Returns:
        Heterostructure
    """
    sub = dict_to_structure(substrate_dict)
    over = dict_to_structure(overlayer_dict)
    
    if sub is None or over is None:
        return {"success": False, "error": "Invalid input structures"}
        
    # Validation: Check lattice match (a, b)
    # This is a 'dumb' stacker. It assumes user knows lattices match.
    # We can warn if they differ significantly.
    
    # Get lattice vectors
    a1, b1, _ = sub.lattice.abc
    a2, b2, _ = over.lattice.abc
    
    mismatch_a = abs(a1 - a2) / a1
    mismatch_b = abs(b1 - b2) / b1
    
    warnings = []
    if mismatch_a > 0.05 or mismatch_b > 0.05:
        warnings.append(f"High lattice mismatch: a={mismatch_a:.1%}, b={mismatch_b:.1%}")
    else:
        scale_a = a1 / a2
        scale_b = b1 / b2
        if abs(scale_a - 1.0) > 1e-3 or abs(scale_b - 1.0) > 1e-3:
            deformation = np.array([[scale_a, 0, 0], [0, scale_b, 0], [0, 0, 1.0]])
            over = apply_deformation(over, deformation)
            warnings.append(f"Applied in-plane scaling to overlayer: scale_a={scale_a:.4f}, scale_b={scale_b:.4f}")
            a2, b2, _ = over.lattice.abc
        
    # Stack layers properly for non-orthogonal cells
    # IMPORTANT: Do NOT assume c-axis is aligned with Cartesian z-axis!
    # Instead, work with projections onto the c lattice vector direction.

    # Get c-axis direction (normalized)
    sub_lattice_matrix = sub.lattice.matrix
    c_vector = sub_lattice_matrix[2]  # Third lattice vector
    c_unit = c_vector / np.linalg.norm(c_vector)

    # Project substrate coordinates onto c-direction
    sub_coords = sub.cart_coords
    sub_proj = sub_coords @ c_unit  # Projection onto c-direction
    max_proj_sub = np.max(sub_proj) if len(sub) > 0 else 0
    min_proj_sub = np.min(sub_proj) if len(sub) > 0 else 0
    sub_height = max_proj_sub - min_proj_sub

    # Project overlayer coordinates onto c-direction
    over_coords = over.cart_coords
    over_proj = over_coords @ c_unit
    max_proj_over = np.max(over_proj) if len(over) > 0 else 0
    min_proj_over = np.min(over_proj) if len(over) > 0 else 0
    over_height = max_proj_over - min_proj_over

    # Calculate shift along c-direction
    # Place overlayer bottom at: substrate_top + interface_distance
    shift_distance = (max_proj_sub + interface_distance) - min_proj_over
    shift_vector = shift_distance * c_unit

    # Shift overlayer coordinates
    new_over_coords = over_coords + shift_vector

    # Calculate new c lattice parameter
    # new_c = substrate_height + interface_distance + overlayer_height + vacuum
    new_c_length = sub_height + interface_distance + over_height + vacuum

    # Create new lattice with same a, b vectors but scaled c vector
    # This preserves non-orthogonal angles
    c_scale = new_c_length / np.linalg.norm(c_vector)
    new_c_vector = c_vector * c_scale

    new_lattice_matrix = np.array([
        sub_lattice_matrix[0],  # Keep a vector
        sub_lattice_matrix[1],  # Keep b vector
        new_c_vector            # Scale c vector
    ])

    new_lattice = Lattice(new_lattice_matrix)
    
    # Collect all species and coords
    all_species = sub.species + over.species
    all_coords = np.vstack([sub.cart_coords, new_over_coords])
    
    hetero = Structure(
        new_lattice,
        all_species,
        all_coords,
        coords_are_cartesian=True
    )
    
    return {
        "success": True,
        "heterostructure": structure_to_dict(hetero),
        "warnings": warnings
    }


def _calculate_surface_normal(structure: Structure, site_index: int, neighbor_cutoff: float = 5.0) -> np.ndarray:
    """
    Calculate local surface normal at a given site by fitting a plane to nearby atoms.

    Args:
        structure: Pymatgen Structure (typically a slab)
        site_index: Index of the site to calculate normal at
        neighbor_cutoff: Maximum distance to consider neighbors (Angstroms)

    Returns:
        Unit normal vector as numpy array [nx, ny, nz]
    """
    target_coords = structure[site_index].coords
    all_coords = structure.cart_coords

    # Find neighbors within cutoff
    distances = np.linalg.norm(all_coords - target_coords, axis=1)
    neighbor_mask = (distances < neighbor_cutoff) & (distances > 0.1)  # Exclude self
    neighbor_coords = all_coords[neighbor_mask]

    if len(neighbor_coords) < 3:
        # Not enough neighbors for plane fitting, fall back to z-axis
        # This happens for isolated molecules or very sparse structures
        return np.array([0.0, 0.0, 1.0])

    # Fit plane using PCA (principal component analysis)
    # The normal is the direction with least variance (smallest eigenvalue)
    centered = neighbor_coords - np.mean(neighbor_coords, axis=0)
    covariance = centered.T @ centered
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    # Normal is eigenvector with smallest eigenvalue (least variance direction)
    normal = eigenvectors[:, 0]

    # Ensure normal points "outward" (positive z-component for typical slabs)
    # This is a heuristic - for slab, we want normal pointing away from bulk
    if normal[2] < 0:
        normal = -normal

    # Normalize to unit vector
    normal = normal / np.linalg.norm(normal)

    return normal


def _rotate_around_axis(coords: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """
    Rotate coordinates around an arbitrary axis using Rodrigues' rotation formula.

    Args:
        coords: Nx3 array of coordinates to rotate
        axis: 3D unit vector defining rotation axis
        angle: Rotation angle in radians

    Returns:
        Rotated coordinates as Nx3 array
    """
    # Rodrigues' rotation formula:
    # v_rot = v*cos(θ) + (k × v)*sin(θ) + k*(k·v)*(1-cos(θ))
    # where k is the unit axis vector

    axis = axis / np.linalg.norm(axis)  # Ensure normalized
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)

    # Rotation matrix using Rodrigues formula
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    R = np.eye(3) + sin_angle * K + (1 - cos_angle) * (K @ K)

    return coords @ R.T


def add_adsorbate(
    structure_dict: Dict[str, Any],
    molecule: Any,  # Molecule dict or name string
    site_index: int,
    distance: float = 2.0
) -> Dict[str, Any]:
    """
    Add an adsorbate molecule to a specific site using surface normal orientation.

    This function calculates the local surface normal at the adsorption site by
    fitting a plane to nearby atoms, then places the molecule along this normal
    direction. This ensures correct placement even for non-(001) surfaces.

    Args:
        structure_dict: Surface structure
        molecule: Molecule structure dict OR name string (e.g., "H2O", "CO", "O2")
        site_index: Index of atom to adsorb on (top site mode)
        distance: Height above site along surface normal (Angstroms)

    Returns:
        Dictionary with success status, structure, and optional warnings
    """
    surface = dict_to_structure(structure_dict)

    if surface is None:
        return {"success": False, "error": {"code": "INVALID_INPUT", "message": "Invalid surface structure"}}

    # Validate site
    if site_index < 0 or site_index >= len(surface):
        return {"success": False, "error": {"code": "INVALID_INPUT", "message": f"Invalid site index {site_index}"}}

    # Handle molecule - support both string names and dicts
    mol_struct = None

    if isinstance(molecule, str):
        # Common molecules database for adsorbates
        common_molecules = {
            "H2O": {"species": ["O", "H", "H"], "coords": [[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]]},
            "CO": {"species": ["C", "O"], "coords": [[0, 0, 0], [1.13, 0, 0]]},
            "CO2": {"species": ["C", "O", "O"], "coords": [[0, 0, 0], [-1.16, 0, 0], [1.16, 0, 0]]},
            "O2": {"species": ["O", "O"], "coords": [[0, 0, 0], [1.21, 0, 0]]},
            "H2": {"species": ["H", "H"], "coords": [[0, 0, 0], [0.74, 0, 0]]},
            "N2": {"species": ["N", "N"], "coords": [[0, 0, 0], [1.10, 0, 0]]},
            "NH3": {"species": ["N", "H", "H", "H"], "coords": [[0, 0, 0.12], [0, 0.94, -0.27], [0.81, -0.47, -0.27], [-0.81, -0.47, -0.27]]},
            "CH4": {"species": ["C", "H", "H", "H", "H"], "coords": [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [0.63, -0.63, -0.63], [-0.63, 0.63, -0.63]]},
            "OH": {"species": ["O", "H"], "coords": [[0, 0, 0], [0.97, 0, 0]]},
            "NO": {"species": ["N", "O"], "coords": [[0, 0, 0], [1.15, 0, 0]]},
            "NO2": {"species": ["N", "O", "O"], "coords": [[0, 0, 0], [1.20, 0, 0], [-0.60, 1.04, 0]]},
            "CH3": {"species": ["C", "H", "H", "H"], "coords": [[0, 0, 0], [1.09, 0, 0], [-0.55, 0.95, 0], [-0.55, -0.95, 0]]},
        }

        mol_key = molecule.upper()

        # First check our common molecules database
        if mol_key in common_molecules:
            mol_data = common_molecules[mol_key]
            from pymatgen.core import Molecule
            mol_struct = Molecule(mol_data["species"], mol_data["coords"])
        else:
            # Check if it's a valid ASE molecule name before calling ASE
            from ase.data import chemical_symbols
            from ase.collections import g2

            # ASE g2 database contains molecule names
            ase_valid_names = set(g2.names)

            if molecule in ase_valid_names:
                from ase.build import molecule as ase_molecule
                from pymatgen.io.ase import AseAtomsAdaptor
                ase_mol = ase_molecule(molecule)
                mol_struct = AseAtomsAdaptor.get_structure(ase_mol)
            else:
                supported = ", ".join(sorted(common_molecules.keys()))
                return {"success": False, "error": {"code": "INVALID_INPUT", "message": f"Unknown molecule: {molecule}. Supported: {supported}. Also supports ASE g2 database molecules."}}

    elif isinstance(molecule, dict):
        if "structure" in molecule:
            mol_struct = dict_to_structure(molecule["structure"])  # If wrapper dict
        else:
            mol_struct = dict_to_structure(molecule)  # If direct dict
    else:
        return {"success": False, "error": {"code": "INVALID_INPUT", "message": "Molecule must be a name string or structure dict"}}

    if mol_struct is None:
        return {"success": False, "error": {"code": "INVALID_INPUT", "message": "Failed to parse molecule structure"}}
    
    target_site = surface[site_index]
    target_coords = target_site.coords

    # Calculate local surface normal at the adsorption site
    # Find nearby atoms and fit a plane to determine normal direction
    surface_normal = _calculate_surface_normal(surface, site_index)

    # Normalize molecule coordinates so lowest point is at origin
    mol_coords = np.array(mol_struct.cart_coords, dtype=float)
    # Project molecule onto surface normal direction to find extent
    mol_proj = mol_coords @ surface_normal
    min_proj = np.min(mol_proj)
    # Shift molecule so its lowest point (along normal) is at zero
    mol_coords = mol_coords - min_proj * surface_normal

    surface_coords = surface.cart_coords
    angles = [0, 90, 180, 270] if len(mol_coords) > 1 else [0]
    best_coords = None
    best_min_dist = -1.0
    best_angle = 0

    for angle in angles:
        # Rotate around the surface normal (not just z-axis)
        theta = np.deg2rad(angle)
        rotated = _rotate_around_axis(mol_coords, surface_normal, theta)

        # Center molecule COM on target site in the plane perpendicular to normal
        mol_com = np.mean(rotated, axis=0)

        # Decompose: component along normal and perpendicular component
        target_perp = target_coords - (target_coords @ surface_normal) * surface_normal
        mol_com_perp = mol_com - (mol_com @ surface_normal) * surface_normal

        # Place molecule: align perpendicular components, then offset along normal
        shift = target_perp - mol_com_perp + (target_coords @ surface_normal + distance) * surface_normal
        placed = rotated + shift

        if len(surface_coords) > 0:
            deltas = placed[:, None, :] - surface_coords[None, :, :]
            min_dist = float(np.min(np.linalg.norm(deltas, axis=2)))
        else:
            min_dist = float("inf")

        if min_dist > best_min_dist:
            best_min_dist = min_dist
            best_coords = placed
            best_angle = angle

    warnings = []
    if best_min_dist < 1.0:
        warnings.append(f"Adsorbate is very close to surface (min distance {best_min_dist:.2f} Angstroms)")
    if len(angles) > 1:
        warnings.append(f"Selected rotation {best_angle} degrees around surface normal")

    new_surface = surface.copy()
    for i, spec in enumerate(mol_struct.species):
        new_surface.append(spec, best_coords[i], coords_are_cartesian=True)

    result = {
        "success": True,
        "structure": structure_to_dict(new_surface)
    }
    if warnings:
        result["warnings"] = warnings

    return result


def main():
    """
    Main entry point for command-line execution.
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python structure_tools.py <input_json_file>"
            }
        }))
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input file
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    # Determine operation
    operation = params.pop("operation", None)  # Remove operation key before passing to functions

    if operation == "supercell":
        result = make_supercell(**params)
    elif operation == "slab":
        result = generate_slab(**params)
    elif operation == "defect":
        result = create_defect(**params)
    elif operation == "strain":
        result = apply_strain(**params)
    elif operation == "alloy":
        result = create_alloy(**params)
    elif operation == "heterostructure":
        result = create_heterostructure(**params)
    elif operation == "adsorbate":
        result = add_adsorbate(**params)
    else:
        result = {
            "success": False,
            "error": {
                "code": "INVALID_OPERATION",
                "message": f"Unknown operation: {operation}",
                "details": {"operation": operation}
            }
        }
    
    # Output result
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
