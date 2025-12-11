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
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator


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
    
    if "lattice" not in structure_dict or "atoms" not in structure_dict:
        return None
    
    lattice_data = structure_dict["lattice"]
    if "matrix" not in lattice_data:
        return None
    
    lattice_matrix = lattice_data["matrix"]
    lattice = Lattice(lattice_matrix)
    
    atoms_data = structure_dict["atoms"]
    if not atoms_data:
        return None
    
    species = []
    coords = []
    
    for atom in atoms_data:
        if "element" not in atom or "coords" not in atom:
            return None
        species.append(atom["element"])
        coords.append(atom["coords"])
    
    structure = Structure(lattice, species, coords, coords_are_cartesian=False)
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
    
    result = {
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
        "metadata": {
            "natoms": len(structure),
            "formula": structure.composition.formula
        }
    }
    
    for site in structure:
        result["atoms"].append({
            "element": str(site.specie),
            "coords": site.frac_coords.tolist(),
            "cartesian": site.coords.tolist()
        })
    
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
                    "message": f"Minimum distance {min_dist:.3f} Å is too small",
                    "details": {"min_distance": min_dist},
                    "suggestion": "Check scaling matrix or original structure"
                }
            }
    
    # Convert to dictionary
    supercell_dict = structure_to_dict(supercell)
    
    # Note: Space group detection would require spglib integration
    # For now, we assume symmetry is lost
    original_spg = structure_dict.get("space_group", {}).get("number", None)
    
    return {
        "success": True,
        "supercell": supercell_dict,
        "transformation_matrix": scaling_array.tolist(),
        "original_space_group": original_spg,
        "new_space_group": 1,  # Typically P1 after supercell
        "volume_multiplier": volume_multiplier,
        "preserved_symmetry": False
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
                "message": f"Thickness must be >= 1, got {thickness}",
                "details": {"thickness": thickness}
            }
        }
    
    if vacuum < 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"Vacuum must be >= 0, got {vacuum}",
                "details": {"vacuum": vacuum}
            }
        }
    
    # Generate slab using Pymatgen
    slabgen = SlabGenerator(
        structure,
        miller_indices,
        min_slab_size=max(min_slab_size, thickness),
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
            "original_space_group": structure_dict.get("space_group", {}).get("number", None),
            "slab_space_group": 1,  # Typically P1
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
        # Add interstitial at nearby position (simplified)
        site = structure[defect_site]
        new_coords = site.frac_coords + np.array([0.1, 0.1, 0.1])
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
    strain_tensor: List[List[float]],
    strain_type: Optional[str] = None,
    strain_value: Optional[float] = None
) -> Dict[str, Any]:
    """
    Apply strain to structure.
    
    Args:
        structure_dict: Original structure
        strain_tensor: 3x3 strain tensor
        strain_type: Predefined strain type
        strain_value: Strain magnitude
    
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
    
    # Apply strain
    strained = structure.copy()
    deformation = np.eye(3) + np.array(strain_tensor)
    strained.apply_strain(deformation)
    
    # Convert to dictionary
    strained_dict = structure_to_dict(strained)
    
    return {
        "success": True,
        "strained_structure": strained_dict,
        "strain_tensor": strain_tensor,
        "deformation_matrix": deformation.tolist()
    }


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
        
    # Stack
    # We essentially put 'over' on top of 'sub'
    # 1. Expand c lattice to hold both + gap + vacuum
    # 2. Shift overlayer
    
    # Get z-range of substrate
    sub_coords = sub.cart_coords
    max_z_sub = np.max(sub_coords[:, 2]) if len(sub) > 0 else 0
    min_z_sub = np.min(sub_coords[:, 2]) if len(sub) > 0 else 0
    
    # Get z-range of overlayer
    over_coords = over.cart_coords
    max_z_over = np.max(over_coords[:, 2]) if len(over) > 0 else 0
    min_z_over = np.min(over_coords[:, 2]) if len(over) > 0 else 0
    over_height = max_z_over - min_z_over
    
    # New z position for overlayer bottom
    # Place overlayer bottom at max_z_sub + interface_distance
    # Shift overlayer: Subtract its min_z first, then add offset
    shift_z = (max_z_sub + interface_distance) - min_z_over
    
    new_over_coords = over_coords + np.array([0, 0, shift_z])
    
    # Combine atoms
    # We use substrate's lattice a, b. And new c.
    new_c = (max_z_sub + interface_distance + over_height + vacuum)
    
    # Ensure c is orthogonal to a,b?
    # Pymatgen Structure expects Lattice.
    # We take sub.lattice but modify c
    # Assuming sub is oriented with c along Z for simplicity?
    # If not, this is complex.
    # We'll assume orthogonal-ish or just set c length?
    # Safer: Create new Lattice with sub's a,b, gamma and new c aligned with Z
    
    # For now, simplistic approach:
    new_lattice = Lattice.from_parameters(
        a1, b1, new_c, 
        sub.lattice.angles[0], sub.lattice.angles[1], sub.lattice.angles[2]
    )
    
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


def add_adsorbate(
    structure_dict: Dict[str, Any],
    molecule: Dict[str, Any], # Molecule dict or name?
    site_index: int,
    distance: float = 2.0
) -> Dict[str, Any]:
    """
    Add an adsorbate molecule to a specific site.
    
    Args:
        structure_dict: Surface structure
        molecule: Molecule structure dict OR name string?
                  Ideally a dict from molecule_generator.
        site_index: Index of atom to adsorb on (top site mode)
        distance: Height above site
    """
    surface = dict_to_structure(structure_dict)
    if "structure" in molecule:
        mol_struct = dict_to_structure(molecule["structure"]) # If wrapper dict
    else:
        mol_struct = dict_to_structure(molecule) # If direct dict
        
    if surface is None:
         return {"success": False, "error": "Invalid surface structure"}
         
    # Validate site
    if site_index < 0 or site_index >= len(surface):
         return {"success": False, "error": f"Invalid site index {site_index}"}
         
    # Handle molecule
    if mol_struct is None:
         # Maybe it's just a species string?
         return {"success": False, "error": "Invalid molecule structure"}
    
    # Strategy:
    # 1. Identify "anchoring atom" in molecule? (First atom?)
    # 2. Translate molecule so anchor is at site_coords + [0, 0, distance]
    #    (Assuming Z is normal)
    
    target_site = surface[site_index]
    target_coords = target_site.coords
    
    # Molecule center/anchor
    # We'll validly assume the molecule is centered or user wants CoM?
    # Let's align CoM to target? Or bottom-most atom?
    # Bottom-most is safer for adsorption.
    mol_coords = mol_struct.cart_coords
    min_z_mol = np.min(mol_coords[:, 2])
    
    # Shift molecule so min_z is at 0
    mol_coords = mol_coords - np.array([0, 0, min_z_mol])
    
    # Now shift to target (+ distance)
    # We assume 'distance' is bond length approx.
    final_origin = target_coords + np.array([0, 0, distance])
    
    # Center molecule horizontally relative to target?
    # Assume molecule CoM (x,y) should match target (x,y)
    mol_com = np.mean(mol_coords, axis=0)
    shift_xy = final_origin[:2] - mol_com[:2]
    
    # Apply shift
    total_shift = np.array([shift_xy[0], shift_xy[1], final_origin[2]])
    # Note: final_origin[2] is target_z + distance. 
    # Current mol min_z is 0. So adding final_origin[2] puts bottom at target+dist.
    
    mol_coords = mol_coords + np.array([shift_xy[0], shift_xy[1], final_origin[2]]) # Wait, shift_xy plus Z move
    
    # Wait, mol_coords was zeroed at Z.
    # We want min_z to be at target_z + distance.
    # So Z shift = (target_z + distance) - (original_min_z=0) => target_z + distance.
    # Yes.
    
    # Add atoms to surface
    new_surface = surface.copy()
    for i, spec in enumerate(mol_struct.species):
        new_surface.append(spec, mol_coords[i], coords_are_cartesian=True)
        
    return {
        "success": True,
        "structure": structure_to_dict(new_surface)
    }


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
    operation = params.get("operation", None)
    
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
