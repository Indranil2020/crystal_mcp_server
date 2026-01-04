"""
bulk/supercells.py - Supercell Generation

Comprehensive supercell generation utilities:
- Diagonal and non-diagonal supercells
- Random alloy supercells
- Special supercells for DFT
"""

from typing import Dict, Any, List, Optional, Union
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from .base import structure_to_dict


# Common supercell sizes for DFT
COMMON_SUPERCELLS = {
    "2x2x2": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    "3x3x3": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
    "4x4x4": [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
    "2x2x1": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
    "3x3x1": [[3, 0, 0], [0, 3, 0], [0, 0, 1]],
    "4x4x1": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    "sqrt2_sqrt2_1": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],  # 45° rotation
    "sqrt3_sqrt3_1": [[2, 1, 0], [-1, 1, 0], [0, 0, 1]],  # 30° rotation
}


# Special cells for calculations
SPECIAL_CELLS = {
    "Monkhorst_Pack_convergence": {
        "description": "Typical sizes for MP grid convergence",
        "sizes": [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]]
    },
    "phonon_supercells": {
        "description": "Common sizes for phonon calculations",
        "cubic": [2, 3, 4],
        "hexagonal": [[3, 3, 2], [4, 4, 2], [5, 5, 3]]
    },
    "defect_supercells": {
        "description": "Minimum sizes for defect isolation",
        "min_size_A": 10.0,
        "recommended": [[3, 3, 3], [4, 4, 4], [5, 5, 5]]
    },
}


def make_supercell(
    structure: Structure,
    scaling: Union[List[int], List[List[int]], str] = None,
    scaling_matrix: Union[List[int], List[List[int]], str] = None
) -> Dict[str, Any]:
    """
    Create supercell from structure.

    Args:
        structure: Input structure (pymatgen Structure)
        scaling: [nx, ny, nz] or transformation matrix or preset name
        scaling_matrix: Alias for scaling (for API compatibility)

    Returns:
        Supercell structure
    """
    # Accept either 'scaling' or 'scaling_matrix' parameter
    if scaling_matrix is not None and scaling is None:
        scaling = scaling_matrix

    if scaling is None:
        return {
            "success": False,
            "error": {"code": "MISSING_SCALING", "message": "Must provide 'scaling' or 'scaling_matrix'"}
        }

    # Handle preset names
    if isinstance(scaling, str):
        if scaling not in COMMON_SUPERCELLS:
            return {
                "success": False,
                "error": {"code": "INVALID_SCALING", "message": f"Unknown preset",
                          "available": list(COMMON_SUPERCELLS.keys())}
            }
        matrix = COMMON_SUPERCELLS[scaling]
    elif len(scaling) == 3 and all(isinstance(s, int) for s in scaling):
        matrix = [[scaling[0], 0, 0], [0, scaling[1], 0], [0, 0, scaling[2]]]
    else:
        matrix = scaling
    
    supercell = structure.copy()
    supercell.make_supercell(matrix)
    
    # Calculate expansion factor
    det = abs(np.linalg.det(matrix))
    
    return {
        "success": True,
        "original_n_atoms": len(structure),
        "supercell_n_atoms": len(supercell),
        "expansion_factor": int(det),
        "transformation_matrix": matrix,
        "supercell_lattice_A": [supercell.lattice.a, supercell.lattice.b, supercell.lattice.c],
        "structure": structure_to_dict(supercell)
    }


def generate_random_alloy_supercell(
    base_structure: Structure,
    substitutions: Dict[str, Dict[str, float]],
    supercell: List[int] = [3, 3, 3],
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate random alloy supercell.
    
    Args:
        base_structure: Base structure
        substitutions: {element: {new_element: fraction}}
        supercell: Supercell size
        seed: Random seed
    
    Returns:
        Random alloy supercell
    """
    np.random.seed(seed)
    
    # Make supercell first
    sc_structure = base_structure.copy()
    sc_structure.make_supercell(supercell)
    
    # Apply substitutions
    for site_idx, site in enumerate(sc_structure):
        elem = str(site.specie)
        
        if elem in substitutions:
            sub_dict = substitutions[elem]
            rand = np.random.random()
            
            cumulative = 0
            for new_elem, frac in sub_dict.items():
                cumulative += frac
                if rand < cumulative:
                    sc_structure.replace(site_idx, new_elem)
                    break
    
    # Calculate actual composition
    composition = {}
    for site in sc_structure:
        elem = str(site.specie)
        composition[elem] = composition.get(elem, 0) + 1
    
    return {
        "success": True,
        "n_atoms": len(sc_structure),
        "substitutions": substitutions,
        "actual_composition": composition,
        "random_seed": seed,
        "structure": structure_to_dict(sc_structure)
    }


def generate_graded_supercell(
    base_structure: Structure,
    element_from: str,
    element_to: str,
    gradient_axis: int = 2,
    supercell: List[int] = [2, 2, 8]
) -> Dict[str, Any]:
    """
    Generate compositionally graded supercell.
    
    Args:
        base_structure: Base structure
        element_from: Starting element
        element_to: Ending element
        gradient_axis: Axis along which to grade (0, 1, or 2)
        supercell: Supercell size
    
    Returns:
        Graded supercell
    """
    sc_structure = base_structure.copy()
    sc_structure.make_supercell(supercell)
    
    n_axis = supercell[gradient_axis]
    
    for site_idx, site in enumerate(sc_structure):
        if str(site.specie) == element_from:
            # Position along gradient axis determines probability
            pos = site.frac_coords[gradient_axis]
            prob_substitute = pos
            
            if np.random.random() < prob_substitute:
                sc_structure.replace(site_idx, element_to)
    
    return {
        "success": True,
        "n_atoms": len(sc_structure),
        "element_from": element_from,
        "element_to": element_to,
        "gradient_axis": gradient_axis,
        "structure": structure_to_dict(sc_structure)
    }


def calculate_supercell_for_size(
    base_structure: Structure,
    target_size_A: float = 15.0
) -> Dict[str, Any]:
    """
    Calculate supercell needed to reach target size.
    
    Args:
        base_structure: Base structure
        target_size_A: Target minimum dimension in Angstrom
    
    Returns:
        Recommended supercell
    """
    a = base_structure.lattice.a
    b = base_structure.lattice.b
    c = base_structure.lattice.c
    
    nx = max(1, int(np.ceil(target_size_A / a)))
    ny = max(1, int(np.ceil(target_size_A / b)))
    nz = max(1, int(np.ceil(target_size_A / c)))
    
    return {
        "success": True,
        "base_lattice_A": [a, b, c],
        "target_size_A": target_size_A,
        "recommended_supercell": [nx, ny, nz],
        "resulting_size_A": [a * nx, b * ny, c * nz],
        "n_atoms_estimate": len(base_structure) * nx * ny * nz
    }


def generate_slab_supercell(
    base_structure: Structure,
    miller_index: List[int] = [1, 1, 1],
    slab_thickness: int = 5,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate slab supercell for surface calculations.

    IMPORTANT: This function properly cuts a slab from the bulk structure,
    preserving the original composition and geometry. The slab is created
    using PyMatGen's SlabGenerator.

    Args:
        base_structure: Bulk structure (composition and geometry will be preserved)
        miller_index: Surface orientation (Miller indices [h, k, l])
        slab_thickness: Number of layers (will be converted to Angstroms via d-spacing)
        vacuum: Vacuum thickness in Angstroms

    Returns:
        Dictionary with success status and slab structure
    """
    # Validate input
    if len(miller_index) != 3:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MILLER_INDEX",
                "message": f"Miller index must have 3 components, got {len(miller_index)}",
                "details": {"miller_index": miller_index}
            }
        }

    if all(x == 0 for x in miller_index):
        return {
            "success": False,
            "error": {
                "code": "INVALID_MILLER_INDEX",
                "message": "Miller index cannot be [0, 0, 0]",
                "details": {}
            }
        }

    if slab_thickness < 1:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"Slab thickness must be >= 1 layer, got {slab_thickness}",
                "details": {"slab_thickness": slab_thickness}
            }
        }

    # Calculate d-spacing for the Miller plane using reciprocal lattice
    # This ensures "5 layers" means 5 actual atomic layers, not 5 Angstroms
    lattice = base_structure.lattice
    h, k, l = miller_index

    recip_lattice = lattice.reciprocal_lattice
    G_vector = h * recip_lattice.matrix[0] + k * recip_lattice.matrix[1] + l * recip_lattice.matrix[2]
    d_spacing = 1.0 / np.linalg.norm(G_vector) if np.linalg.norm(G_vector) > 1e-10 else 3.0

    # Convert layers to Angstroms
    min_slab_size = slab_thickness * d_spacing

    # Generate slab using PyMatGen's SlabGenerator
    # This preserves the original structure's composition and geometry
    slabgen = SlabGenerator(
        base_structure,
        miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=vacuum,
        center_slab=True,
        primitive=False,
        max_normal_search=1
    )

    slabs = slabgen.get_slabs(symmetrize=True)

    if not slabs or len(slabs) == 0:
        return {
            "success": False,
            "error": {
                "code": "SLAB_GENERATION_FAILED",
                "message": "Failed to generate slab from bulk structure",
                "details": {
                    "miller_index": miller_index,
                    "slab_thickness": slab_thickness,
                    "min_slab_size_angstrom": min_slab_size
                }
            }
        }

    # Use the first slab (typically the most symmetric)
    slab = slabs[0]

    # Calculate actual number of layers for reporting
    # (may differ slightly from requested due to periodicity)
    slab_coords_z = sorted([site.coords[2] for site in slab])
    unique_z = []
    tol = 0.1  # Angstrom tolerance for grouping layers
    for z in slab_coords_z:
        if not unique_z or abs(z - unique_z[-1]) > tol:
            unique_z.append(z)
    actual_layers = len(unique_z)

    return {
        "success": True,
        "miller_index": miller_index,
        "requested_thickness_layers": slab_thickness,
        "actual_thickness_layers": actual_layers,
        "d_spacing_angstrom": float(d_spacing),
        "slab_thickness_angstrom": float(max(slab_coords_z) - min(slab_coords_z)),
        "vacuum_angstrom": vacuum,
        "n_atoms": len(slab),
        "composition": slab.composition.formula,
        "structure": structure_to_dict(slab)
    }
