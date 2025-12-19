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


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def make_supercell(
    structure: Structure,
    scaling: Union[List[int], List[List[int]], str]
) -> Dict[str, Any]:
    """
    Create supercell from structure.
    
    Args:
        structure: Input structure
        scaling: [nx, ny, nz] or transformation matrix or preset name
    
    Returns:
        Supercell structure
    """
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
    
    Args:
        base_structure: Bulk structure
        miller_index: Surface orientation
        slab_thickness: Number of layers
        vacuum: Vacuum thickness in Angstrom
    
    Returns:
        Slab supercell
    """
    a = base_structure.lattice.a
    
    # Simplified slab generation
    h, k, l = miller_index
    
    # Surface lattice vectors
    if miller_index == [1, 1, 1]:
        a_surf = a / np.sqrt(2)
        interlayer = a / np.sqrt(3)
    elif miller_index == [1, 1, 0]:
        a_surf = a
        interlayer = a / np.sqrt(8)
    else:
        a_surf = a
        interlayer = a / np.sqrt(h**2 + k**2 + l**2)
    
    c_total = slab_thickness * interlayer + vacuum
    
    lattice = Lattice.orthorhombic(a_surf * 2, a_surf * 2 * np.sqrt(3), c_total)
    
    species = []
    coords = []
    
    base_elem = str(base_structure[0].specie)
    
    for layer in range(slab_thickness):
        z = (layer + 0.5) * interlayer / c_total
        for i in range(4):
            for j in range(4):
                x = (i + 0.5 * (layer % 2)) / 4
                y = (j + (layer % 3) * 0.33) / 4
                species.append(base_elem)
                coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "miller_index": miller_index,
        "slab_thickness_layers": slab_thickness,
        "vacuum_A": vacuum,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
