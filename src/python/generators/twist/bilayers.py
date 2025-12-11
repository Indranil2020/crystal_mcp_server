"""
twist/bilayers.py - Twisted Bilayer Structures

Generates twisted homo- and hetero-bilayers:
- Twisted bilayer graphene (TBG)
- TMD bilayers (MoS2/MoS2, WS2/WS2)
- Graphene/h-BN heterostructures
- TMD/TMD heterostructures

Scientific basis:
- Lattice matching algorithms
- Strain distribution options
- AB stacking → twisted transition
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import rotation_matrix_2d, structure_to_dict, calculate_moire_angle


# Common bilayer combinations
BILAYER_DATABASE = {
    "graphene/graphene": {"a1": 2.46, "a2": 2.46, "z_sep": 3.35, "e1": "C", "e2": "C"},
    "graphene/hBN": {"a1": 2.46, "a2": 2.50, "z_sep": 3.30, "e1": "C", "e2": "B", "e2b": "N"},
    "MoS2/MoS2": {"a1": 3.16, "a2": 3.16, "z_sep": 6.15, "e1": "Mo", "e2": "Mo"},
    "MoS2/WS2": {"a1": 3.16, "a2": 3.18, "z_sep": 6.20, "e1": "Mo", "e2": "W"},
    "WSe2/WSe2": {"a1": 3.28, "a2": 3.28, "z_sep": 6.50, "e1": "W", "e2": "W"},
    "MoSe2/WSe2": {"a1": 3.29, "a2": 3.28, "z_sep": 6.45, "e1": "Mo", "e2": "W"},
}


def generate_twisted_bilayer(
    material: str = "graphene",
    twist_angle: float = 5.0,
    stacking: str = "AA",
    supercell_size: int = 5,
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate twisted homo-bilayer structure.
    
    Args:
        material: Material type (graphene, MoS2, WS2, etc.)
        twist_angle: Twist angle in degrees
        stacking: Initial stacking (AA, AB)
        supercell_size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Twisted bilayer structure
    """
    # Material parameters
    if material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep = 3.35
        is_tmd = False
    elif material.lower() in ["mos2", "ws2", "wse2", "mose2"]:
        a = 3.16 if "mo" in material.lower() and "s2" in material.lower() else 3.28
        element = "Mo" if "mo" in material.lower() else "W"
        z_sep = 6.5
        is_tmd = True
    else:
        a = 2.5
        element = "C"
        z_sep = 3.5
        is_tmd = False
    
    # Calculate Moiré info
    moire_info = calculate_moire_angle(twist_angle, a)
    
    # Limit supercell for practicality
    supercell_size = min(supercell_size, 30)
    
    c = z_sep + vacuum
    lattice = Lattice.hexagonal(a * supercell_size, c)
    
    species = []
    coords = []
    z_center = 0.5
    dz = z_sep / c / 2
    
    # Bottom layer (unrotated)
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size
            y1 = (2/3 + j) / supercell_size
            x2 = (2/3 + i) / supercell_size
            y2 = (1/3 + j) / supercell_size
            
            if is_tmd:
                # TMD: M-X2 sandwich
                species.extend([element, "S", "S"])
                t = 1.5 / c  # TMD thickness
                coords.extend([
                    [x1, y1, z_center - dz],
                    [x1, y1, z_center - dz - t],
                    [x1, y1, z_center - dz + t]
                ])
            else:
                species.extend([element, element])
                coords.extend([[x1, y1, z_center - dz], [x2, y2, z_center - dz]])
    
    # Top layer (rotated)
    rot = rotation_matrix_2d(twist_angle)
    ab_shift = [0, 0] if stacking == "AA" else [1/3, 0]  # AB shift
    
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size + ab_shift[0] / supercell_size
            y1 = (2/3 + j) / supercell_size + ab_shift[1] / supercell_size
            x2 = (2/3 + i) / supercell_size + ab_shift[0] / supercell_size
            y2 = (1/3 + j) / supercell_size + ab_shift[1] / supercell_size
            
            # Apply rotation
            for (x, y) in [(x1, y1), (x2, y2)] if not is_tmd else [(x1, y1)]:
                pos = np.array([x - 0.5, y - 0.5])
                pos_rot = rot @ pos + 0.5
                
                if is_tmd:
                    t = 1.5 / c
                    species.extend([element, "S", "S"])
                    coords.extend([
                        [pos_rot[0] % 1, pos_rot[1] % 1, z_center + dz],
                        [pos_rot[0] % 1, pos_rot[1] % 1, z_center + dz - t],
                        [pos_rot[0] % 1, pos_rot[1] % 1, z_center + dz + t]
                    ])
                else:
                    species.append(element)
                    coords.append([pos_rot[0] % 1, pos_rot[1] % 1, z_center + dz])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "twist_angle_deg": twist_angle,
        "initial_stacking": stacking,
        "supercell_size": supercell_size,
        "n_atoms": len(structure),
        "moire_period_angstrom": moire_info["moire_period_angstrom"],
        "interlayer_distance_angstrom": z_sep,
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_twisted_heterobilayer(
    layer1: str = "graphene",
    layer2: str = "hBN",
    twist_angle: float = 0.0,
    supercell_size: int = 5,
    strain_distribution: str = "equal",
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate twisted hetero-bilayer (different materials).
    
    Args:
        layer1: Bottom layer material
        layer2: Top layer material
        twist_angle: Twist angle in degrees
        supercell_size: Supercell size
        strain_distribution: 'equal', 'layer1', 'layer2'
        vacuum: Vacuum padding
    
    Returns:
        Twisted heterobilayer structure
    """
    # Material parameters
    params = {
        "graphene": {"a": 2.46, "e": "C", "z": 0},
        "hBN": {"a": 2.50, "e": ["B", "N"], "z": 0},
        "MoS2": {"a": 3.16, "e": "Mo", "z": 1.5},
        "WS2": {"a": 3.18, "e": "W", "z": 1.5},
        "WSe2": {"a": 3.28, "e": "W", "z": 1.6},
    }
    
    p1 = params.get(layer1, {"a": 2.5, "e": "C", "z": 0})
    p2 = params.get(layer2, {"a": 2.5, "e": "C", "z": 0})
    
    # Handle lattice mismatch
    a1, a2 = p1["a"], p2["a"]
    mismatch = (a2 - a1) / a1 * 100
    
    if strain_distribution == "equal":
        a_avg = (a1 + a2) / 2
        a1_strained = a2_strained = a_avg
    elif strain_distribution == "layer1":
        a1_strained = a2
        a2_strained = a2
    else:  # layer2
        a1_strained = a1
        a2_strained = a1
    
    z_sep = 3.3 + (p1["z"] + p2["z"]) / 2
    
    # Build structure
    a = (a1_strained + a2_strained) / 2
    supercell_size = min(supercell_size, 20)
    c = z_sep + vacuum
    
    lattice = Lattice.hexagonal(a * supercell_size, c)
    
    species = []
    coords = []
    z_center = 0.5
    dz = z_sep / c / 2
    
    # Bottom layer (layer1)
    e1 = p1["e"]
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size
            y1 = (2/3 + j) / supercell_size
            x2 = (2/3 + i) / supercell_size
            y2 = (1/3 + j) / supercell_size
            
            if isinstance(e1, list):
                species.extend(e1)
                coords.extend([[x1, y1, z_center - dz], [x2, y2, z_center - dz]])
            else:
                species.extend([e1, e1])
                coords.extend([[x1, y1, z_center - dz], [x2, y2, z_center - dz]])
    
    # Top layer (layer2, rotated)
    rot = rotation_matrix_2d(twist_angle)
    e2 = p2["e"]
    
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size
            y1 = (2/3 + j) / supercell_size
            x2 = (2/3 + i) / supercell_size
            y2 = (1/3 + j) / supercell_size
            
            for (x, y), elem in [((x1, y1), e2[0] if isinstance(e2, list) else e2), 
                                  ((x2, y2), e2[1] if isinstance(e2, list) else e2)]:
                pos = np.array([x - 0.5, y - 0.5])
                pos_rot = rot @ pos + 0.5
                species.append(elem)
                coords.append([pos_rot[0] % 1, pos_rot[1] % 1, z_center + dz])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "layer1": layer1,
        "layer2": layer2,
        "twist_angle_deg": twist_angle,
        "lattice_mismatch_percent": round(mismatch, 2),
        "strain_distribution": strain_distribution,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }
