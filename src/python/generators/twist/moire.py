"""
twist/moire.py - Moiré Superlattice Generation

Generates Moiré superlattices for twisted 2D materials:
- Magic-angle graphene (TBG)
- TMD Moiré patterns
- Commensurate vs incommensurate cells

Scientific basis:
- Bistritzer-MacDonald model for flat bands
- Magic angles: 1.08°, 0.5°, 0.35°, etc.
- Commensurate cell algorithms
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import calculate_moire_angle, find_commensurate_cell, rotation_matrix_2d, structure_to_dict


# Magic angles for different materials
MAGIC_ANGLES = {
    "graphene": {
        "first": 1.08,
        "second": 0.50,
        "third": 0.35,
        "description": "Bistritzer-MacDonald flat band angles"
    },
    "MoS2": {
        "first": 2.0,
        "description": "TMD magic angles (estimated)"
    },
    "WSe2": {
        "first": 1.8,
        "description": "TMD magic angles"
    }
}


def generate_moire_superlattice(
    material: str = "graphene",
    twist_angle: float = 1.08,
    n_layers: int = 2,
    vacuum: float = 20.0,
    use_commensurate: bool = True
) -> Dict[str, Any]:
    """
    Generate Moiré superlattice structure.
    
    Args:
        material: Base 2D material (graphene, MoS2, etc.)
        twist_angle: Twist angle in degrees
        n_layers: Number of layers (2 for bilayer)
        vacuum: Vacuum padding
        use_commensurate: Use commensurate approximation
    
    Returns:
        Moiré superlattice structure
    """
    # Material parameters
    if material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep = 3.35
    elif material.lower() in ["mos2", "ws2", "wse2", "mose2"]:
        a = 3.16
        element = "Mo" if "mo" in material.lower() else "W"
        z_sep = 6.5
    else:
        a = 2.5
        element = "C"
        z_sep = 3.5
    
    # Calculate Moiré period
    moire_info = calculate_moire_angle(twist_angle, a)
    L_M = moire_info["moire_period_angstrom"]
    
    if L_M > 500:  # Too large
        return {
            "success": False,
            "error": {
                "code": "TOO_LARGE",
                "message": f"Moiré period {L_M:.1f} Å too large for atomistic model",
                "suggestion": "Use larger twist angle or continuum model"
            }
        }
    
    # Find commensurate cell
    if use_commensurate:
        comm_cell = find_commensurate_cell(twist_angle, a)
        if comm_cell:
            m, r, nx, ny = comm_cell
            # Use commensurate supercell
            actual_size = max(nx, ny)
        else:
            actual_size = max(1, int(L_M / a))
    else:
        actual_size = max(1, int(L_M / a))
    
    # Limit size for practical computation
    actual_size = min(actual_size, 30)
    
    c = z_sep * n_layers + vacuum
    lattice = Lattice.hexagonal(a * actual_size, c)
    
    species = []
    coords = []
    
    for layer_idx in range(n_layers):
        z_layer = (layer_idx * z_sep + vacuum/2) / c
        layer_angle = 0 if layer_idx == 0 else twist_angle * layer_idx
        rot = rotation_matrix_2d(layer_angle)
        
        # Generate atoms in supercell
        for i in range(actual_size):
            for j in range(actual_size):
                # Honeycomb positions
                for sub in range(2):
                    if sub == 0:
                        x0, y0 = 1/3, 2/3
                    else:
                        x0, y0 = 2/3, 1/3
                    
                    x_sc = (x0 + i) / actual_size
                    y_sc = (y0 + j) / actual_size
                    
                    # Apply rotation
                    pos_2d = np.array([x_sc - 0.5, y_sc - 0.5])
                    pos_rot = rot @ pos_2d + 0.5
                    
                    species.append(element)
                    coords.append([pos_rot[0] % 1, pos_rot[1] % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    # Check for magic angle
    is_magic = False
    if material.lower() == "graphene":
        for magic_name, magic_val in [("first", 1.08), ("second", 0.50), ("third", 0.35)]:
            if abs(twist_angle - magic_val) < 0.05:
                is_magic = True
                break
    
    return {
        "success": True,
        "material": material,
        "twist_angle_deg": twist_angle,
        "n_layers": n_layers,
        "moire_period_angstrom": L_M,
        "supercell_size": actual_size,
        "n_atoms": len(structure),
        "is_magic_angle": is_magic,
        "flat_bands_expected": is_magic,
        "structure": structure_to_dict(structure, vacuum)
    }


def calculate_moire_periodicity(
    a1: float,
    a2: float,
    theta: float
) -> Dict[str, Any]:
    """
    Calculate Moiré periodicity for hetero-bilayer.
    
    Args:
        a1: Lattice constant of layer 1
        a2: Lattice constant of layer 2
        theta: Twist angle
    
    Returns:
        Moiré parameters
    """
    delta = (a2 - a1) / a1
    theta_rad = np.radians(theta)
    
    # General formula for heterobilayer Moiré
    L_M = a1 / np.sqrt(delta**2 + 4 * np.sin(theta_rad/2)**2 * (1 + delta/2)**2)
    
    # Rotation angle of Moiré pattern
    phi = np.arctan2(np.sin(theta_rad), 1 + delta - np.cos(theta_rad))
    
    return {
        "moire_period_angstrom": round(L_M, 2),
        "moire_rotation_deg": round(np.degrees(phi), 2),
        "lattice_mismatch_percent": round(abs(delta) * 100, 3),
        "twist_angle_deg": theta
    }


def generate_magic_angle_graphene(
    magic_order: int = 1,
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate magic-angle twisted bilayer graphene (TBG).
    
    Args:
        magic_order: 1=first magic angle (1.08°), 2=second (0.5°), etc.
        vacuum: Vacuum padding
    
    Returns:
        Magic-angle TBG structure
    """
    magic_angles_list = [1.08, 0.50, 0.35, 0.24]
    
    if magic_order < 1 or magic_order > len(magic_angles_list):
        return {
            "success": False,
            "error": {
                "code": "INVALID_ORDER",
                "message": f"Magic order must be 1-{len(magic_angles_list)}",
                "available": magic_angles_list
            }
        }
    
    theta = magic_angles_list[magic_order - 1]
    
    result = generate_moire_superlattice("graphene", theta, 2, vacuum)
    
    if result["success"]:
        result["magic_order"] = magic_order
        result["expected_properties"] = {
            "flat_bands": True,
            "correlated_insulator": True,
            "superconductivity": magic_order == 1
        }
    
    return result
