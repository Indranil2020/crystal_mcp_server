"""
two_d/grain_boundaries.py - 2D Grain Boundaries

Generates grain boundaries in 2D materials:
- Coincident Site Lattice (CSL) boundaries
- Σ3, Σ5, Σ7 and higher
- Dislocation cores
- Tilt boundaries

Scientific basis:
- Misorientation angle and CSL theory2D GB formation during CVD growth
- Electronic properties of GBs in graphene
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Common 2D grain boundary parameters
GB_TYPES_2D = {
    "sigma3": {
        "sigma": 3,
        "angle": 60.0,
        "description": "60° tilt boundary (armchair-armchair)",
        "periodicity": "short"
    },
    "sigma5": {
        "sigma": 5,
        "angle": 36.87,
        "description": "36.87° tilt boundary",
        "periodicity": "medium"
    },
    "sigma7": {
        "sigma": 7,
        "angle": 21.79,
        "description": "21.79° tilt boundary",
        "periodicity": "long"
    },
    "sigma13": {
        "sigma": 13,
        "angle": 32.2,
        "description": "32.2° tilt boundary",
        "periodicity": "long"
    },
    "5-7": {
        "sigma": None,
        "angle": 0,
        "description": "5-7 ring defect (dislocation core)",
        "periodicity": "short"
    },
    "5-8-5": {
        "sigma": None,
        "angle": 0,
        "description": "Double pentagon-octagon pair",
        "periodicity": "short"
    }
}


def generate_2d_grain_boundary(
    base_material: str = "graphene",
    gb_type: str = "sigma5",
    length: int = 5,
    width: int = 5,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate 2D grain boundary structure.
    
    Args:
        base_material: Base 2D material (graphene, MoS2, etc.)
        gb_type: Type of grain boundary (sigma3, sigma5, sigma7, 5-7, 5-8-5)
        length: Length along GB in unit cells
        width: Width perpendicular to GB in unit cells
        vacuum: Vacuum padding
    
    Returns:
        Structure with grain boundary
    """
    if gb_type not in GB_TYPES_2D:
        return {
            "success": False,
            "error": {
                "code": "INVALID_GB_TYPE",
                "message": f"Unknown GB type '{gb_type}'",
                "available": list(GB_TYPES_2D.keys())
            }
        }
    
    gb_info = GB_TYPES_2D[gb_type]
    
    if base_material.lower() == "graphene":
        a = 2.46
        element = "C"
    elif base_material.lower() in ["mos2", "ws2"]:
        a = 3.16
        element = "Mo" if "mo" in base_material.lower() else "W"
    else:
        a = 2.5
        element = "C"
    
    # Build bicrystal with GB
    if gb_type in ["5-7", "5-8-5"]:
        # Local defect
        return _generate_dislocation_core(base_material, gb_type, length, width, vacuum)
    
    theta = np.radians(gb_info["angle"])
    
    # Create two grains with different orientations
    # GB along y-axis
    total_width = width * 2 * a
    total_length = length * a * np.sqrt(3)
    c = vacuum
    
    lattice = Lattice.orthorhombic(total_width, total_length, c)
    
    species = []
    coords = []
    z = 0.5
    
    # Left grain (unrotated)
    for i in range(width):
        for j in range(length * 2):
            x = (i * a + a/3) / total_width
            y = (j * a * np.sqrt(3) / 4) / total_length
            
            if x < 0.5:  # Left half
                species.append(element)
                coords.append([x, y, z])
    
    # Right grain (rotated by theta)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    for i in range(width):
        for j in range(length * 2):
            x0 = i * a + a/3
            y0 = j * a * np.sqrt(3) / 4
            
            # Rotate around origin
            x_rot = x0 * cos_t - y0 * sin_t
            y_rot = x0 * sin_t + y0 * cos_t
            
            x = (x_rot + total_width/2) / total_width
            y = (y_rot % total_length) / total_length
            
            if 0.5 < x < 1 and 0 < y < 1:
                species.append(element)
                coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "base_material": base_material,
        "gb_type": gb_type,
        "sigma": gb_info["sigma"],
        "misorientation_angle": gb_info["angle"],
        "description": gb_info["description"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }


def _generate_dislocation_core(
    base_material: str,
    core_type: str,
    length: int,
    width: int,
    vacuum: float
) -> Dict[str, Any]:
    """
    Generate dislocation core (5-7 or 5-8-5 ring defect).
    """
    if base_material.lower() == "graphene":
        a = 2.46
        element = "C"
    else:
        a = 3.0
        element = "C"
    
    # Create graphene supercell
    total_width = width * a * 2
    total_length = length * a * np.sqrt(3)
    c = vacuum
    
    lattice = Lattice.orthorhombic(total_width, total_length, c)
    
    species = []
    coords = []
    z = 0.5
    
    center_x = 0.5
    center_y = 0.5
    
    for i in range(width * 2):
        for j in range(length * 2):
            x = (1/3 + i) / (width * 2)
            y = (2/3 + j) / (length * 2)
            
            # Check if near defect center
            dist = ((x - center_x)**2 + (y - center_y)**2)**0.5
            
            if dist > 0.15:  # Outside defect region
                species.append(element)
                coords.append([x, y, z])
    
    # Add defect atoms
    if core_type == "5-7":
        # Pentagon-heptagon pair
        defect_coords = [
            [center_x - 0.05, center_y, z],
            [center_x + 0.05, center_y, z],
            [center_x, center_y - 0.08, z],
            [center_x, center_y + 0.08, z],
        ]
    else:  # 5-8-5
        defect_coords = [
            [center_x - 0.08, center_y, z],
            [center_x + 0.08, center_y, z],
            [center_x, center_y - 0.05, z],
            [center_x, center_y + 0.05, z],
            [center_x - 0.04, center_y - 0.04, z],
            [center_x + 0.04, center_y + 0.04, z],
        ]
    
    for dc in defect_coords:
        species.append(element)
        coords.append(dc)
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "base_material": base_material,
        "gb_type": core_type,
        "description": GB_TYPES_2D[core_type]["description"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_polycrystalline_2d(
    base_material: str = "graphene",
    n_grains: int = 4,
    box_size: float = 50.0,  # Angstroms
    vacuum: float = 15.0,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate polycrystalline 2D material with multiple grains.
    
    Args:
        base_material: Base 2D material
        n_grains: Number of grains
        box_size: Simulation box size
        vacuum: Vacuum padding
        seed: Random seed
    
    Returns:
        Polycrystalline structure
    """
    if seed is not None:
        np.random.seed(seed)
    
    if base_material.lower() == "graphene":
        a = 2.46
        element = "C"
    else:
        a = 3.0
        element = "C"
    
    # Generate random grain centers and orientations
    grain_centers = np.random.rand(n_grains, 2) * box_size
    grain_angles = np.random.rand(n_grains) * 60  # 0-60 degrees
    
    lattice = Lattice.orthorhombic(box_size, box_size, vacuum)
    
    species = []
    coords = []
    z = 0.5
    
    # For each lattice point, assign to nearest grain
    nx = int(box_size / a)
    ny = int(box_size / (a * np.sqrt(3) / 2))
    
    for i in range(nx):
        for j in range(ny):
            x = i * a
            y = j * a * np.sqrt(3) / 2
            
            # Find nearest grain
            dists = np.linalg.norm(grain_centers - np.array([x, y]), axis=1)
            grain_idx = np.argmin(dists)
            theta = np.radians(grain_angles[grain_idx])
            
            # Rotate position by grain angle
            cos_t = np.cos(theta)
            sin_t = np.sin(theta)
            
            x_rot = x * cos_t - y * sin_t
            y_rot = x * sin_t + y * cos_t
            
            # Wrap to box
            if 0 <= x_rot < box_size and 0 <= y_rot < box_size:
                species.append(element)
                coords.append([x_rot / box_size, y_rot / box_size, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "base_material": base_material,
        "n_grains": n_grains,
        "grain_angles": [round(a, 1) for a in grain_angles],
        "box_size_angstrom": box_size,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }
