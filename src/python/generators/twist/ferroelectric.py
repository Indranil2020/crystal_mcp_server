"""
twist/ferroelectric.py - 2D Ferroelectric Materials

Generates layer-locked ferroelectric 2D materials:
- CuInP2S6 (CIPS)
- In2Se3
- MoTe2 (distorted)
- SnTe
- Other polar 2D materials

Scientific basis:
- Out-of-plane spontaneous polarization
- Layer-dependent switching
- Sliding ferroelectricity
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# 2D Ferroelectric materials database
FERROELECTRIC_2D_MATERIALS = {
    "CuInP2S6": {
        "formula": "CuInP2S6",
        "a": 6.10,
        "thickness": 6.5,
        "polarization_direction": "out-of-plane",
        "tc_kelvin": 315,
        "polarization_uC_cm2": 4.0,
        "description": "Most studied 2D layered ferroelectric"
    },
    "In2Se3": {
        "formula": "In2Se3",
        "a": 4.00,
        "thickness": 9.5,
        "polarization_direction": "out-of-plane",
        "tc_kelvin": 450,
        "polarization_uC_cm2": 1.5,
        "description": "Room-temperature 2D ferroelectric"
    },
    "SnTe": {
        "formula": "SnTe",
        "a": 4.50,
        "thickness": 3.0,
        "polarization_direction": "in-plane",
        "tc_kelvin": 270,
        "polarization_uC_cm2": 8.0,
        "description": "In-plane ferroelectric"
    },
    "1T'-MoTe2": {
        "formula": "MoTe2",
        "a": 3.52,
        "b": 6.33,
        "thickness": 3.6,
        "polarization_direction": "in-plane",
        "tc_kelvin": 200,
        "polarization_uC_cm2": 0.5,
        "description": "Distorted TMD ferroelectric"
    },
    "WTe2": {
        "formula": "WTe2",
        "a": 3.48,
        "b": 6.25,
        "thickness": 7.0,  # Bilayer
        "polarization_direction": "out-of-plane",
        "tc_kelvin": 350,
        "polarization_uC_cm2": 0.2,
        "description": "Sliding ferroelectricity in bilayer"
    },
    "BN_bilayer": {
        "formula": "BN",
        "a": 2.50,
        "thickness": 6.6,
        "polarization_direction": "out-of-plane",
        "tc_kelvin": None,  # Stacking-dependent
        "polarization_uC_cm2": 0.1,
        "description": "Moiré ferroelectricity in twisted h-BN"
    }
}


def generate_ferroelectric_2d(
    material: str,
    polarization_state: str = "up",
    n_layers: int = 1,
    size: List[int] = [1, 1, 1],
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate 2D ferroelectric structure.
    
    Args:
        material: Material name (CuInP2S6, In2Se3, etc.)
        polarization_state: 'up', 'down', or 'paraelectric'
        n_layers: Number of layers
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Ferroelectric 2D structure
    """
    if material not in FERROELECTRIC_2D_MATERIALS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MATERIAL",
                "message": f"Unknown ferroelectric '{material}'",
                "available": list(FERROELECTRIC_2D_MATERIALS.keys())
            }
        }
    
    params = FERROELECTRIC_2D_MATERIALS[material]
    a = params["a"]
    b = params.get("b", a)
    thickness = params["thickness"]
    
    c = thickness * n_layers + vacuum
    
    if abs(a - b) < 0.01:
        lattice = Lattice.hexagonal(a, c)
    else:
        lattice = Lattice.orthorhombic(a, b, c)
    
    z_center = 0.5
    dz = thickness / c / 2
    
    # Build structure based on material
    if material == "CuInP2S6":
        species = ["Cu", "In", "P", "P", "S", "S", "S", "S", "S", "S"]
        
        # Polarization displacement
        pol_shift = 0.05 if polarization_state == "up" else (-0.05 if polarization_state == "down" else 0)
        
        coords = [
            [0.0, 0.0, z_center + pol_shift],  # Cu (mobile)
            [0.33, 0.33, z_center],             # In
            [0.17, 0.5, z_center + 0.08],       # P dimer
            [0.5, 0.17, z_center - 0.08],       # P dimer
            [0.1, 0.1, z_center + dz], [0.9, 0.1, z_center + dz],
            [0.1, 0.9, z_center - dz], [0.5, 0.5, z_center + dz],
            [0.3, 0.7, z_center - dz], [0.7, 0.3, z_center - dz],
        ]
    
    elif material == "In2Se3":
        species = ["In", "In", "Se", "Se", "Se"]
        
        pol_shift = 0.03 if polarization_state == "up" else (-0.03 if polarization_state == "down" else 0)
        
        coords = [
            [0.33, 0.67, z_center + pol_shift],
            [0.67, 0.33, z_center - pol_shift],
            [0.0, 0.0, z_center + dz],
            [0.33, 0.67, z_center],
            [0.67, 0.33, z_center - dz],
        ]
    
    elif material == "SnTe":
        species = ["Sn", "Te"]
        
        # In-plane polarization
        pol_shift = 0.02 if polarization_state == "up" else (-0.02 if polarization_state == "down" else 0)
        
        coords = [
            [0.0 + pol_shift, 0.0, z_center],
            [0.5, 0.5, z_center],
        ]
    
    elif "MoTe2" in material or "WTe2" in material:
        species = ["Mo" if "Mo" in material else "W", "Te", "Te"]
        
        pol_shift = 0.02 if polarization_state == "up" else (-0.02 if polarization_state == "down" else 0)
        
        coords = [
            [0.0 + pol_shift, 0.25, z_center],
            [0.33, 0.0, z_center + dz],
            [0.67, 0.5, z_center - dz],
        ]
    
    else:
        # Generic
        species = ["B", "N"]
        coords = [[0.33, 0.67, z_center + dz], [0.67, 0.33, z_center - dz]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": material,
        "formula": params["formula"],
        "polarization_state": polarization_state,
        "polarization_direction": params["polarization_direction"],
        "polarization_uC_cm2": params["polarization_uC_cm2"],
        "tc_kelvin": params["tc_kelvin"],
        "n_layers": n_layers,
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_sliding_ferroelectric(
    material: str = "hBN",
    stacking: str = "AB",
    twist_angle: float = 0.0,
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate sliding ferroelectric bilayer.
    
    Args:
        material: Base material (hBN, graphene)
        stacking: Initial stacking (AB, BA, AA)
        twist_angle: Twist angle for Moiré ferroelectricity
        vacuum: Vacuum padding
    
    Returns:
        Sliding ferroelectric structure
    """
    if material.lower() == "hbn":
        a = 2.50
        elements = ["B", "N"]
        z_sep = 3.33
    else:
        a = 2.46
        elements = ["C", "C"]
        z_sep = 3.35
    
    c = z_sep + vacuum
    lattice = Lattice.hexagonal(a, c)
    
    z_center = 0.5
    dz = z_sep / c / 2
    
    # Stacking-dependent positions
    if stacking == "AB":
        shift = [1/3, 0]
    elif stacking == "BA":
        shift = [-1/3, 0]
    else:  # AA
        shift = [0, 0]
    
    species = elements + elements
    coords = [
        [1/3, 2/3, z_center - dz],
        [2/3, 1/3, z_center - dz],
        [(1/3 + shift[0]) % 1, (2/3 + shift[1]) % 1, z_center + dz],
        [(2/3 + shift[0]) % 1, (1/3 + shift[1]) % 1, z_center + dz],
    ]
    
    structure = Structure(lattice, species, coords)
    
    # Determine polarization from stacking
    if stacking == "AB":
        polarization = "up"
    elif stacking == "BA":
        polarization = "down"
    else:
        polarization = "none"
    
    return {
        "success": True,
        "material": material,
        "stacking": stacking,
        "twist_angle": twist_angle,
        "polarization_direction": polarization,
        "ferroelectric_type": "sliding",
        "structure": structure_to_dict(structure, vacuum)
    }
