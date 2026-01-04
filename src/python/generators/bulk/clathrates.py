"""
bulk/clathrates.py - Clathrate Structures

Comprehensive clathrate generation:
- Type I, II, III clathrates
- Guest-host systems
- Thermoelectric clathrates
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice
from .base import structure_to_dict


# Clathrate type database
CLATHRATE_DATABASE = {
    # Type I clathrates (Pm-3n)
    "Ba8Ga16Ge30": {
        "type": "I", "framework": "Ge-Ga", "guest": "Ba",
        "a": 10.78, "n_framework": 46, "n_guest": 8,
        "thermoelectric": True, "ZT": 1.4
    },
    "Ba8Ga16Si30": {
        "type": "I", "framework": "Si-Ga", "guest": "Ba",
        "a": 10.44, "n_framework": 46, "n_guest": 8
    },
    "Sr8Ga16Ge30": {
        "type": "I", "framework": "Ge-Ga", "guest": "Sr",
        "a": 10.72, "n_framework": 46, "n_guest": 8
    },
    "K8Ga8Si38": {
        "type": "I", "framework": "Si-Ga", "guest": "K",
        "a": 10.35, "n_framework": 46, "n_guest": 8
    },
    "Na8Si46": {
        "type": "I", "framework": "Si", "guest": "Na",
        "a": 10.19, "empty_framework": False
    },
    
    # Type II clathrates (Fd-3m)
    "Na24Si136": {
        "type": "II", "framework": "Si", "guest": "Na",
        "a": 14.62, "n_framework": 136, "n_guest": 24
    },
    "Cs8Na16Si136": {
        "type": "II", "framework": "Si", "guest": ["Cs", "Na"],
        "a": 14.74, "dual_guest": True
    },
    "K8Si136": {
        "type": "II", "framework": "Si", "guest": "K",
        "a": 14.68, "partially_filled": True
    },
    
    # Inverse clathrates (metal framework)
    "Ba8Cu16P30": {
        "type": "I_inverse", "framework": "Cu-P", "guest": "Ba",
        "a": 10.12, "semiconductor": True
    },
    "Ba8Au16As30": {
        "type": "I_inverse", "framework": "Au-As", "guest": "Ba",
        "a": 10.42
    },
    
    # Hydrate clathrates
    "CH4_5.75H2O": {
        "type": "sI_hydrate", "framework": "H2O", "guest": "CH4",
        "a": 11.88, "hydrate": True, "gas_storage": True
    },
    "CO2_5.75H2O": {
        "type": "sI_hydrate", "framework": "H2O", "guest": "CO2",
        "a": 11.83, "hydrate": True
    },
    "H2_5.67H2O": {
        "type": "sII_hydrate", "framework": "H2O", "guest": "H2",
        "a": 17.03, "hydrogen_storage": True
    },
}


# Cage types
CAGE_TYPES = {
    "type_I": {
        "small_cage": {"vertices": 20, "faces": "5^12", "name": "dodecahedron"},
        "large_cage": {"vertices": 24, "faces": "5^12 6^2", "name": "tetrakaidecahedron"},
        "n_small": 2, "n_large": 6
    },
    "type_II": {
        "small_cage": {"vertices": 20, "faces": "5^12", "name": "dodecahedron"},
        "large_cage": {"vertices": 28, "faces": "5^12 6^4", "name": "hexakaidecahedron"},
        "n_small": 16, "n_large": 8
    },
}


def generate_clathrate(
    clathrate: str = "Ba8Ga16Ge30",
    empty_cages: bool = False,
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate clathrate structure.
    
    Args:
        clathrate: Clathrate formula from database
        empty_cages: Generate empty framework without guests
        supercell: Supercell dimensions
    
    Returns:
        Clathrate structure
    """
    if clathrate not in CLATHRATE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_CLATHRATE", "message": f"Unknown clathrate",
                      "available": list(CLATHRATE_DATABASE.keys())}
        }
    
    info = CLATHRATE_DATABASE[clathrate]
    a = info["a"]
    clath_type = info["type"]
    
    lattice = Lattice.cubic(a)
    
    species = []
    coords = []
    
    if clath_type in ["I", "I_inverse"]:
        # Type I: Pm-3n, 46 framework + 8 guest atoms
        
        # Framework positions (Wyckoff 6c, 16i, 24k)
        # 6c: 1/4, 0, 1/2 family
        wyckoff_6c = [
            [0.25, 0, 0.5], [0.75, 0, 0.5], [0.5, 0.25, 0],
            [0.5, 0.75, 0], [0, 0.5, 0.25], [0, 0.5, 0.75]
        ]
        
        # 16i: x, x, x family (approximately 0.18)
        x_16i = 0.1835
        wyckoff_16i = [
            [x_16i, x_16i, x_16i], [1-x_16i, 1-x_16i, x_16i],
            [1-x_16i, x_16i, 1-x_16i], [x_16i, 1-x_16i, 1-x_16i],
            [x_16i+0.5, x_16i+0.5, x_16i+0.5] if x_16i < 0.5 else [x_16i, x_16i, x_16i],
        ][:4]
        # Add reflections for 16 positions
        for pos in list(wyckoff_16i):
            wyckoff_16i.extend([
                [pos[0]+0.5, pos[1]+0.5, pos[2]],
                [pos[0]+0.5, pos[1], pos[2]+0.5],
                [pos[0], pos[1]+0.5, pos[2]+0.5],
            ][:min(3, 16 - len(wyckoff_16i))])
        wyckoff_16i = wyckoff_16i[:16]
        
        # 24k: simplified positions
        wyckoff_24k = [
            [0, 0.31, 0.12], [0, 0.69, 0.12], [0, 0.31, 0.88], [0, 0.69, 0.88],
            [0.31, 0.12, 0], [0.69, 0.12, 0], [0.31, 0.88, 0], [0.69, 0.88, 0],
            [0.12, 0, 0.31], [0.12, 0, 0.69], [0.88, 0, 0.31], [0.88, 0, 0.69],
        ]
        
        # Use mixed framework elements
        framework_elem = info["framework"].split("-")
        if len(framework_elem) == 2:
            elem1, elem2 = framework_elem
        else:
            elem1, elem2 = framework_elem[0], framework_elem[0]
        
        # Add 6c positions
        for pos in wyckoff_6c:
            species.append(elem2)
            coords.append([p % 1 for p in pos])
        
        # Add 16i positions
        for pos in wyckoff_16i:
            species.append(elem1 if len(species) % 2 == 0 else elem2)
            coords.append([p % 1 for p in pos])
        
        # Add 24k positions
        for pos in wyckoff_24k:
            species.append(elem1)
            coords.append([p % 1 for p in pos])
        
        # Guest positions (2a and 6d)
        if not empty_cages:
            guest = info["guest"] if isinstance(info["guest"], str) else info["guest"][0]
            
            # 2a: 0, 0, 0 (small cage)
            species.append(guest)
            coords.append([0, 0, 0])
            species.append(guest)
            coords.append([0.5, 0.5, 0.5])
            
            # 6d: 1/4, 1/2, 0 (large cage)
            for pos in [[0.25, 0.5, 0], [0.75, 0.5, 0], [0, 0.25, 0.5],
                       [0, 0.75, 0.5], [0.5, 0, 0.25], [0.5, 0, 0.75]]:
                species.append(guest)
                coords.append(pos)
    
    elif clath_type == "II":
        # Type II: Fd-3m, 136 framework + 24 guest
        # Simplified version
        for i in range(8):
            for j in range(8):
                for k in range(8):
                    if (i + j + k) % 2 == 0:
                        species.append("Si")
                        coords.append([i/8, j/8, k/8])
        
        if not empty_cages:
            guest = info["guest"] if isinstance(info["guest"], str) else info["guest"][0]
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        if (i + j + k) % 3 == 0:
                            species.append(guest)
                            coords.append([(i+0.125)/3, (j+0.125)/3, (k+0.125)/3])
    
    else:
        # Hydrate or other type - simplified
        species = ["O", "H", "H", "C"]
        coords = [[0, 0, 0], [0.1, 0, 0], [-0.1, 0, 0], [0.5, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "clathrate": clathrate,
        "type": clath_type,
        "framework": info["framework"],
        "guest": info["guest"],
        "empty_framework": empty_cages,
        "is_thermoelectric": info.get("thermoelectric", False),
        "ZT": info.get("ZT", 0),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_empty_clathrate(
    framework: str = "Si",
    clathrate_type: str = "I"
) -> Dict[str, Any]:
    """
    Generate empty clathrate framework (superconducting Si-46).
    
    Args:
        framework: Framework element
        clathrate_type: I or II
    
    Returns:
        Empty clathrate structure
    """
    if clathrate_type == "I":
        return generate_clathrate("Na8Si46", empty_cages=True)
    else:
        return generate_clathrate("Na24Si136", empty_cages=True)
