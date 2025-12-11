"""
thermoelectric/materials.py - Thermoelectric Material Generation

Comprehensive thermoelectric materials with ZT values.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


THERMOELECTRIC_DATABASE = {
    # Bismuth Tellurides
    "Bi2Te3": {"a": 4.385, "c": 30.497, "structure": "rhombohedral", "zt_max": 1.0, 
               "t_optimal": 300, "type": "n/p", "description": "Commercial standard"},
    "Sb2Te3": {"a": 4.264, "c": 30.458, "structure": "rhombohedral", "zt_max": 0.8,
               "t_optimal": 300, "type": "p", "description": "p-type telluride"},
    "Bi2Se3": {"a": 4.143, "c": 28.636, "structure": "rhombohedral", "zt_max": 0.5,
               "t_optimal": 400, "type": "n", "description": "n-type selenide"},
    
    # Lead Chalcogenides
    "PbTe": {"a": 6.462, "structure": "rocksalt", "zt_max": 2.2, 
             "t_optimal": 800, "type": "n/p", "description": "Mid-temp champion"},
    "PbSe": {"a": 6.128, "structure": "rocksalt", "zt_max": 1.6,
             "t_optimal": 850, "type": "n/p", "description": "Mid-temp material"},
    "PbS": {"a": 5.936, "structure": "rocksalt", "zt_max": 1.2,
            "t_optimal": 900, "type": "n/p", "description": "Earth-abundant"},
    
    # Skutterudites
    "CoSb3": {"a": 9.039, "structure": "skutterudite", "zt_max": 1.0,
              "t_optimal": 800, "type": "n", "description": "Unfilled skutterudite"},
    "Ba0.3Co4Sb12": {"a": 9.048, "structure": "filled_skutterudite", "zt_max": 1.4,
                     "t_optimal": 800, "type": "n", "description": "Ba-filled"},
    "Yb0.2Co4Sb12": {"a": 9.043, "structure": "filled_skutterudite", "zt_max": 1.6,
                     "t_optimal": 800, "type": "n", "description": "Yb-filled, high ZT"},
    
    # Half-Heuslers
    "ZrNiSn": {"a": 6.115, "structure": "half_heusler", "zt_max": 1.0,
               "t_optimal": 900, "type": "n", "description": "n-type half-Heusler"},
    "HfNiSn": {"a": 6.085, "structure": "half_heusler", "zt_max": 1.0,
               "t_optimal": 900, "type": "n", "description": "Hf-based"},
    "TiNiSn": {"a": 5.930, "structure": "half_heusler", "zt_max": 0.7,
               "t_optimal": 800, "type": "n", "description": "Ti-based"},
    "ZrCoSb": {"a": 6.074, "structure": "half_heusler", "zt_max": 0.8,
               "t_optimal": 1000, "type": "p", "description": "p-type half-Heusler"},
    
    # Clathrates
    "Ba8Ga16Ge30": {"a": 10.78, "structure": "clathrate_I", "zt_max": 1.35,
                    "t_optimal": 900, "type": "n", "description": "Type-I clathrate"},
    "Sr8Ga16Ge30": {"a": 10.72, "structure": "clathrate_I", "zt_max": 1.0,
                    "t_optimal": 800, "type": "n", "description": "Sr clathrate"},
    
    # High-entropy
    "GeSe": {"a": 4.40, "b": 3.83, "c": 10.83, "structure": "orthorhombic", "zt_max": 2.5,
             "t_optimal": 700, "type": "p", "description": "Layered, record ZT"},
    "SnSe": {"a": 4.44, "b": 4.15, "c": 11.49, "structure": "orthorhombic", "zt_max": 2.6,
             "t_optimal": 923, "type": "p", "description": "Record holder at 923K"},
    "Cu2Se": {"a": 5.76, "structure": "fcc", "zt_max": 2.1,
              "t_optimal": 1000, "type": "p", "description": "Superionic conductor"},
    
    # Oxides
    "SrTiO3": {"a": 3.905, "structure": "perovskite", "zt_max": 0.4,
               "t_optimal": 1000, "type": "n", "description": "Oxide thermoelectric"},
    "CaMnO3": {"a": 5.28, "b": 5.27, "c": 7.46, "structure": "perovskite", "zt_max": 0.3,
               "t_optimal": 1000, "type": "n", "description": "Mn oxide"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_thermoelectric(
    material: str = "Bi2Te3",
    size: List[int] = [1, 1, 1],
    doping_type: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate thermoelectric material structure.
    
    Args:
        material: Material name
        size: Supercell size
        doping_type: 'n' or 'p' type
    
    Returns:
        Thermoelectric structure with ZT
    """
    if material not in THERMOELECTRIC_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {material}",
                "available": list(THERMOELECTRIC_DATABASE.keys())}}
    
    info = THERMOELECTRIC_DATABASE[material]
    a = info["a"]
    struct_type = info["structure"]
    
    if struct_type == "rhombohedral":
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        X = material.split("2")[1][:-1]
        M = material.split("2")[0]
        species = [M, M, X, X, X]
        coords = [[0, 0, 0.4], [0, 0, 0.6], [0, 0, 0], [0, 0, 0.21], [0, 0, 0.79]]
    
    elif struct_type == "rocksalt":
        lattice = Lattice.cubic(a)
        M, X = material[:-2], material[-2:]
        if X == "Te":
            X = "Te"
        elif X == "Se":
            X = "Se"
        else:
            X = "S"
        species = [M, M, M, M, X, X, X, X]
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
    
    elif "skutterudite" in struct_type:
        lattice = Lattice.cubic(a)
        species = ["Co"] * 8 + ["Sb"] * 24
        cobalt_coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.25, 0.75, 0.75],
                         [0.75, 0.25, 0.75], [0.75, 0.75, 0.25], [0.25, 0.25, 0.75],
                         [0.25, 0.75, 0.25], [0.75, 0.25, 0.25]]
        coords = cobalt_coords + [[i/24, (i+5)/24, (i+10)/24] for i in range(24)]
    
    elif struct_type == "half_heusler":
        lattice = Lattice.cubic(a)
        # XYZ structure
        X, Y, Z = material[:-4], material[-4:-2], material[-2:]
        species = [X, Y, Z, Z]
        coords = [[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]]
    
    elif struct_type == "orthorhombic":
        b = info["b"]
        c = info["c"]
        lattice = Lattice.orthorhombic(a, b, c)
        M, X = material[:-2], material[-2:]
        species = [M, M, X, X]
        coords = [[0, 0, 0.12], [0.5, 0, 0.38], [0, 0.5, 0.15], [0.5, 0.5, 0.35]]
    
    else:
        lattice = Lattice.cubic(a)
        species = [material[0], material[0]]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": material,
        "zt_max": info["zt_max"],
        "t_optimal_K": info["t_optimal"],
        "carrier_type": info["type"],
        "description": info["description"],
        "structure": structure_to_dict(structure)
    }
