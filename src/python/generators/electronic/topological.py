"""
electronic/topological.py - Topological Materials

Topological insulators, Weyl/Dirac semimetals, nodal-line materials.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


TOPOLOGICAL_DATABASE = {
    # 3D Topological Insulators
    "Bi2Se3": {"a": 4.143, "c": 28.636, "structure": "rhombohedral", "type": "TI", 
               "gap_eV": 0.3, "description": "Prototypical 3D TI"},
    "Bi2Te3": {"a": 4.385, "c": 30.497, "structure": "rhombohedral", "type": "TI",
               "gap_eV": 0.15, "description": "3D TI, thermoelectric"},
    "Sb2Te3": {"a": 4.264, "c": 30.458, "structure": "rhombohedral", "type": "TI",
               "gap_eV": 0.21, "description": "3D TI"},
    
    # Weyl Semimetals
    "TaAs": {"a": 3.437, "c": 11.656, "structure": "tetragonal", "type": "Weyl-I",
             "gap_eV": 0, "n_weyl_points": 24, "description": "Type-I Weyl semimetal"},
    "NbAs": {"a": 3.452, "c": 11.680, "structure": "tetragonal", "type": "Weyl-I",
             "gap_eV": 0, "n_weyl_points": 24, "description": "Type-I Weyl semimetal"},
    "TaP": {"a": 3.318, "c": 11.348, "structure": "tetragonal", "type": "Weyl-I",
            "gap_eV": 0, "n_weyl_points": 24, "description": "Type-I Weyl semimetal"},
    "WTe2": {"a": 3.477, "b": 6.249, "c": 14.018, "structure": "orthorhombic", "type": "Weyl-II",
             "gap_eV": 0, "description": "Type-II Weyl semimetal"},
    "MoTe2-Td": {"a": 3.469, "b": 6.310, "c": 13.861, "structure": "orthorhombic", "type": "Weyl-II",
                 "gap_eV": 0, "description": "Type-II Weyl semimetal"},
    
    # Dirac Semimetals
    "Cd3As2": {"a": 12.633, "c": 25.427, "structure": "tetragonal", "type": "Dirac",
               "gap_eV": 0, "description": "3D Dirac semimetal"},
    "Na3Bi": {"a": 5.448, "c": 9.655, "structure": "hexagonal", "type": "Dirac",
              "gap_eV": 0, "description": "3D Dirac semimetal"},
    
    # Nodal-line Semimetals
    "ZrSiS": {"a": 3.544, "c": 8.055, "structure": "tetragonal", "type": "nodal-line",
              "gap_eV": 0, "description": "Nodal-line semimetal"},
    "PbTaSe2": {"a": 3.439, "c": 9.416, "structure": "hexagonal", "type": "nodal-line",
                "gap_eV": 0, "description": "Nodal-line semimetal with SOC"},
    
    # Higher-order TIs
    "Bi": {"a": 4.546, "c": 11.862, "structure": "rhombohedral", "type": "HOTI",
           "gap_eV": 0.01, "description": "Higher-order TI"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_topological_insulator(
    material: str = "Bi2Se3",
    size: List[int] = [1, 1, 1],
    n_quintuple_layers: int = 1
) -> Dict[str, Any]:
    """
    Generate topological material structure.
    
    Args:
        material: Material name
        size: Supercell size
        n_quintuple_layers: Number of QL for films
    
    Returns:
        Topological material structure
    """
    if material not in TOPOLOGICAL_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {material}",
                "available": list(TOPOLOGICAL_DATABASE.keys())}}
    
    info = TOPOLOGICAL_DATABASE[material]
    a = info["a"]
    c = info["c"]
    struct_type = info["structure"]
    
    if struct_type == "rhombohedral" and "Bi2" in material:
        # Bi2X3 quintuple layer structure
        lattice = Lattice.hexagonal(a, c)
        X = material.split("2")[1][:-1]  # Se or Te
        # Simplified 1 QL
        species = ["Bi", "Bi", X, X, X]
        coords = [[0, 0, 0.4], [0, 0, 0.6], [0, 0, 0], [0, 0, 0.21], [0, 0, 0.79]]
    
    elif struct_type == "tetragonal" and info["type"].startswith("Weyl"):
        lattice = Lattice.tetragonal(a, c)
        if "Ta" in material or "Nb" in material:
            M = material[:2]
            X = material[2:]
            species = [M, M, M, M, X, X, X, X]
            coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0, 0.5, 0.25], [0.5, 0, 0.75],
                      [0, 0, 0.41], [0, 0, 0.59], [0.5, 0.5, 0.09], [0.5, 0.5, 0.91]]
        else:
            species = ["Zr", "Si", "S", "S"]
            coords = [[0, 0, 0], [0, 0, 0.5], [0, 0.5, 0.25], [0, 0.5, 0.75]]
    
    elif struct_type == "orthorhombic":
        b = info["b"]
        lattice = Lattice.orthorhombic(a, b, c)
        M = "W" if "W" in material else "Mo"
        species = [M, M, "Te", "Te", "Te", "Te"]
        coords = [[0, 0.25, 0.5], [0, 0.75, 0.5],
                  [0.25, 0.35, 0.6], [0.75, 0.65, 0.4], [0.25, 0.65, 0.6], [0.75, 0.35, 0.4]]
    
    elif struct_type == "hexagonal":
        lattice = Lattice.hexagonal(a, c)
        if material == "Na3Bi":
            species = ["Na", "Na", "Na", "Bi"]
            coords = [[0, 0, 0], [0.33, 0.67, 0.25], [0.67, 0.33, 0.75], [0.33, 0.67, 0.75]]
        else:
            species = ["Pb", "Ta", "Se", "Se"]
            coords = [[0, 0, 0], [0.33, 0.67, 0.5], [0.33, 0.67, 0.15], [0.33, 0.67, 0.85]]
    
    else:
        lattice = Lattice.hexagonal(a, c)
        species = ["Bi", "Bi"]
        coords = [[0, 0, 0.23], [0, 0, 0.77]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": material,
        "topological_type": info["type"],
        "gap_eV": info["gap_eV"],
        "n_weyl_points": info.get("n_weyl_points"),
        "description": info["description"],
        "structure": structure_to_dict(structure)
    }
