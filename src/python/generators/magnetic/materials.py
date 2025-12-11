"""
magnetic/materials.py - Magnetic Material Structures

Ferromagnets, antiferromagnets, ferrimagnets with magnetic moments.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


MAGNETIC_DATABASE = {
    # Elemental Ferromagnets
    "Fe-bcc": {"a": 2.87, "structure": "bcc", "ordering": "FM", "moment": 2.2, "tc": 1043},
    "Co-hcp": {"a": 2.51, "c": 4.07, "structure": "hcp", "ordering": "FM", "moment": 1.7, "tc": 1388},
    "Ni-fcc": {"a": 3.52, "structure": "fcc", "ordering": "FM", "moment": 0.6, "tc": 627},
    "Gd-hcp": {"a": 3.64, "c": 5.78, "structure": "hcp", "ordering": "FM", "moment": 7.6, "tc": 293},
    
    # Antiferromagnets
    "Cr-bcc": {"a": 2.88, "structure": "bcc", "ordering": "AFM", "moment": 0.6, "tn": 311},
    "MnO": {"a": 4.44, "structure": "rocksalt", "ordering": "AFM", "moment": 4.6, "tn": 118},
    "NiO": {"a": 4.18, "structure": "rocksalt", "ordering": "AFM", "moment": 1.9, "tn": 523},
    "FeO": {"a": 4.30, "structure": "rocksalt", "ordering": "AFM", "moment": 4.6, "tn": 198},
    "CoO": {"a": 4.26, "structure": "rocksalt", "ordering": "AFM", "moment": 3.8, "tn": 289},
    "alpha-Fe2O3": {"a": 5.04, "c": 13.75, "structure": "corundum", "ordering": "AFM", "moment": 4.6, "tn": 955},
    
    # Ferrimagnets
    "Fe3O4": {"a": 8.396, "structure": "spinel", "ordering": "FiM", "moment": 4.1, "tc": 858},
    "CoFe2O4": {"a": 8.391, "structure": "spinel", "ordering": "FiM", "moment": 3.7, "tc": 793},
    "NiFe2O4": {"a": 8.339, "structure": "spinel", "ordering": "FiM", "moment": 2.4, "tc": 858},
    "YIG": {"a": 12.376, "structure": "garnet", "ordering": "FiM", "moment": 5.0, "tc": 559,
            "formula": "Y3Fe5O12", "description": "Yttrium iron garnet"},
    
    # Hard Magnets
    "Nd2Fe14B": {"a": 8.80, "c": 12.20, "structure": "tetragonal", "ordering": "FM", 
                 "moment": 2.15, "tc": 585, "BHmax": 512, "description": "Strongest permanent magnet"},
    "SmCo5": {"a": 5.01, "c": 3.97, "structure": "hexagonal", "ordering": "FM",
              "moment": 1.0, "tc": 1020, "BHmax": 240, "description": "High-temperature magnet"},
    "Sm2Co17": {"a": 8.40, "c": 12.20, "structure": "hexagonal", "ordering": "FM",
                "moment": 1.2, "tc": 1190, "BHmax": 280, "description": "High-coercivity"},
    
    # Soft Magnets
    "Fe-Si": {"a": 2.87, "structure": "bcc", "ordering": "FM", "moment": 2.0, "tc": 1000,
              "description": "Electrical steel"},
    "Permalloy": {"a": 3.55, "structure": "fcc", "ordering": "FM", "moment": 1.0, "tc": 870,
                  "formula": "Ni80Fe20", "description": "High permeability"},
    
    # Skyrmion Hosts
    "MnSi": {"a": 4.558, "structure": "B20", "ordering": "helical", "moment": 0.4, "tc": 29,
             "skyrmion": True, "description": "Chiral magnet with skyrmions"},
    "FeGe": {"a": 4.700, "structure": "B20", "ordering": "helical", "moment": 1.0, "tc": 278,
             "skyrmion": True, "description": "Room-temperature skyrmions"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_magnetic_material(
    material: str = "Fe-bcc",
    size: List[int] = [1, 1, 1],
    spin_configuration: str = "ground"
) -> Dict[str, Any]:
    """
    Generate magnetic material with moments.
    
    Args:
        material: Material name
        size: Supercell size
        spin_configuration: 'ground', 'FM', 'AFM', or 'random'
    
    Returns:
        Magnetic structure with moments
    """
    if material not in MAGNETIC_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {material}",
                "available": list(MAGNETIC_DATABASE.keys())}}
    
    info = MAGNETIC_DATABASE[material]
    a = info["a"]
    struct_type = info["structure"]
    ordering = info["ordering"]
    moment = info["moment"]
    
    if struct_type == "bcc":
        lattice = Lattice.bcc(a)
        elem = material.split("-")[0]
        species = [elem, elem]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        
    elif struct_type == "fcc":
        lattice = Lattice.fcc(a)
        elem = material.split("-")[0]
        species = [elem] * 4
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        
    elif struct_type == "hcp":
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        elem = material.split("-")[0]
        species = [elem, elem]
        coords = [[1/3, 2/3, 0.25], [2/3, 1/3, 0.75]]
        
    elif struct_type == "rocksalt":
        lattice = Lattice.cubic(a)
        M = material[:-1]
        species = [M, M, M, M, "O", "O", "O", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        
    elif struct_type == "spinel":
        lattice = Lattice.cubic(a)
        species = ["Fe"] * 8 + ["O"] * 4
        coords = [[i/8, i/8, i/8] for i in range(8)] + [[0.25, 0.25, 0.25 + i/4] for i in range(4)]
        
    elif struct_type == "B20":
        lattice = Lattice.cubic(a)
        M, X = material[:2], material[2:]
        species = [M, M, M, M, X, X, X, X]
        coords = [[0, 0, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0],
                  [0.137, 0.137, 0.137], [0.637, 0.363, 0.863], [0.363, 0.863, 0.637], [0.863, 0.637, 0.363]]
        
    else:
        c = info.get("c", a)
        lattice = Lattice.tetragonal(a, c) if c != a else Lattice.cubic(a)
        species = ["Fe"] * 4
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Add magnetic moments
    struct_dict = structure_to_dict(structure)
    
    for i, atom in enumerate(struct_dict["atoms"]):
        if spin_configuration == "ground":
            if ordering == "FM":
                atom["magmom"] = [0, 0, moment]
            elif ordering == "AFM":
                sign = 1 if i % 2 == 0 else -1
                atom["magmom"] = [0, 0, sign * moment]
            elif ordering == "FiM":
                sign = 1 if i < len(struct_dict["atoms"]) // 2 else -0.5
                atom["magmom"] = [0, 0, sign * moment]
            else:
                atom["magmom"] = [0, 0, moment]
        elif spin_configuration == "FM":
            atom["magmom"] = [0, 0, moment]
        elif spin_configuration == "AFM":
            sign = 1 if i % 2 == 0 else -1
            atom["magmom"] = [0, 0, sign * moment]
        else:
            np.random.seed(i)
            atom["magmom"] = list(np.random.randn(3) * moment)
    
    return {
        "success": True,
        "material": material,
        "ordering": ordering,
        "moment_uB": moment,
        "tc_kelvin": info.get("tc"),
        "tn_kelvin": info.get("tn"),
        "is_skyrmion_host": info.get("skyrmion", False), 
        "structure": struct_dict
    }
