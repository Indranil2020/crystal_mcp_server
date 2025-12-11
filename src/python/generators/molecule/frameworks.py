"""
molecule/frameworks.py - MOF and COF Generation

Generates framework structures:
- Metal-Organic Frameworks (MOFs)
- Covalent Organic Frameworks (COFs)
- Zeolitic Imidazolate Frameworks (ZIFs)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


MOF_DATABASE = {
    "MOF-5": {"formula": "Zn4O(BDC)3", "a": 25.83, "topology": "pcu", "surface_area": 3800},
    "HKUST-1": {"formula": "Cu3(BTC)2", "a": 26.34, "topology": "tbo", "surface_area": 1850},
    "UiO-66": {"formula": "Zr6O4(OH)4(BDC)6", "a": 20.74, "topology": "fcu", "surface_area": 1200},
    "ZIF-8": {"formula": "Zn(mIm)2", "a": 16.99, "topology": "sod", "surface_area": 1947},
    "MIL-101": {"formula": "Cr3O(BDC)3", "a": 89.0, "topology": "mtn", "surface_area": 4100},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_mof(
    mof_name: str,
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate Metal-Organic Framework structure.
    
    Args:
        mof_name: MOF name (MOF-5, HKUST-1, UiO-66, ZIF-8)
        size: Supercell size
    
    Returns:
        MOF structure
    """
    if mof_name not in MOF_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown MOF: {mof_name}", "available": list(MOF_DATABASE.keys())}}
    
    info = MOF_DATABASE[mof_name]
    a = info["a"]
    
    lattice = Lattice.cubic(a)
    
    # Simplified representation (metal nodes only)
    if mof_name == "MOF-5":
        species = ["Zn"] * 4 + ["O"]
        coords = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.25, 0.25, 0.25]]
    elif mof_name == "HKUST-1":
        species = ["Cu"] * 3
        coords = [[0, 0, 0], [0.33, 0.33, 0.33], [0.67, 0.67, 0.67]]
    elif mof_name == "ZIF-8":
        species = ["Zn"] * 2 + ["N"] * 4
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0, 0.5], [0.75, 0, 0.5], [0, 0.25, 0.5], [0, 0.75, 0.5]]
    else:
        species = ["Zr"] * 6
        coords = [[i/6, i/6, i/6] for i in range(6)]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "mof_name": mof_name,
        "formula": info["formula"],
        "topology": info["topology"],
        "surface_area_m2_g": info["surface_area"],
        "structure": structure_to_dict(structure)
    }


def generate_cof(
    topology: str = "hcb",
    linker: str = "BDC",
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate Covalent Organic Framework.
    
    Args:
        topology: COF topology (hcb, sql, kgm)
        linker: Organic linker
        size: Supercell size
    
    Returns:
        COF structure
    """
    topologies = {
        "hcb": {"a": 20.0, "description": "Honeycomb"},
        "sql": {"a": 15.0, "description": "Square lattice"},
        "kgm": {"a": 25.0, "description": "Kagome"},
    }
    
    if topology not in topologies:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown topology"}}
    
    info = topologies[topology]
    a = info["a"]
    c = 3.5  # Layer separation
    
    lattice = Lattice.hexagonal(a, c)
    
    # Simplified node positions
    species = ["C"] * 6
    coords = [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5], [0, 0, 0.5], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "topology": topology,
        "linker": linker,
        "description": info["description"],
        "structure": structure_to_dict(structure)
    }
