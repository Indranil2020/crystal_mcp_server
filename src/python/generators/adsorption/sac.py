"""
adsorption/sac.py - Single Atom Catalyst Generation

Generates SAC structures on various supports.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


SAC_SUPPORTS = {
    "graphene": {"a": 2.46, "element": "C", "type": "2D"},
    "g-C3N4": {"a": 7.13, "element": "C", "type": "2D"},
    "C2N": {"a": 8.33, "element": "C", "type": "2D"},
    "MoS2": {"a": 3.16, "element": "Mo", "type": "2D"},
    "TiO2": {"a": 4.59, "c": 2.96, "element": "Ti", "type": "3D"},
    "CeO2": {"a": 5.41, "element": "Ce", "type": "3D"},
}


def generate_sac(
    metal: str = "Pt",
    support: str = "graphene",
    site: str = "hollow",
    size: List[int] = [3, 3, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate Single Atom Catalyst structure.
    
    Args:
        metal: SAC metal atom
        support: Support material
        site: Adsorption site type
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        SAC structure
    """
    if support not in SAC_SUPPORTS:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown support: {support}"}}
    
    info = SAC_SUPPORTS[support]
    a = info["a"]
    
    c = vacuum
    lattice = Lattice.hexagonal(a * size[0], c)
    
    species = []
    coords = []
    z_support = 0.4
    
    # Build support
    for i in range(size[0]):
        for j in range(size[1]):
            x = (1/3 + i) / size[0]
            y = (2/3 + j) / size[1]
            species.append(info["element"])
            coords.append([x, y, z_support])
    
    # Add SAC metal
    if site == "hollow":
        x_sac, y_sac = 0.5, 0.5
    elif site == "top":
        x_sac, y_sac = 1/(3*size[0]), 2/(3*size[1])
    else:
        x_sac, y_sac = 0.25, 0.25
    
    species.append(metal)
    coords.append([x_sac, y_sac, z_support + 2.5/c])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "support": support,
        "site": site,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"a": lattice.a, "c": lattice.c},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }
