"""
catalyst/saa.py - Single Atom Alloy Catalysts
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


SAA_DATABASE = {
    "Pt/Cu": {"host": "Cu", "guest": "Pt", "host_a": 3.61, "guest_a": 3.92},
    "Pd/Au": {"host": "Au", "guest": "Pd", "host_a": 4.08, "guest_a": 3.89},
    "Rh/Ag": {"host": "Ag", "guest": "Rh", "host_a": 4.09, "guest_a": 3.80},
    "Ni/Cu": {"host": "Cu", "guest": "Ni", "host_a": 3.61, "guest_a": 3.52},
    "Pt/Ag": {"host": "Ag", "guest": "Pt", "host_a": 4.09, "guest_a": 3.92},
}


def generate_saa(
    guest_metal: str = "Pt",
    host_metal: str = "Cu",
    surface: str = "111",
    size: List[int] = [3, 3, 4],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate Single Atom Alloy structure.
    
    Args:
        guest_metal: Isolated guest atom
        host_metal: Host metal
        surface: Surface facet
        size: Slab size
        vacuum: Vacuum padding
    
    Returns:
        SAA structure
    """
    key = f"{guest_metal}/{host_metal}"
    if key in SAA_DATABASE:
        info = SAA_DATABASE[key]
        a = info["host_a"]
    else:
        a = 3.6  # Default FCC
    
    # Build FCC (111) slab
    a_surf = a / np.sqrt(2)
    c = size[2] * a / np.sqrt(3) + vacuum
    lattice = Lattice.hexagonal(a_surf * size[0], c)
    
    species = []
    coords = []
    
    n_layers = size[2]
    z_spacing = a / np.sqrt(3) / c
    
    for layer in range(n_layers):
        z = layer * z_spacing + vacuum / 2 / c
        shift = (layer % 3) / 3
        
        for i in range(size[0]):
            for j in range(size[1]):
                x = (i + shift) / size[0]
                y = (j + shift) / size[1]
                species.append(host_metal)
                coords.append([x % 1, y % 1, z])
    
    # Replace one surface atom with guest
    surface_atoms = [i for i, c in enumerate(coords) if c[2] > (n_layers - 1.5) * z_spacing]
    if surface_atoms:
        species[surface_atoms[len(surface_atoms) // 2]] = guest_metal
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "guest_metal": guest_metal,
        "host_metal": host_metal,
        "surface": surface,
        "n_atoms": len(structure),
        "guest_concentration": 1 / len(surface_atoms) if surface_atoms else 0,
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }
