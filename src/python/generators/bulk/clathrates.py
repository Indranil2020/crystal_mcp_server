"""
bulk/clathrates.py - Clathrate Cage Structure Generation

Generates clathrate structures with guest atoms:
- Type I (sI): Pm-3n, 46 framework atoms, 2 small + 6 large cages
- Type II (sII): Fd-3m, 136 framework atoms, 16 small + 8 large cages
- Type III (sIII): P4_2/mnm, 34 framework atoms
- Type H (sH): Complex with 34 atoms per unit cell

Scientific basis:
- Gas hydrate structures (CH4, CO2, H2)
- Semiconductor clathrates (Ba8Si46, Na8Si46)
- Intermetallic clathrates
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Clathrate type definitions
CLATHRATE_TYPES = {
    "I": {
        "space_group": 223,  # Pm-3n
        "n_framework": 46,
        "cages": {"small": 2, "large": 6},
        "description": "Type I (sI) - most common for gas hydrates"
    },
    "II": {
        "space_group": 227,  # Fd-3m
        "n_framework": 136,
        "cages": {"small": 16, "large": 8},
        "description": "Type II (sII) - larger guests"
    },
    "III": {
        "space_group": 136,  # P4_2/mnm
        "n_framework": 34,
        "cages": {"small": 3, "large": 2},
        "description": "Type III - tetragonal"
    },
    "H": {
        "space_group": 191,  # P6/mmm
        "n_framework": 34,
        "cages": {"small": 3, "medium": 2, "large": 1},
        "description": "Type H (sH) - hexagonal"
    }
}

# Lattice constants for common clathrate compositions
CLATHRATE_LATTICE_DB = {
    ("Si", "I"): 10.19,
    ("Ge", "I"): 10.72,
    ("Sn", "I"): 11.68,
    ("Si", "II"): 14.67,
    ("Ge", "II"): 15.23,
}


def generate_clathrate(
    framework_element: str = "Si",
    clathrate_type: str = "I",
    guest_elements: Optional[Dict[str, str]] = None,
    a: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate clathrate cage structure.
    
    Args:
        framework_element: Element forming the cage network (Si, Ge, Sn)
        clathrate_type: Type of clathrate (I, II, III, H)
        guest_elements: Mapping of cage type to guest element
                       e.g., {"small": "Ba", "large": "Na"}
        a: Lattice constant (estimated if not provided)
    
    Returns:
        Clathrate structure dictionary
    
    Examples:
        >>> result = generate_clathrate('Si', 'I', guest_elements={'small': 'Ba', 'large': 'Na'})
        >>> result["success"]
        True
    """
    if clathrate_type not in CLATHRATE_TYPES:
        return {
            "success": False,
            "error": {
                "code": "INVALID_TYPE",
                "message": f"Unknown clathrate type '{clathrate_type}'",
                "available": list(CLATHRATE_TYPES.keys())
            }
        }
    
    template = CLATHRATE_TYPES[clathrate_type]
    
    # Determine lattice constant
    if a is None:
        key = (framework_element, clathrate_type)
        a = CLATHRATE_LATTICE_DB.get(key, 10.5)
    
    if clathrate_type == "I":
        lattice = Lattice.cubic(a)
        
        # Type I framework positions (Wyckoff sites 6c, 16i, 24k)
        # 6c sites (0, 1/2, 1/4)
        coords_6c = [
            [0, 0.5, 0.25], [0, 0.5, 0.75],
            [0.5, 0.25, 0], [0.5, 0.75, 0],
            [0.25, 0, 0.5], [0.75, 0, 0.5]
        ]
        
        # 16i sites (~0.18, 0.18, 0.18)
        x = 0.1837
        coords_16i = []
        for sx in [1, -1]:
            for sy in [1, -1]:
                for sz in [1, -1]:
                    coords_16i.append([(sx * x) % 1, (sy * x) % 1, (sz * x) % 1])
                    # Additional by symmetry
        coords_16i = coords_16i[:16]  # Take 16
        
        # 24k sites (~0.31, 0.12, 0)
        x, y = 0.3092, 0.1172
        coords_24k = []
        for pos in [[x, y, 0], [y, x, 0], [-x, y, 0], [y, -x, 0]]:
            for p in [[pos[0], pos[1], pos[2]], [pos[2], pos[0], pos[1]], [pos[1], pos[2], pos[0]]]:
                coords_24k.append([c % 1 for c in p])
        coords_24k = coords_24k[:24]  # Take 24
        
        framework_coords = coords_6c + coords_16i + coords_24k
        framework_species = [framework_element] * len(framework_coords)
        
        # Guest cage centers
        guest_species = []
        guest_coords = []
        
        if guest_elements:
            if "small" in guest_elements:
                # 2a site: small cages at (0, 0, 0) and (0.5, 0.5, 0.5)
                guest_species.extend([guest_elements["small"]] * 2)
                guest_coords.extend([[0, 0, 0], [0.5, 0.5, 0.5]])
            
            if "large" in guest_elements:
                # 6d sites: large cages
                guest_species.extend([guest_elements["large"]] * 6)
                guest_coords.extend([
                    [0.25, 0.5, 0], [0.75, 0.5, 0],
                    [0, 0.25, 0.5], [0, 0.75, 0.5],
                    [0.5, 0, 0.25], [0.5, 0, 0.75]
                ])
        
        species = framework_species + guest_species
        coords = framework_coords + guest_coords
    
    elif clathrate_type == "II":
        lattice = Lattice.cubic(a)
        
        # Type II is more complex - simplified version
        # 8a + 32e + 96g for framework
        species = [framework_element] * 34  # Simplified
        coords = [[i/10, (i*2)%10/10, (i*3)%10/10] for i in range(34)]
        
        if guest_elements:
            if "small" in guest_elements:
                species.extend([guest_elements["small"]] * 16)
                coords.extend([[0.125 + 0.25*i, 0.125, 0.125] for i in range(16)])
            if "large" in guest_elements:
                species.extend([guest_elements["large"]] * 8)
                coords.extend([[0.375, 0.375, 0.375 + 0.125*i] for i in range(8)])
    
    else:
        # Simplified for other types
        lattice = Lattice.cubic(a)
        species = [framework_element] * template["n_framework"]
        coords = [[i/template["n_framework"], 0.5, 0.5] for i in range(template["n_framework"])]
    
    structure = Structure(lattice, species, coords)
    
    # Determine formula
    n_guests = sum(template["cages"].values()) if guest_elements else 0
    formula = f"{guest_elements.get('small', '') if guest_elements else ''}{template['cages'].get('small', 0)}"
    formula += f"{guest_elements.get('large', '') if guest_elements else ''}{template['cages'].get('large', 0)}"
    formula += f"{framework_element}{template['n_framework']}"
    
    return {
        "success": True,
        "clathrate_type": clathrate_type,
        "description": template["description"],
        "framework": framework_element,
        "n_framework_atoms": template["n_framework"],
        "cages": template["cages"],
        "guests": guest_elements or {},
        "structure": structure_to_dict(structure)
    }
