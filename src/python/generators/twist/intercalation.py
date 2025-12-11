"""
twist/intercalation.py - Intercalated Bilayer Structures

Generates intercalated 2D bilayer structures:
- Li, Na, K intercalated graphite
- Stage-1, stage-2 intercalation
- Ca, Fe, Co, Cu, Ag, Au intercalation
- Organic molecule intercalation

Scientific basis:
- Daumas-Hérold staging model
- Intercalant-layer interactions
- Battery electrode models
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Intercalant database
INTERCALANTS = {
    "Li": {"radius": 0.76, "charge": 1, "spacing_increase": 0.35, "sites": "hollow"},
    "Na": {"radius": 1.02, "charge": 1, "spacing_increase": 0.72, "sites": "hollow"},
    "K": {"radius": 1.38, "charge": 1, "spacing_increase": 1.35, "sites": "hollow"},
    "Ca": {"radius": 1.00, "charge": 2, "spacing_increase": 0.65, "sites": "hollow"},
    "Fe": {"radius": 0.78, "charge": 2, "spacing_increase": 0.40, "sites": "above_C"},
    "Co": {"radius": 0.75, "charge": 2, "spacing_increase": 0.35, "sites": "above_C"},
    "Cu": {"radius": 0.73, "charge": 1, "spacing_increase": 0.30, "sites": "bridge"},
    "Ag": {"radius": 1.15, "charge": 1, "spacing_increase": 0.85, "sites": "bridge"},
    "Au": {"radius": 1.37, "charge": 1, "spacing_increase": 1.10, "sites": "bridge"},
}


def generate_intercalated_bilayer(
    host_material: str = "graphene",
    intercalant: str = "Li",
    stage: int = 1,
    concentration: float = 1.0,
    supercell_size: int = 3,
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate intercalated bilayer structure.
    
    Args:
        host_material: Host 2D material (graphene, MoS2)
        intercalant: Intercalant element
        stage: Staging number (1 = every gallery, 2 = every other)
        concentration: Intercalant concentration (0-1)
        supercell_size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Intercalated bilayer structure
    """
    if intercalant not in INTERCALANTS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_INTERCALANT",
                "message": f"Unknown intercalant '{intercalant}'",
                "available": list(INTERCALANTS.keys())
            }
        }
    
    int_info = INTERCALANTS[intercalant]
    
    if host_material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep_base = 3.35
    else:
        a = 3.16
        element = "Mo"
        z_sep_base = 6.5
    
    z_sep = z_sep_base + int_info["spacing_increase"]
    
    supercell_size = min(supercell_size, 10)
    c = z_sep + vacuum
    lattice = Lattice.hexagonal(a * supercell_size, c)
    
    species = []
    coords = []
    z_center = 0.5
    dz = z_sep / c / 2
    
    # Bottom layer
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size
            y1 = (2/3 + j) / supercell_size
            x2 = (2/3 + i) / supercell_size
            y2 = (1/3 + j) / supercell_size
            
            species.extend([element, element])
            coords.extend([[x1, y1, z_center - dz], [x2, y2, z_center - dz]])
    
    # Top layer
    for i in range(supercell_size):
        for j in range(supercell_size):
            x1 = (1/3 + i) / supercell_size
            y1 = (2/3 + j) / supercell_size
            x2 = (2/3 + i) / supercell_size
            y2 = (1/3 + j) / supercell_size
            
            species.extend([element, element])
            coords.extend([[x1, y1, z_center + dz], [x2, y2, z_center + dz]])
    
    # Intercalant layer
    n_intercalant_sites = supercell_size ** 2
    n_intercalants = int(concentration * n_intercalant_sites)
    
    site_indices = list(range(n_intercalant_sites))
    np.random.seed(42)
    np.random.shuffle(site_indices)
    
    for idx in site_indices[:n_intercalants]:
        i = idx // supercell_size
        j = idx % supercell_size
        
        # Place at hollow site (center of hexagon)
        if int_info["sites"] == "hollow":
            x = (i + 0.5) / supercell_size
            y = (j + 0.5) / supercell_size
        else:  # above_C or bridge
            x = (1/3 + i) / supercell_size
            y = (2/3 + j) / supercell_size
        
        species.append(intercalant)
        coords.append([x, y, z_center])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate composition
    formula = f"{element}{supercell_size**2 * 4}{intercalant}{n_intercalants}"
    
    return {
        "success": True,
        "host_material": host_material,
        "intercalant": intercalant,
        "stage": stage,
        "concentration": round(n_intercalants / n_intercalant_sites, 3),
        "n_intercalants": n_intercalants,
        "interlayer_expansion_angstrom": int_info["spacing_increase"],
        "formula": formula,
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_staged_graphite(
    intercalant: str = "Li",
    stage: int = 1,
    n_galleries: int = 4,
    vacuum: float = 20.0
) -> Dict[str, Any]:
    """
    Generate staged intercalated graphite.
    
    Args:
        intercalant: Intercalant element
        stage: Staging number
        n_galleries: Number of graphene galleries
        vacuum: Vacuum padding
    
    Returns:
        Staged graphite structure
    """
    if intercalant not in INTERCALANTS:
        return {"success": False, "error": {"code": "INVALID", "message": f"Unknown: {intercalant}"}}
    
    int_info = INTERCALANTS[intercalant]
    
    a = 2.46
    z_sep_empty = 3.35
    z_sep_filled = z_sep_empty + int_info["spacing_increase"]
    
    # Calculate total height
    n_filled = n_galleries // stage
    n_empty = n_galleries - n_filled
    total_height = n_filled * z_sep_filled + n_empty * z_sep_empty
    
    c = total_height + vacuum
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    z = vacuum / 2 / c
    
    for gallery_idx in range(n_galleries + 1):
        # Graphene layer
        species.extend(["C", "C"])
        coords.extend([[1/3, 2/3, z], [2/3, 1/3, z]])
        
        if gallery_idx < n_galleries:
            # Check if this gallery is filled
            is_filled = (gallery_idx % stage == 0)
            if is_filled:
                z_inc = z_sep_filled / c
                # Add intercalant
                species.append(intercalant)
                coords.append([0, 0, z + z_inc/2])
            else:
                z_inc = z_sep_empty / c
            z += z_inc
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "intercalant": intercalant,
        "stage": stage,
        "n_galleries": n_galleries,
        "n_intercalants": n_filled,
        "structure": structure_to_dict(structure, vacuum)
    }
