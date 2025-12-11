"""
two_d/porous.py - Porous 2D Materials

Generates porous 2D materials:
- g-C3N4 (graphitic carbon nitride)
- Holey graphene
- Biphenylene network
- Covalent organic frameworks (2D COFs)
- Graphyne variants (α, β, γ, δ)

Scientific basis:
- Regular pore patterns
- sp-sp2 carbon networks
- Photocatalytic materials
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Porous 2D materials database
POROUS_2D_MATERIALS = {
    "g-C3N4": {
        "formula": "C3N4",
        "a": 7.13,
        "description": "Graphitic carbon nitride (tri-s-triazine)",
        "pore_size": 3.2,
        "applications": ["photocatalysis", "CO2 reduction"]
    },
    "C2N": {
        "formula": "C2N",
        "a": 8.33,
        "description": "Holey graphene-like, Nitrogenated",
        "pore_size": 5.5,
        "applications": ["gas separation", "Li-ion battery"]
    },
    "biphenylene": {
        "formula": "C",
        "a": 3.73,
        "b": 4.52,
        "description": "4-6-8 carbon ring network",
        "pore_size": 0,  # No true pores
        "applications": ["electronics", "Li storage"]
    },
    "graphyne-alpha": {
        "formula": "C",
        "a": 6.86,
        "description": "α-graphyne (sp-sp2 carbon)",
        "pore_size": 3.8,
        "applications": ["gas separation"]
    },
    "graphyne-gamma": {
        "formula": "C",
        "a": 6.86,
        "description": "γ-graphyne (graphdiyne-like)",
        "pore_size": 5.2,
        "applications": ["gas separation", "electronics"]
    },
    "graphdiyne": {
        "formula": "C",
        "a": 9.46,
        "description": "Graphdiyne with diacetylene links",
        "pore_size": 5.5,
        "applications": ["Li/Na storage", "catalysis"]
    },
    "holey_graphene": {
        "formula": "C",
        "a": 7.38,
        "description": "Graphene with periodic holes",
        "pore_size": 2.5,
        "applications": ["desalination", "gas separation"]
    }
}


def generate_porous_2d(
    material: str,
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate porous 2D material structure.
    
    Args:
        material: Material name (g-C3N4, C2N, biphenylene, etc.)
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Porous 2D structure
    """
    material_key = material.replace(" ", "_")
    
    # Match by key or formula
    matched = None
    for key, info in POROUS_2D_MATERIALS.items():
        if key.lower() == material_key.lower() or key.lower() == material.lower() or info["formula"].lower() == material.lower():
            matched = key
            break
    
    if not matched:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MATERIAL",
                "message": f"Unknown porous 2D material '{material}'",
                "available": list(POROUS_2D_MATERIALS.keys())
            }
        }
    
    params = POROUS_2D_MATERIALS[matched]
    a = params["a"]
    b = params.get("b", a)
    c = vacuum
    
    if matched == "g-C3N4":
        # Graphitic C3N4 (tri-s-triazine / heptazine)
        lattice = Lattice.hexagonal(a, c)
        # Simplified structure with C3N4 stoichiometry
        species = ["C", "C", "C", "N", "N", "N", "N"]
        z = 0.5
        coords = [
            [0.0, 0.0, z],
            [1/3, 0.0, z],
            [2/3, 0.0, z],
            [0.17, 0.17, z],
            [0.5, 0.17, z],
            [0.83, 0.17, z],
            [0.5, 0.5, z],
        ]
    
    elif matched == "C2N":
        lattice = Lattice.hexagonal(a, c)
        # C2N with six-membered C-N rings around holes
        species = ["C"] * 6 + ["N"] * 3
        z = 0.5
        coords = [
            [0.1, 0.0, z], [0.2, 0.1, z], [0.3, 0.0, z],
            [0.4, 0.1, z], [0.5, 0.0, z], [0.6, 0.1, z],
            [0.15, 0.15, z], [0.35, 0.15, z], [0.55, 0.15, z],
        ]
    
    elif matched == "biphenylene":
        lattice = Lattice.orthorhombic(a, b, c)
        # 4-6-8 ring network
        species = ["C"] * 6
        z = 0.5
        coords = [
            [0.0, 0.0, z], [0.25, 0.0, z], [0.5, 0.0, z],
            [0.0, 0.5, z], [0.25, 0.5, z], [0.5, 0.5, z],
        ]
    
    elif "graphyne" in matched or matched == "graphdiyne":
        lattice = Lattice.hexagonal(a, c)
        # Simplified graphyne with sp-sp2 pattern
        species = ["C"] * 6
        z = 0.5
        coords = [
            [0.0, 0.0, z], [0.17, 0.17, z], [0.33, 0.0, z],
            [0.5, 0.17, z], [0.67, 0.0, z], [0.83, 0.17, z],
        ]
    
    elif matched == "holey_graphene":
        # Graphene with periodic vacancies
        lattice = Lattice.hexagonal(a, c)
        species = ["C"] * 5  # One vacancy per unit cell
        z = 0.5
        coords = [
            [0.0, 0.0, z], [0.25, 0.0, z], [0.5, 0.0, z],
            [0.125, 0.25, z], [0.375, 0.25, z],
        ]
    
    else:
        return {"success": False, "error": {"code": "NOT_IMPLEMENTED", "message": f"Structure not implemented: {matched}"}}
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": matched,
        "formula": params["formula"],
        "a_angstrom": a,
        "pore_size_angstrom": params["pore_size"],
        "description": params["description"],
        "applications": params.get("applications", []),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_holey_graphene(
    hole_size: int = 6,  # C6 hole
    hole_spacing: int = 4,
    base_size: List[int] = [6, 6, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate holey graphene with custom hole pattern.
    
    Args:
        hole_size: Atoms to remove per hole (6, 12, 24, ...)
        hole_spacing: Spacing between holes in unit cells
        base_size: Base graphene supercell
        vacuum: Vacuum padding
    
    Returns:
        Holey graphene structure
    """
    # Start with graphene
    a = 2.46
    c = vacuum
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    z = 0.5
    
    # Generate graphene supercell
    for i in range(base_size[0]):
        for j in range(base_size[1]):
            x1 = (1/3 + i) / base_size[0]
            y1 = (2/3 + j) / base_size[1]
            x2 = (2/3 + i) / base_size[0]
            y2 = (1/3 + j) / base_size[1]
            
            # Check if in hole region
            hole_i = i // hole_spacing
            hole_j = j // hole_spacing
            in_hole = (i % hole_spacing == 0 and j % hole_spacing == 0)
            
            if not in_hole:
                species.extend(["C", "C"])
                coords.extend([[x1, y1, z], [x2, y2, z]])
    
    lattice = Lattice.hexagonal(a * base_size[0], c)
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": "holey_graphene",
        "hole_size": hole_size,
        "hole_spacing": hole_spacing,
        "n_atoms": len(structure),
        "porosity": 1 - len(structure) / (base_size[0] * base_size[1] * 2),
        "structure": structure_to_dict(structure, vacuum)
    }
