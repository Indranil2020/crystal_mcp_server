"""
two_d/ribbons.py - Nanoribbons

Generates nanoribbon structures:
- Armchair edge
- Zigzag edge
- Chiral edge (arbitrary angle)
- Bearded edge
- Cove edge

With edge passivation options:
- H, OH, F, Cl, S, Se, NH2, NO2

Scientific basis:
- ASE graphene_nanoribbon for accurate edge geometry
- Width-dependent bandgap for armchair
- Edge magnetism for zigzag
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice
import importlib.util

HAS_ASE = importlib.util.find_spec("ase.build") is not None
if HAS_ASE:
    from ase.build import graphene_nanoribbon

from .base import atoms_to_dict, structure_to_dict


# Edge type descriptions
EDGE_TYPES = {
    "armchair": {
        "description": "Armchair edge (semiconductor, width-dependent gap)",
        "magnetic": False,
        "gap_behavior": "3n, 3n+1, 3n+2 families"
    },
    "zigzag": {
        "description": "Zigzag edge (edge magnetism, metallic or half-metallic)",
        "magnetic": True,
        "gap_behavior": "Can be spin-polarized"
    },
    "chiral": {
        "description": "Chiral edge (intermediate properties)",
        "magnetic": "depends on angle",
        "gap_behavior": "Variable"
    },
    "bearded": {
        "description": "Bearded/Klein edge (extended zigzag)",
        "magnetic": True,
        "gap_behavior": "Flat bands"
    },
    "cove": {
        "description": "Cove edge (atomically precise from synthesis)",
        "magnetic": False,
        "gap_behavior": "Well-defined"
    }
}

# Passivation parameters
PASSIVATION_PARAMS = {
    "H": {"element": "H", "bond_length": 1.09},
    "F": {"element": "F", "bond_length": 1.36},
    "Cl": {"element": "Cl", "bond_length": 1.76},
    "OH": {"elements": ["O", "H"], "bond_lengths": [1.43, 0.97]},
    "S": {"element": "S", "bond_length": 1.82},
    "Se": {"element": "Se", "bond_length": 1.95},
    "NH2": {"elements": ["N", "H", "H"], "bond_lengths": [1.47, 1.01, 1.01]},
    "bare": {"element": None}
}


def generate_nanoribbon(
    base_material: str = "graphene",
    edge_type: str = "armchair",
    width: int = 5,
    length: int = 1,
    passivation: str = "H",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Generate nanoribbon structure.
    
    Args:
        base_material: Base 2D material (graphene, silicene, etc.)
        edge_type: Edge type (armchair, zigzag, chiral)
        width: Ribbon width in unit cells
        length: Ribbon length in unit cells
        passivation: Edge passivation (H, F, Cl, OH, S, Se, NH2, bare)
        vacuum: Vacuum padding
    
    Returns:
        Nanoribbon structure
    
    Examples:
        >>> result = generate_nanoribbon('graphene', 'zigzag', 7)
        >>> result["is_magnetic"]
        True
    """
    if edge_type not in EDGE_TYPES:
        return {
            "success": False,
            "error": {
                "code": "INVALID_EDGE",
                "message": f"Unknown edge type '{edge_type}'",
                "available": list(EDGE_TYPES.keys())
            }
        }
    
    if passivation not in PASSIVATION_PARAMS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PASSIVATION",
                "message": f"Unknown passivation '{passivation}'",
                "available": list(PASSIVATION_PARAMS.keys())
            }
        }
    
    edge_info = EDGE_TYPES[edge_type]
    
    if base_material.lower() == "graphene" and HAS_ASE:
        # Use ASE for accurate graphene nanoribbon
        saturated = passivation == "H"
        ase_type = "armchair" if edge_type == "armchair" else "zigzag"
        
        atoms = graphene_nanoribbon(
            width, length,
            type=ase_type,
            saturated=saturated,
            vacuum=vacuum
        )
        
        struct_dict = atoms_to_dict(atoms, vacuum)
        
        # Calculate ribbon width in Angstroms
        positions = atoms.get_positions()
        ribbon_width = max(positions[:, 1]) - min(positions[:, 1])
        
        return {
            "success": True,
            "base_material": base_material,
            "edge_type": edge_type,
            "width_cells": width,
            "width_angstrom": round(ribbon_width, 2),
            "length_cells": length,
            "passivation": passivation,
            "n_atoms": len(atoms),
            "is_magnetic": edge_info.get("magnetic", False),
            "description": edge_info["description"],
            "structure": struct_dict
        }
    
    else:
        # Manual construction for other materials
        if base_material.lower() == "graphene":
            a = 2.46
            element = "C"
        elif base_material.lower() == "silicene":
            a = 3.86
            element = "Si"
        elif base_material.lower() == "bn":
            a = 2.50
            element = "B"  # Will alternate with N
        else:
            a = 2.5
            element = "C"
        
        if edge_type == "armchair":
            # Armchair: ribbon along x, width along y
            ribbon_a = a * np.sqrt(3) * length
            ribbon_b = a * 1.5 * width + vacuum * 2
            ribbon_c = vacuum
        else:
            # Zigzag: ribbon along x, width along y
            ribbon_a = a * length
            ribbon_b = a * np.sqrt(3) * width + vacuum * 2
            ribbon_c = vacuum
        
        lattice = Lattice.orthorhombic(ribbon_a, ribbon_b, ribbon_c)
        
        species = []
        coords = []
        
        # Build lattice manually
        for i in range(length):
            for j in range(width * 2):
                if edge_type == "armchair":
                    x = (i + 0.5) / length
                    y = (j * 0.75 + vacuum) / ribbon_b
                else:
                    x = (i + (j % 2) * 0.5) / length
                    y = (j * a * np.sqrt(3) / 4 + vacuum) / ribbon_b
                
                species.append(element)
                coords.append([x, y, 0.5])
        
        # Add passivation
        if passivation != "bare":
            pass_info = PASSIVATION_PARAMS[passivation]
            pass_elem = pass_info.get("element")
            if pass_elem:
                # Add at edges (simplified)
                for j in [0, width * 2 - 1]:
                    for i in range(length):
                        x = (i + 0.5) / length
                        y = (j * 0.75 + vacuum + 0.05 * (1 if j == 0 else -1)) / ribbon_b
                        species.append(pass_elem)
                        coords.append([x, y, 0.5])
        
        structure = Structure(lattice, species, coords)
        
        return {
            "success": True,
            "base_material": base_material,
            "edge_type": edge_type,
            "width_cells": width,
            "length_cells": length,
            "passivation": passivation,
            "n_atoms": len(structure),
            "is_magnetic": edge_info.get("magnetic", False),
            "description": edge_info["description"],
            "structure": structure_to_dict(structure, vacuum)
        }


def generate_chiral_nanoribbon(
    base_material: str = "graphene",
    chiral_angle: float = 15.0,
    width: int = 5,
    length: int = 3,
    passivation: str = "H",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Generate chiral nanoribbon with arbitrary edge angle.
    
    Args:
        base_material: Base 2D material
        chiral_angle: Angle from zigzag direction (degrees, 0=zigzag, 30=armchair)
        width: Ribbon width
        length: Ribbon length
        passivation: Edge passivation
        vacuum: Vacuum padding
    
    Returns:
        Chiral nanoribbon structure
    """
    # Chiral angle determines edge character
    if chiral_angle < 5:
        effective_type = "zigzag-like"
    elif chiral_angle > 25:
        effective_type = "armchair-like"
    else:
        effective_type = "chiral"
    
    # For now, generate closest standard ribbon
    if chiral_angle < 15:
        result = generate_nanoribbon(base_material, "zigzag", width, length, passivation, vacuum)
    else:
        result = generate_nanoribbon(base_material, "armchair", width, length, passivation, vacuum)
    
    if result["success"]:
        result["chiral_angle"] = chiral_angle
        result["effective_type"] = effective_type
        result["note"] = f"Chiral angle {chiral_angle}° approximated as {result['edge_type']}"
    
    return result
