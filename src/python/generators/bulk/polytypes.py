"""
bulk/polytypes.py - Polytype Stacking Sequences

Generates polytypes with different stacking sequences (2H, 4H, 6H, 15R, etc.)
for materials like SiC, ZnS, GaN, and close-packed metals.

Scientific basis:
- ABC stacking notation for close-packed layers
- Ramsdell notation (nX where n=layers, X=symmetry)
- Jagodzinski notation (h/c for hexagonal/cubic stacking)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice
from ase import Atoms

from .base import structure_to_dict, atoms_to_dict


# Standard polytype stacking sequences
POLYTYPE_SEQUENCES = {
    # Wurtzite-type (2H)
    "2H": ["A", "B"],
    # Zincblende-type (3C)
    "3C": ["A", "B", "C"],
    # Common SiC polytypes
    "4H": ["A", "B", "C", "B"],
    "6H": ["A", "B", "C", "A", "C", "B"],
    "15R": ["A", "B", "C", "A", "C", "B", "C", "A", "B", "A", "C", "B", "A", "B", "C"],
    # Additional polytypes
    "8H": ["A", "B", "A", "C", "A", "B", "A", "C"],
    "10H": ["A", "B", "C", "B", "A", "B", "C", "B", "A", "B"],
    "21R": list("ABCACBCABACACBABCBCAB"),
    # Simple metal stacking
    "fcc": ["A", "B", "C"],
    "hcp": ["A", "B"],
    "dhcp": ["A", "B", "A", "C"],  # Double hexagonal close-packed
}

# Layer-to-fractional coordinate mapping
LAYER_POSITIONS = {
    "A": (0.0, 0.0),
    "B": (1/3, 2/3),
    "C": (2/3, 1/3)
}


def generate_polytype(
    element: str,
    polytype: str,
    a: Optional[float] = None,
    c_layer: Optional[float] = None,
    bilayer_element: Optional[str] = None,
    vacuum: float = 0.0
) -> Dict[str, Any]:
    """
    Generate a polytype structure with specified stacking sequence.
    
    Args:
        element: Primary element symbol (e.g., 'Si' for SiC)
        polytype: Polytype designation (2H, 3C, 4H, 6H, 15R, etc.)
        a: In-plane lattice constant (Å) - estimated if not provided
        c_layer: Layer spacing (Å) - estimated if not provided
        bilayer_element: Second element for binary compounds (e.g., 'C' for SiC)
        vacuum: Vacuum to add (for slab calculations)
    
    Returns:
        Dictionary with polytype structure
    
    Examples:
        >>> result = generate_polytype('Si', '6H', a=3.08, c_layer=2.52, bilayer_element='C')
        >>> result["success"]
        True
        >>> result["n_layers"]
        6
    """
    if polytype not in POLYTYPE_SEQUENCES:
        return {
            "success": False,
            "error": {
                "code": "INVALID_POLYTYPE",
                "message": f"Unknown polytype '{polytype}'",
                "available": list(POLYTYPE_SEQUENCES.keys())
            }
        }
    
    sequence = POLYTYPE_SEQUENCES[polytype]
    n_layers = len(sequence)
    
    # Estimate lattice parameters if not provided
    # Based on tabulated values for common materials
    LATTICE_DATABASE = {
        "Si": {"a": 3.84, "c_layer": 2.35},
        "C": {"a": 2.52, "c_layer": 2.10},
        "Zn": {"a": 2.66, "c_layer": 2.12},
        "Ga": {"a": 3.19, "c_layer": 2.60},
        "Cd": {"a": 2.98, "c_layer": 2.81},
        "Co": {"a": 2.51, "c_layer": 2.02},
    }
    
    # For bilayer compounds like SiC
    BILAYER_DATABASE = {
        ("Si", "C"): {"a": 3.08, "c_layer": 2.52},
        ("Ga", "N"): {"a": 3.19, "c_layer": 2.59},
        ("Zn", "O"): {"a": 3.25, "c_layer": 2.60},
        ("Cd", "S"): {"a": 4.14, "c_layer": 3.37},
        ("Zn", "S"): {"a": 3.82, "c_layer": 3.12},
    }
    
    if a is None or c_layer is None:
        if bilayer_element:
            key = (element, bilayer_element)
            if key in BILAYER_DATABASE:
                a = a or BILAYER_DATABASE[key]["a"]
                c_layer = c_layer or BILAYER_DATABASE[key]["c_layer"]
            else:
                a = a or 3.0
                c_layer = c_layer or 2.5
        else:
            if element in LATTICE_DATABASE:
                a = a or LATTICE_DATABASE[element]["a"]
                c_layer = c_layer or LATTICE_DATABASE[element]["c_layer"]
            else:
                a = a or 3.0
                c_layer = c_layer or 2.5
    
    # Calculate total c parameter
    c = c_layer * n_layers + vacuum
    
    # Build hexagonal lattice
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    for i, layer_type in enumerate(sequence):
        xy = LAYER_POSITIONS[layer_type]
        z = (i + 0.25) / n_layers * (1 - vacuum / c)  # Normalize to unit cell
        
        # Add primary atom
        species.append(element)
        coords.append([xy[0], xy[1], z])
        
        # Add bilayer atom if specified
        if bilayer_element:
            z_offset = 0.125 / n_layers  # Typically 1/4 of layer spacing
            species.append(bilayer_element)
            coords.append([xy[0], xy[1], z + z_offset])
    
    structure = Structure(lattice, species, coords)
    
    # Determine crystal system
    crystal_system = "hexagonal" if polytype.endswith("H") else "rhombohedral" if polytype.endswith("R") else "cubic"
    
    return {
        "success": True,
        "polytype": polytype,
        "sequence": sequence,
        "n_layers": n_layers,
        "crystal_system": crystal_system,
        "lattice_params": {"a": a, "c": c, "c_layer": c_layer},
        "structure": structure_to_dict(structure)
    }


def generate_custom_stacking(
    element: str,
    sequence: List[str],
    a: float = 3.0,
    c_layer: float = 2.5,
    bilayer_element: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a structure with custom stacking sequence.
    
    Args:
        element: Primary element
        sequence: Custom stacking sequence, e.g., ["A", "B", "A", "C", "A", "B"]
        a: In-plane lattice constant
        c_layer: Layer spacing
        bilayer_element: Optional second element
    
    Returns:
        Structure with custom stacking
    """
    # Validate sequence
    valid_layers = set(LAYER_POSITIONS.keys())
    for layer in sequence:
        if layer.upper() not in valid_layers:
            return {"success": False, "error": {"code": "INVALID_LAYER", 
                    "message": f"Layer must be A, B, or C, got '{layer}'"}}
    
    n_layers = len(sequence)
    c = c_layer * n_layers
    
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    for i, layer_type in enumerate(sequence):
        xy = LAYER_POSITIONS[layer_type.upper()]
        z = (i + 0.25) / n_layers
        
        species.append(element)
        coords.append([xy[0], xy[1], z])
        
        if bilayer_element:
            species.append(bilayer_element)
            coords.append([xy[0], xy[1], z + 0.125 / n_layers])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "sequence": sequence,
        "n_layers": n_layers,
        "structure": structure_to_dict(structure)
    }
