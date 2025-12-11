"""
two_d/mxenes.py - MXene 2D Materials

Generates MXene structures (M_{n+1}X_n T_x):
- M = Ti, V, Nb, Mo, Ta, Cr, Zr, Hf, Sc
- X = C, N
- T = O, OH, F, Cl, S, bare
- n = 1, 2, 3 (M2X, M3X2, M4X3)

Scientific basis:
- Etched MAX phases
- Surface termination chemistry
- Experimental and DFT lattice constants
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# MXene database
MXENE_DATABASE = {
    # Ti-based
    "Ti2C": {"metal": "Ti", "x": "C", "n": 1, "a": 3.06, "n_metal_layers": 2},
    "Ti3C2": {"metal": "Ti", "x": "C", "n": 2, "a": 3.07, "n_metal_layers": 3},
    "Ti4C3": {"metal": "Ti", "x": "C", "n": 3, "a": 3.08, "n_metal_layers": 4},
    "Ti2N": {"metal": "Ti", "x": "N", "n": 1, "a": 3.00, "n_metal_layers": 2},
    "Ti4N3": {"metal": "Ti", "x": "N", "n": 3, "a": 3.04, "n_metal_layers": 4},
    # V-based
    "V2C": {"metal": "V", "x": "C", "n": 1, "a": 2.91, "n_metal_layers": 2},
    "V4C3": {"metal": "V", "x": "C", "n": 3, "a": 2.93, "n_metal_layers": 4},
    # Nb-based
    "Nb2C": {"metal": "Nb", "x": "C", "n": 1, "a": 3.13, "n_metal_layers": 2},
    "Nb4C3": {"metal": "Nb", "x": "C", "n": 3, "a": 3.16, "n_metal_layers": 4},
    # Mo-based
    "Mo2C": {"metal": "Mo", "x": "C", "n": 1, "a": 2.86, "n_metal_layers": 2},
    # Cr-based
    "Cr2C": {"metal": "Cr", "x": "C", "n": 1, "a": 2.85, "n_metal_layers": 2},
    # Sc-based
    "Sc2C": {"metal": "Sc", "x": "C", "n": 1, "a": 3.37, "n_metal_layers": 2},
    # Double (ordered)
    "Mo2TiC2": {"metals": ["Mo", "Ti", "Mo"], "x": "C", "n": 2, "a": 3.00, "n_metal_layers": 3, "ordered": True},
    "Cr2TiC2": {"metals": ["Cr", "Ti", "Cr"], "x": "C", "n": 2, "a": 3.05, "n_metal_layers": 3, "ordered": True},
}

# Termination parameters
TERMINATIONS = {
    "O": {"element": "O", "z_offset": 1.0, "description": "Oxygen termination"},
    "OH": {"element": "O", "z_offset": 1.0, "h_offset": 0.97, "description": "Hydroxyl termination"},
    "F": {"element": "F", "z_offset": 1.1, "description": "Fluorine termination"},
    "Cl": {"element": "Cl", "z_offset": 1.5, "description": "Chlorine termination"},
    "S": {"element": "S", "z_offset": 1.3, "description": "Sulfur termination"},
    "bare": {"element": None, "z_offset": 0, "description": "No termination"},
}


def generate_mxene(
    formula: str,
    termination: str = "O",
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate MXene structure.
    
    Args:
        formula: MXene formula (Ti3C2, V2C, Mo2TiC2, etc.)
        termination: Surface termination (O, OH, F, Cl, S, bare)
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        MXene structure dictionary
    
    Examples:
        >>> result = generate_mxene('Ti3C2', 'O')
        >>> result["formula_with_termination"]
        'Ti3C2O2'
    """
    if formula not in MXENE_DATABASE:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MXENE",
                "message": f"Unknown MXene '{formula}'",
                "available": list(MXENE_DATABASE.keys())
            }
        }
    
    if termination not in TERMINATIONS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_TERMINATION",
                "message": f"Unknown termination '{termination}'",
                "available": list(TERMINATIONS.keys())
            }
        }
    
    params = MXENE_DATABASE[formula]
    term_info = TERMINATIONS[termination]
    
    a = params["a"]
    n_metal = params["n_metal_layers"]
    x_element = params["x"]
    
    # Calculate MXene thickness
    layer_spacing = 2.3  # Å between metal layers
    total_thickness = (n_metal - 1) * layer_spacing + 2 * term_info["z_offset"] if termination != "bare" else (n_metal - 1) * layer_spacing
    
    c = vacuum
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    z_center = 0.5
    z_spacing = layer_spacing / c
    
    # Get metal elements
    if "ordered" in params and params["ordered"]:
        metals = params["metals"]
    else:
        metals = [params["metal"]] * n_metal
    
    # Metal layers
    for i in range(n_metal):
        offset = i - (n_metal - 1) / 2
        z = z_center + offset * z_spacing
        
        # Alternate ABC stacking
        if i % 3 == 0:
            xy = [1/3, 2/3]
        elif i % 3 == 1:
            xy = [2/3, 1/3]
        else:
            xy = [0, 0]
        
        species.append(metals[i])
        coords.append([xy[0], xy[1], z])
    
    # X-element (C or N) layers between metals
    for i in range(n_metal - 1):
        offset = i - (n_metal - 2) / 2
        z = z_center + offset * z_spacing
        species.append(x_element)
        coords.append([0, 0, z])
    
    # Terminations
    if termination != "bare":
        z_top = z_center + ((n_metal - 1) / 2 + 0.4) * z_spacing
        z_bot = z_center - ((n_metal - 1) / 2 + 0.4) * z_spacing
        
        term_elem = term_info["element"]
        species.extend([term_elem, term_elem])
        coords.extend([[1/3, 2/3, z_top], [2/3, 1/3, z_bot]])
        
        # Add H for OH termination
        if termination == "OH":
            h_z_top = z_top + 0.04
            h_z_bot = z_bot - 0.04
            species.extend(["H", "H"])
            coords.extend([[1/3, 2/3, h_z_top], [2/3, 1/3, h_z_bot]])
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Build formula with termination
    if termination == "bare":
        full_formula = formula
    elif termination == "OH":
        full_formula = f"{formula}(OH)2"
    else:
        full_formula = f"{formula}{term_info['element']}2"
    
    return {
        "success": True,
        "base_formula": formula,
        "termination": termination,
        "formula_with_termination": full_formula,
        "n_metal_layers": n_metal,
        "a_angstrom": a,
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_mxene_multilayer(
    formula: str,
    termination: str = "O",
    n_layers: int = 2,
    interlayer_spacing: float = 9.0,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate multi-layer MXene structure.
    
    Args:
        formula: MXene formula
        termination: Surface termination
        n_layers: Number of stacked MXene layers
        interlayer_spacing: Spacing between layers
        vacuum: Vacuum padding
    
    Returns:
        Multilayer MXene structure
    """
    # First generate single layer
    single = generate_mxene(formula, termination, vacuum=0)
    
    if not single["success"]:
        return single
    
    # Stack layers
    single_struct = single["structure"]
    atoms = single_struct["atoms"]
    lattice = single_struct["lattice"]
    
    a = lattice["a"]
    single_thickness = 6.0  # Approximate
    total_c = n_layers * (single_thickness + interlayer_spacing) - interlayer_spacing + vacuum
    
    new_lattice = Lattice.hexagonal(a, total_c)
    
    all_species = []
    all_coords = []
    
    for layer_idx in range(n_layers):
        z_offset = layer_idx * (single_thickness + interlayer_spacing) / total_c
        
        for atom in atoms:
            all_species.append(atom["element"])
            orig_z = atom["coords"][2]
            new_z = (orig_z * single_thickness / total_c) + z_offset
            all_coords.append([atom["coords"][0], atom["coords"][1], new_z])
    
    structure = Structure(new_lattice, all_species, all_coords)
    
    return {
        "success": True,
        "formula": formula,
        "termination": termination,
        "n_layers": n_layers,
        "interlayer_spacing": interlayer_spacing,
        "structure": structure_to_dict(structure, vacuum)
    }
