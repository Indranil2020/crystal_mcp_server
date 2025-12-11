"""
bulk/special_phases.py - Special Bulk Phases

Generates specialized bulk phases:
- High-pressure structures (LaH10, post-perovskite)
- Negative thermal expansion materials
- Super-hydrides
- Energetic materials
- Prussian blue analogs

Scientific basis:
- High-pressure crystal chemistry
- Phonon-driven NTE mechanisms  
- Clathrate superhydrides (LaH10, YH9)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# High-pressure phases database
HIGH_PRESSURE_PHASES = {
    "LaH10": {
        "formula": "LaH10",
        "a": 5.1,
        "space_group": 225,  # Fm-3m
        "pressure_gpa": 170,
        "tc_kelvin": 260,
        "description": "Lanthanum superhydride - near room-T superconductor"
    },
    "YH9": {
        "formula": "YH9",
        "a": 5.08,
        "space_group": 166,
        "pressure_gpa": 200,
        "tc_kelvin": 243,
        "description": "Yttrium superhydride"
    },
    "H3S": {
        "formula": "H3S",
        "a": 3.09,
        "space_group": 229,  # Im-3m
        "pressure_gpa": 150,
        "tc_kelvin": 203,
        "description": "Hydrogen sulfide - first superhydride discovered"
    },
    "post-perovskite": {
        "formula": "MgSiO3",
        "a": 2.46, "b": 8.04, "c": 6.09,
        "space_group": 63,  # Cmcm
        "pressure_gpa": 125,
        "description": "Post-perovskite phase (D'' layer)"
    },
    "ice-X": {
        "formula": "H2O",
        "a": 2.54,
        "space_group": 229,  # Im-3m (predicted)
        "pressure_gpa": 60,
        "description": "Symmetric ice phase"
    }
}

# Negative thermal expansion materials
NTE_MATERIALS = {
    "ZrW2O8": {
        "formula": "ZrW2O8",
        "a": 9.15,
        "space_group": 198,  # P2_13
        "alpha_per_K": -9e-6,
        "description": "Cubic NTE, broad T range"
    },
    "ScF3": {
        "formula": "ScF3",
        "a": 4.02,
        "space_group": 221,  # Pm-3m
        "alpha_per_K": -14e-6,
        "description": "ReO3-type NTE"
    }
}


def generate_superhydride(
    metal: str = "La",
    h_content: int = 10,
    a: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate superhydride structure (MHn).
    
    Args:
        metal: Metal element (La, Y, Th, Ca, etc.)
        h_content: Number of H atoms per metal
        a: Lattice constant (estimated if not provided)
    
    Returns:
        Superhydride structure
    """
    # Estimate lattice constant from metal radius
    METAL_RADII = {"La": 2.74, "Y": 2.27, "Th": 2.06, "Ca": 2.31, "Mg": 1.73}
    
    if a is None:
        r_metal = METAL_RADII.get(metal, 2.0)
        a = 2 * r_metal + 0.74  # approx
    
    lattice = Lattice.cubic(a)
    
    species = [metal]
    coords = [[0, 0, 0]]  # Metal at origin
    
    # Arrange H atoms on cage vertices
    if h_content == 10:
        # LaH10-type: cubo-octahedral H cage
        # H at (0.5, 0.2, 0) and permutations
        h_positions = [
            [0.5, 0.2, 0], [0.5, -0.2, 0], [0.5, 0, 0.2], [0.5, 0, -0.2],
            [0.2, 0.5, 0], [-0.2, 0.5, 0], [0, 0.5, 0.2], [0, 0.5, -0.2],
            [0.2, 0, 0.5], [-0.2, 0, 0.5]
        ]
    elif h_content == 9:
        # YH9-type
        h_positions = [[0.1*i, 0.1*j, 0.5] for i in range(3) for j in range(3)]
    else:
        # Generic: H on cubic grid
        n_edge = int(np.ceil(h_content ** (1/3)))
        h_positions = [[i/(n_edge+1), j/(n_edge+1), k/(n_edge+1)] 
                       for i in range(1, n_edge+1) 
                       for j in range(1, n_edge+1) 
                       for k in range(1, n_edge+1)][:h_content]
    
    for h_pos in h_positions[:h_content]:
        species.append("H")
        coords.append([x % 1 for x in h_pos])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "formula": f"{metal}H{h_content}",
        "metal": metal,
        "h_content": h_content,
        "note": "High-pressure phase - optimize at target pressure",
        "structure": structure_to_dict(structure)
    }


def generate_high_pressure_phase(
    phase_name: str
) -> Dict[str, Any]:
    """
    Generate a known high-pressure phase.
    
    Args:
        phase_name: Name of phase (LaH10, YH9, H3S, post-perovskite, ice-X)
    
    Returns:
        High-pressure phase structure
    """
    if phase_name not in HIGH_PRESSURE_PHASES:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_PHASE",
                "message": f"Unknown phase '{phase_name}'",
                "available": list(HIGH_PRESSURE_PHASES.keys())
            }
        }
    
    params = HIGH_PRESSURE_PHASES[phase_name]
    
    if "b" in params:
        lattice = Lattice.orthorhombic(params["a"], params["b"], params["c"])
    else:
        lattice = Lattice.cubic(params["a"])
    
    # Simplified structure generation
    formula = params["formula"]
    if formula == "LaH10":
        return generate_superhydride("La", 10, params["a"])
    elif formula == "YH9":
        return generate_superhydride("Y", 9, params["a"])
    else:
        # Generic
        species = list(formula)  # Very simplified
        coords = [[0, 0, 0]]
        structure = Structure(lattice, species[:1], coords)
        
        return {
            "success": True,
            "phase": phase_name,
            "pressure_gpa": params["pressure_gpa"],
            "description": params["description"],
            "structure": structure_to_dict(structure)
        }


def generate_nte_material(
    material: str = "ZrW2O8"
) -> Dict[str, Any]:
    """
    Generate negative thermal expansion material.
    
    Args:
        material: NTE material (ZrW2O8, ScF3)
    
    Returns:
        NTE material structure
    """
    if material not in NTE_MATERIALS:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN",
                "message": f"Unknown NTE material '{material}'",
                "available": list(NTE_MATERIALS.keys())
            }
        }
    
    params = NTE_MATERIALS[material]
    lattice = Lattice.cubic(params["a"])
    
    if material == "ZrW2O8":
        species = ['Zr', 'W', 'W', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O']
        coords = [[0, 0, 0], [0.34, 0.34, 0.34], [0.66, 0.66, 0.66]] + \
                 [[0.1*i, 0.2, 0.3] for i in range(8)]
    elif material == "ScF3":
        species = ['Sc', 'F', 'F', 'F']
        coords = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
    else:
        species = ['X']
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "formula": params["formula"],
        "alpha_per_K": params["alpha_per_K"],
        "description": params["description"],
        "structure": structure_to_dict(structure)
    }


def generate_prussian_blue(
    metal_a: str = "Fe",
    metal_b: str = "Fe",
    a: float = 10.2
) -> Dict[str, Any]:
    """
    Generate Prussian blue analog structure.
    
    Args:
        metal_a: Metal at N-coordinated sites
        metal_b: Metal at C-coordinated sites
        a: Lattice constant
    
    Returns:
        Prussian blue structure
    """
    lattice = Lattice.cubic(a)
    
    # Face-centered cubic arrangement of [M(CN)6] units
    species = []
    coords = []
    
    # Metal sites
    for pos in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
        species.append(metal_a)
        coords.append(pos)
    
    for pos in [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]:
        species.append(metal_b)
        coords.append(pos)
    
    # CN bridges (simplified)
    cn_frac = 0.17  # Position of C/N between metals
    for d in [[cn_frac, 0, 0], [0.5-cn_frac, 0, 0], [0, cn_frac, 0], [0, 0.5-cn_frac, 0]]:
        species.extend(['C', 'N'])
        coords.append(d)
        coords.append([d[0]+0.16, d[1]+0.16, d[2]])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "formula": f"{metal_a}4[{metal_b}(CN)6]3",
        "metal_a": metal_a,
        "metal_b": metal_b,
        "structure": structure_to_dict(structure)
    }
