"""
two_d/janus.py - Janus 2D Materials (Broken Inversion)

Generates Janus monolayer structures with different top/bottom surfaces:
- MoSSe, MoSeS, MoSTe, MoTeSe
- WSSe, WSeS, WSTe, WTeSe
- PtSSe, SnSSe

Scientific basis:
- Broken inversion symmetry → out-of-plane dipole
- Piezoelectric response
- Rashba spin-orbit coupling
- DFT-optimized structural parameters
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Janus materials database
JANUS_MATERIALS = {
    # Mo-based Janus
    "MoSSe": {
        "metal": "Mo",
        "top": "S", "bottom": "Se",
        "a": 3.22, "d_MS": 2.42, "d_MSe": 2.54,
        "dipole": 0.02,  # Debye/Å²
        "description": "Most studied Janus TMD"
    },
    "MoSeS": {
        "metal": "Mo",
        "top": "Se", "bottom": "S",
        "a": 3.22, "d_MS": 2.42, "d_MSe": 2.54,
        "dipole": -0.02,
        "description": "MoSSe with flipped surfaces"
    },
    "MoSTe": {
        "metal": "Mo",
        "top": "S", "bottom": "Te",
        "a": 3.42, "d_MS": 2.42, "d_MTe": 2.73,
        "dipole": 0.08,
        "description": "Large asymmetry"
    },
    # W-based Janus
    "WSSe": {
        "metal": "W",
        "top": "S", "bottom": "Se",
        "a": 3.25, "d_WS": 2.42, "d_WSe": 2.55,
        "dipole": 0.02,
        "description": "W-based Janus TMD"
    },
    "WSeS": {
        "metal": "W",
        "top": "Se", "bottom": "S",
        "a": 3.25, "d_WS": 2.42, "d_WSe": 2.55,
        "dipole": -0.02,
        "description": "WSSe with flipped surfaces"
    },
    "WSTe": {
        "metal": "W",
        "top": "S", "bottom": "Te",
        "a": 3.43, "d_WS": 2.42, "d_WTe": 2.74,
        "dipole": 0.09,
        "description": "Large asymmetry W-based"
    },
    # Other Janus
    "PtSSe": {
        "metal": "Pt",
        "top": "S", "bottom": "Se",
        "a": 3.60, "d_PtS": 2.38, "d_PtSe": 2.50,
        "dipole": 0.04,
        "description": "Pt-based Janus"
    },
    "SnSSe": {
        "metal": "Sn",
        "top": "S", "bottom": "Se",
        "a": 3.95, "d_SnS": 2.60, "d_SnSe": 2.75,
        "dipole": 0.03,
        "description": "Group-IV Janus"
    },
}


def generate_janus_2d(
    formula: str = None,
    metal: str = None,
    top_chalcogen: str = None,
    bottom_chalcogen: str = None,
    a: Optional[float] = None,
    vacuum: float = 15.0,
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate Janus 2D structure.
    
    Can be called with either:
    - Named formula: generate_janus_2d("MoSSe")
    - Components: generate_janus_2d(metal="Mo", top="S", bottom="Se")
    
    Args:
        formula: Janus formula (MoSSe, WSSe, etc.)
        metal: Transition metal element
        top_chalcogen: Top chalcogen element
        bottom_chalcogen: Bottom chalcogen element
        a: Lattice constant (estimated if not provided)
        vacuum: Vacuum padding
        size: Supercell size
    
    Returns:
        Janus 2D structure dictionary
    
    Examples:
        >>> result = generate_janus_2d("MoSSe")
        >>> result["out_of_plane_dipole"]
        True
    """
    # Parse inputs
    if formula and formula in JANUS_MATERIALS:
        params = JANUS_MATERIALS[formula]
        metal = params["metal"]
        top = params["top"]
        bottom = params["bottom"]
        a_val = a or params["a"]
        d_top = params.get(f"d_{metal}{top}", 2.5)
        d_bottom = params.get(f"d_{metal}{bottom}", 2.5)
        dipole = params.get("dipole", 0.0)
        description = params.get("description", "")
    elif metal and top_chalcogen and bottom_chalcogen:
        top = top_chalcogen
        bottom = bottom_chalcogen
        a_val = a or 3.25
        d_top = 2.45
        d_bottom = 2.55
        dipole = 0.02 if top != bottom else 0.0
        formula = f"{metal}{top}{bottom}"
        description = f"Custom Janus {formula}"
    else:
        return {
            "success": False,
            "error": {
                "code": "INVALID_INPUT",
                "message": "Provide either 'formula' or 'metal+top_chalcogen+bottom_chalcogen'",
                "available": list(JANUS_MATERIALS.keys())
            }
        }
    
    c = vacuum
    lattice = Lattice.hexagonal(a_val, c)
    
    z_center = 0.5
    dz_top = d_top / c
    dz_bottom = d_bottom / c
    
    species = [metal, top, bottom]
    coords = [
        [1/3, 2/3, z_center],               # Metal
        [1/3, 2/3, z_center + dz_top],     # Top chalcogen
        [1/3, 2/3, z_center - dz_bottom],  # Bottom chalcogen
    ]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Calculate asymmetry metrics
    thickness = d_top + d_bottom
    asymmetry = abs(d_top - d_bottom) / thickness if thickness > 0 else 0
    
    return {
        "success": True,
        "formula": formula,
        "metal": metal,
        "top_chalcogen": top,
        "bottom_chalcogen": bottom,
        "a_angstrom": a_val,
        "thickness_angstrom": thickness,
        "asymmetry": round(asymmetry, 3),
        "out_of_plane_dipole": True,
        "dipole_debye_per_A2": dipole,
        "description": description,
        "properties": {
            "broken_inversion": True,
            "piezoelectric": True,
            "rashba_soc": True
        },
        "structure": structure_to_dict(structure, vacuum)
    }


def list_janus_materials() -> Dict[str, Any]:
    """
    List all available Janus materials.
    
    Returns:
        Dictionary with Janus material information
    """
    materials = []
    for name, info in JANUS_MATERIALS.items():
        materials.append({
            "formula": name,
            "metal": info["metal"],
            "top": info["top"],
            "bottom": info["bottom"],
            "a": info["a"],
            "description": info.get("description", "")
        })
    
    return {
        "success": True,
        "count": len(materials),
        "materials": materials
    }
