"""
two_d/tmds.py - Transition Metal Dichalcogenides

Generates TMD monolayer structures:
- MoS2, MoSe2, MoTe2
- WS2, WSe2, WTe2
- NbS2, NbSe2, TaS2, TaSe2
- Phases: 1H (2H), 1T, 1T', Td, 1T‴

Scientific basis:
- Trigonal prismatic (H) vs octahedral (T) coordination
- Phase transitions and charge density waves
- Experimental lattice constants
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# TMD database with experimental parameters
TMD_DATABASE = {
    # Mo-based
    "MoS2": {"metal": "Mo", "chalcogen": "S", "a": 3.16, "thickness": 3.13, "phases": ["1H", "1T", "1T'"]},
    "MoSe2": {"metal": "Mo", "chalcogen": "Se", "a": 3.29, "thickness": 3.34, "phases": ["1H", "1T"]},
    "MoTe2": {"metal": "Mo", "chalcogen": "Te", "a": 3.52, "thickness": 3.60, "phases": ["1H", "1T'", "Td"]},
    # W-based
    "WS2": {"metal": "W", "chalcogen": "S", "a": 3.18, "thickness": 3.14, "phases": ["1H", "1T'"]},
    "WSe2": {"metal": "W", "chalcogen": "Se", "a": 3.28, "thickness": 3.35, "phases": ["1H", "1T"]},
    "WTe2": {"metal": "W", "chalcogen": "Te", "a": 3.48, "thickness": 3.53, "phases": ["Td"]},  # Semimetal
    # Nb-based
    "NbS2": {"metal": "Nb", "chalcogen": "S", "a": 3.33, "thickness": 3.10, "phases": ["1H", "3R"]},
    "NbSe2": {"metal": "Nb", "chalcogen": "Se", "a": 3.44, "thickness": 3.30, "phases": ["1H", "2H"]},
    # Ta-based
    "TaS2": {"metal": "Ta", "chalcogen": "S", "a": 3.36, "thickness": 3.10, "phases": ["1H", "1T"]},
    "TaSe2": {"metal": "Ta", "chalcogen": "Se", "a": 3.48, "thickness": 3.30, "phases": ["1H", "1T"]},
    # Others
    "HfS2": {"metal": "Hf", "chalcogen": "S", "a": 3.63, "thickness": 2.91, "phases": ["1T"]},
    "HfSe2": {"metal": "Hf", "chalcogen": "Se", "a": 3.77, "thickness": 3.17, "phases": ["1T"]},
    "ZrS2": {"metal": "Zr", "chalcogen": "S", "a": 3.68, "thickness": 2.97, "phases": ["1T"]},
    "PtS2": {"metal": "Pt", "chalcogen": "S", "a": 3.54, "thickness": 2.53, "phases": ["1T"]},
    "PtSe2": {"metal": "Pt", "chalcogen": "Se", "a": 3.73, "thickness": 2.70, "phases": ["1T"]},
}

# Phase coordination descriptions
PHASE_INFO = {
    "1H": {"coordination": "trigonal_prismatic", "symmetry": "D3h", "metallic": False},
    "2H": {"coordination": "trigonal_prismatic", "symmetry": "D3h", "metallic": False, "layers": 2},
    "1T": {"coordination": "octahedral", "symmetry": "D3d", "metallic": True},
    "1T'": {"coordination": "distorted_octahedral", "symmetry": "C2h", "metallic": True},
    "Td": {"coordination": "twisted_octahedral", "symmetry": "C2v", "metallic": True, "weyl": True},
    "3R": {"coordination": "trigonal_prismatic", "symmetry": "R3m", "layers": 3},
}


def generate_tmd(
    formula: str,
    phase: str = "1H",
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate TMD monolayer structure.
    
    Args:
        formula: TMD formula (MoS2, WS2, MoSe2, etc.)
        phase: Crystal phase (1H, 1T, 1T', Td)
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        TMD structure dictionary
    
    Examples:
        >>> result = generate_tmd('MoS2', '1H')
        >>> result["coordination"]
        'trigonal_prismatic'
    """
    if formula not in TMD_DATABASE:
        return {
            "success": False,
            "error": {
                "code": "INVALID_TMD",
                "message": f"Unknown TMD '{formula}'",
                "available": list(TMD_DATABASE.keys())
            }
        }
    
    params = TMD_DATABASE[formula]
    available_phases = params["phases"]
    
    if phase not in available_phases:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PHASE",
                "message": f"Phase '{phase}' not available for {formula}",
                "available": available_phases
            }
        }
    
    metal = params["metal"]
    chalcogen = params["chalcogen"]
    a = params["a"]
    thickness = params["thickness"]
    
    c = vacuum
    dz = thickness / c / 2 if c > 0 else 0.1
    z_center = 0.5
    
    if phase in ["1H", "2H"]:
        # Trigonal prismatic: X-M-X sandwich, X atoms eclipsed
        lattice = Lattice.hexagonal(a, c)
        species = [metal, chalcogen, chalcogen]
        coords = [
            [1/3, 2/3, z_center],              # Metal
            [1/3, 2/3, z_center + dz],         # Top chalcogen
            [1/3, 2/3, z_center - dz],         # Bottom chalcogen
        ]
    
    elif phase == "1T":
        # Octahedral: X atoms staggered
        lattice = Lattice.hexagonal(a, c)
        species = [metal, chalcogen, chalcogen]
        coords = [
            [0, 0, z_center],                  # Metal
            [1/3, 2/3, z_center + dz],         # Top chalcogen
            [2/3, 1/3, z_center - dz],         # Bottom chalcogen (rotated)
        ]
    
    elif phase == "1T'":
        # Distorted octahedral: doubled unit cell with dimerization
        lattice = Lattice.orthorhombic(a, a * np.sqrt(3), c)
        species = [metal, metal, chalcogen, chalcogen, chalcogen, chalcogen]
        
        # Metal dimerization
        dm = 0.03  # Dimerization displacement
        coords = [
            [0.25 - dm, 0.25, z_center],
            [0.75 + dm, 0.25, z_center],
            [0.25, 0.0, z_center + dz], [0.75, 0.0, z_center + dz],
            [0.25, 0.5, z_center - dz], [0.75, 0.5, z_center - dz],
        ]
    
    elif phase == "Td":
        # Twisted 1T': broken inversion, Weyl semimetal
        lattice = Lattice.orthorhombic(a, a * np.sqrt(3), c)
        species = [metal, metal, chalcogen, chalcogen, chalcogen, chalcogen]
        
        # Asymmetric distortion
        dm = 0.03
        twist = 0.02
        coords = [
            [0.25 - dm, 0.25 + twist, z_center],
            [0.75 + dm, 0.25 - twist, z_center],
            [0.25, 0.0, z_center + dz], [0.75, 0.0, z_center + dz],
            [0.25, 0.5, z_center - dz], [0.75, 0.5, z_center - dz],
        ]
    
    else:
        # Default to 1H
        lattice = Lattice.hexagonal(a, c)
        species = [metal, chalcogen, chalcogen]
        coords = [[1/3, 2/3, z_center], [1/3, 2/3, z_center + dz], [1/3, 2/3, z_center - dz]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    phase_info = PHASE_INFO.get(phase, {})
    
    return {
        "success": True,
        "formula": formula,
        "phase": phase,
        "metal": metal,
        "chalcogen": chalcogen,
        "a_angstrom": a,
        "thickness_angstrom": thickness,
        "coordination": phase_info.get("coordination", "unknown"),
        "is_metallic": phase_info.get("metallic", False),
        "is_weyl": phase_info.get("weyl", False),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_tmd_heterostructure(
    tmd_1: str,
    tmd_2: str,
    alignment: str = "AA",
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate TMD vertical heterostructure.
    
    Args:
        tmd_1: Bottom TMD
        tmd_2: Top TMD
        alignment: Stacking (AA, AB)
        vacuum: Vacuum padding
    
    Returns:
        TMD heterostructure
    """
    if tmd_1 not in TMD_DATABASE or tmd_2 not in TMD_DATABASE:
        return {"success": False, "error": {"code": "INVALID_TMD", "message": "Unknown TMD"}}
    
    p1 = TMD_DATABASE[tmd_1]
    p2 = TMD_DATABASE[tmd_2]
    
    # Use average lattice constant
    a_avg = (p1["a"] + p2["a"]) / 2
    mismatch = abs(p1["a"] - p2["a"]) / a_avg * 100
    
    interlayer_sep = 6.5  # Å
    t1 = p1["thickness"]
    t2 = p2["thickness"]
    
    c = t1/2 + interlayer_sep + t2/2 + vacuum
    lattice = Lattice.hexagonal(a_avg, c)
    
    z1 = (t1/2) / c
    z2 = (t1/2 + interlayer_sep + t2/2) / c
    dz1 = (t1/2) / c
    dz2 = (t2/2) / c
    
    species = [
        p1["metal"], p1["chalcogen"], p1["chalcogen"],
        p2["metal"], p2["chalcogen"], p2["chalcogen"]
    ]
    
    coords = [
        [1/3, 2/3, z1], [1/3, 2/3, z1 + dz1], [1/3, 2/3, z1 - dz1],
        [1/3, 2/3, z2], [1/3, 2/3, z2 + dz2], [1/3, 2/3, z2 - dz2],
    ]
    
    if alignment == "AB":
        coords[3][0] = 2/3
        coords[3][1] = 1/3
        coords[4][0] = 2/3
        coords[4][1] = 1/3
        coords[5][0] = 2/3
        coords[5][1] = 1/3
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "tmd_1": tmd_1,
        "tmd_2": tmd_2,
        "alignment": alignment,
        "lattice_mismatch_percent": round(mismatch, 2),
        "interlayer_separation": interlayer_sep,
        "structure": structure_to_dict(structure, vacuum)
    }
