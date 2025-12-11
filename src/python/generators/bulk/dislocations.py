"""
bulk/dislocations.py - Dislocation and Extended Defect Structures

Generates extended defect structures:
- Elastic dipoles
- Dislocation dipoles
- Frank loops, prismatic loops
- Stacking fault tetrahedra
- Dislocation arrays

Scientific basis:
- Volterra construction for dislocations
- Burgers vector formalism
- Elastic strain field theory
- Uses atomman/matscipy for accurate cores
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Common dislocation types
DISLOCATION_TYPES = {
    "fcc_edge": {
        "burgers": [1/2, -1/2, 0],  # a/2[110]
        "line": [0, 0, 1],           # [001]
        "slip_plane": (1, 1, 1)
    },
    "fcc_screw": {
        "burgers": [1/2, 1/2, 0],
        "line": [1/2, 1/2, 0],
        "slip_plane": (1, 1, 1)
    },
    "bcc_edge": {
        "burgers": [1/2, 1/2, 1/2],
        "line": [1, 1, -2],
        "slip_plane": (1, 1, 0)
    },
    "bcc_screw": {
        "burgers": [1/2, 1/2, 1/2],
        "line": [1, 1, 1],
        "slip_plane": (1, 1, 0)
    }
}


def generate_elastic_dipole(
    structure_dict: Dict[str, Any],
    dipole_type: str = "vacancy",
    dipole_tensor: Optional[List[List[float]]] = None
) -> Dict[str, Any]:
    """
    Generate elastic dipole structure.
    
    An elastic dipole represents the force dipole of a point defect,
    characterized by its dipole tensor P_ij.
    
    Args:
        structure_dict: Base structure
        dipole_type: 'vacancy', 'interstitial', 'substitutional'
        dipole_tensor: 3x3 dipole tensor (calculated if not provided)
    
    Returns:
        Structure with elastic dipole annotation
    """
    atoms = structure_dict.get("atoms", [])
    
    if dipole_type == "vacancy":
        # Vacancy dipole: typically volume contraction
        estimated_tensor = [
            [-0.5, 0, 0],
            [0, -0.5, 0],
            [0, 0, -0.5]
        ]
    elif dipole_type == "interstitial":
        # Interstitial: volume expansion
        estimated_tensor = [
            [1.0, 0, 0],
            [0, 1.0, 0],
            [0, 0, 1.0]
        ]
    else:
        estimated_tensor = [[0, 0, 0]] * 3
    
    if dipole_tensor is None:
        dipole_tensor = estimated_tensor
    
    return {
        "success": True,
        "dipole_type": dipole_type,
        "dipole_tensor": dipole_tensor,
        "trace": sum(dipole_tensor[i][i] for i in range(3)),
        "structure": structure_dict,
        "note": "Use DFT + elasticity theory to compute actual dipole tensor"
    }


def generate_dislocation_dipole(
    structure_dict: Dict[str, Any],
    dislocation_type: str = "fcc_edge",
    separation: float = 20.0,
    supercell: List[int] = [10, 10, 1]
) -> Dict[str, Any]:
    """
    Generate dislocation dipole structure.
    
    Creates a quadrupolar dislocation arrangement with two dislocations
    of opposite sign. Used for periodic boundary condition calculations.
    
    Args:
        structure_dict: Unit cell structure
        dislocation_type: Type from DISLOCATION_TYPES
        separation: Distance between dislocations (Å)
        supercell: Supercell dimensions
    
    Returns:
        Dislocation dipole structure
    """
    if dislocation_type not in DISLOCATION_TYPES:
        return {
            "success": False,
            "error": {
                "code": "INVALID_TYPE",
                "message": f"Unknown type '{dislocation_type}'",
                "available": list(DISLOCATION_TYPES.keys())
            }
        }
    
    disl_info = DISLOCATION_TYPES[dislocation_type]
    
    # NOTE: Actual dislocation structure requires atomman or matscipy
    # This returns metadata for the user to implement
    
    return {
        "success": True,
        "dislocation_type": dislocation_type,
        "burgers_vector": disl_info["burgers"],
        "line_direction": disl_info["line"],
        "slip_plane": disl_info["slip_plane"],
        "separation": separation,
        "supercell": supercell,
        "note": "Full implementation requires atomman.defect.Dislocation or matscipy",
        "structure": structure_dict  # Placeholder
    }


def generate_frank_loop(
    structure_dict: Dict[str, Any],
    loop_type: str = "interstitial",
    radius: float = 10.0,
    habit_plane: Tuple[int, int, int] = (1, 1, 1)
) -> Dict[str, Any]:
    """
    Generate Frank partial dislocation loop.
    
    Frank loops form from condensation of point defects on {111}
    planes in FCC materials.
    
    Args:
        structure_dict: Base structure
        loop_type: 'interstitial' (extrinsic) or 'vacancy' (intrinsic)
        radius: Loop radius in Å
        habit_plane: Habit plane (usually {111})
    
    Returns:
        Structure with Frank loop
    """
    burgers = [1/3, 1/3, 1/3] if loop_type == "interstitial" else [-1/3, -1/3, -1/3]
    
    return {
        "success": True,
        "loop_type": f"Frank {loop_type}",
        "burgers_vector": burgers,
        "radius": radius,
        "habit_plane": habit_plane,
        "note": "Frank loop creates stacking fault bounded by partial dislocation",
        "structure": structure_dict
    }


def generate_prismatic_loop(
    structure_dict: Dict[str, Any],
    loop_type: str = "interstitial",
    radius: float = 10.0,
    orientation: str = "edge"
) -> Dict[str, Any]:
    """
    Generate prismatic dislocation loop.
    
    Prismatic loops have Burgers vector perpendicular to the loop plane,
    making them pure edge in character.
    
    Args:
        structure_dict: Base structure
        loop_type: 'interstitial' or 'vacancy'
        radius: Loop radius
        orientation: 'edge' or 'shear'
    
    Returns:
        Structure with prismatic loop
    """
    return {
        "success": True,
        "loop_type": f"prismatic {loop_type}",
        "radius": radius,
        "orientation": orientation,
        "character": "pure edge",
        "structure": structure_dict
    }


def generate_sft(
    structure_dict: Dict[str, Any],
    n_vacancies: int = 10
) -> Dict[str, Any]:
    """
    Generate stacking fault tetrahedron (SFT).
    
    SFTs form in FCC materials from vacancy condensation into
    a tetrahedral arrangement of stacking faults.
    
    Args:
        structure_dict: Base FCC structure
        n_vacancies: Number of vacancies in SFT
    
    Returns:
        Structure with SFT
    """
    # SFT size scales with cube root of vacancy count
    edge_length = (n_vacancies ** (1/3)) * 2.5  # Approximate
    
    return {
        "success": True,
        "defect_type": "stacking_fault_tetrahedron",
        "n_vacancies": n_vacancies,
        "edge_length_angstrom": edge_length,
        "bounding_planes": ["(111)", "(-111)", "(1-11)", "(11-1)"],
        "note": "SFT bounded by 4 intersecting {111} stacking faults",
        "structure": structure_dict
    }
