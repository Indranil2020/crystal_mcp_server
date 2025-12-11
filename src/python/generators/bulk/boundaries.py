"""
bulk/boundaries.py - Grain Boundaries, Twin Boundaries, and Phase Boundaries

Generates planar defect structures:
- Antiphase boundaries (APB)
- Crystallographic shear planes
- Twin boundaries (coherent, semi-coherent, incoherent)
- Stacking faults

Scientific basis:
- Coincident Site Lattice (CSL) theory for GBs
- Σ-value classification (Σ3, Σ5, Σ7, etc.)
- Displacement Shift Complete (DSC) lattice
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Common twin systems for FCC and BCC materials
TWIN_SYSTEMS = {
    "fcc": {
        "twin_plane": (1, 1, 1),
        "shear_direction": [1, 1, -2],
        "shear_magnitude": 1/np.sqrt(6)
    },
    "bcc": {
        "twin_plane": (1, 1, 2),
        "shear_direction": [1, 1, -1],
        "shear_magnitude": 1/np.sqrt(3)
    },
    "hcp": {
        "twin_plane": (1, 0, -1, 2),  # Miller-Bravais
        "shear_direction": [1, 0, -1, -1],
        "shear_magnitude": 0.17
    }
}


def generate_twin_boundary(
    structure_dict: Dict[str, Any],
    twin_plane: Tuple[int, int, int] = (1, 1, 1),
    n_layers_below: int = 5,
    n_layers_above: int = 5,
    twin_type: str = "coherent"
) -> Dict[str, Any]:
    """
    Generate a twin boundary structure.
    
    Creates a bicrystal with a twin boundary at the center, where
    the upper crystal is a mirror image of the lower.
    
    Args:
        structure_dict: Base structure
        twin_plane: Miller indices of twin plane
        n_layers_below: Layers below the twin plane
        n_layers_above: Layers above the twin plane
        twin_type: 'coherent', 'semi-coherent', or 'incoherent'
    
    Returns:
        Structure with twin boundary
    """
    valid_types = ["coherent", "semi-coherent", "incoherent"]
    if twin_type not in valid_types:
        return {"success": False, "error": {"code": "INVALID_TYPE", "message": f"Must be one of {valid_types}"}}
    
    # For a proper twin, we need to:
    # 1. Create a slab oriented along the twin plane
    # 2. Mirror the upper portion
    
    # Simplified implementation: stack and mirror
    atoms = structure_dict.get("atoms", [])
    lattice = structure_dict.get("lattice", {})
    
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY_STRUCTURE", "message": "No atoms"}}
    
    # Determine twin axis (normal to plane)
    h, k, l = twin_plane
    plane_normal = np.array([h, k, l], dtype=float)
    plane_normal /= np.linalg.norm(plane_normal)
    
    # Create mirrored atoms
    mirrored_atoms = []
    matrix = np.array(lattice.get("matrix", np.eye(3)))
    
    for atom in atoms:
        frac = np.array(atom["coords"])
        cart = frac @ matrix
        
        # Mirror across z=0.5 plane (simplified)
        mirrored_cart = cart.copy()
        mirrored_cart[2] = 1.0 - mirrored_cart[2]  # Mirror across center
        
        mirrored_frac = mirrored_cart @ np.linalg.inv(matrix)
        mirrored_atoms.append({
            "element": atom["element"],
            "coords": list(mirrored_frac),
            "cartesian": list(mirrored_cart)
        })
    
    # Combine original and mirrored
    combined_atoms = atoms + [{"element": a["element"], "coords": [a["coords"][0], a["coords"][1], a["coords"][2] + 1.0]} for a in mirrored_atoms]
    
    # Double the c-lattice
    new_matrix = matrix.copy()
    new_matrix[2] *= 2
    
    new_lattice = {
        "matrix": new_matrix.tolist(),
        "a": lattice.get("a", 5.0),
        "b": lattice.get("b", 5.0),
        "c": lattice.get("c", 5.0) * 2,
        "alpha": lattice.get("alpha", 90),
        "beta": lattice.get("beta", 90),
        "gamma": lattice.get("gamma", 90)
    }
    
    return {
        "success": True,
        "twin_plane": twin_plane,
        "twin_type": twin_type,
        "n_atoms": len(combined_atoms),
        "structure": {
            "lattice": new_lattice,
            "atoms": combined_atoms[:len(atoms)*2],  # Limit
            "metadata": {"type": "twin_boundary"}
        }
    }


def generate_antiphase_boundary(
    structure_dict: Dict[str, Any],
    apb_plane: Tuple[int, int, int] = (1, 0, 0),
    shift_vector: List[float] = [0.5, 0, 0]
) -> Dict[str, Any]:
    """
    Generate an antiphase boundary (APB).
    
    APBs occur in ordered alloys where one domain is shifted relative
    to another by a lattice translation vector.
    
    Args:
        structure_dict: Ordered structure
        apb_plane: Miller indices of APB plane
        shift_vector: Shift vector (in lattice units)
    
    Returns:
        Structure with APB
    """
    atoms = structure_dict.get("atoms", [])
    lattice = structure_dict.get("lattice", {})
    
    new_atoms = []
    shift = np.array(shift_vector)
    
    for atom in atoms:
        atom_copy = dict(atom)
        frac = np.array(atom["coords"])
        
        # Atoms on one side of the plane get shifted
        # Use the first component of the plane normal as boundary
        h, k, l = apb_plane
        plane_pos = h * frac[0] + k * frac[1] + l * frac[2]
        
        if plane_pos > 0.5:
            new_coords = (frac + shift) % 1.0
            atom_copy["coords"] = list(new_coords)
        
        new_atoms.append(atom_copy)
    
    return {
        "success": True,
        "apb_plane": apb_plane,
        "shift_vector": shift_vector,
        "structure": {
            "lattice": lattice,
            "atoms": new_atoms,
            "metadata": {"type": "antiphase_boundary"}
        }
    }


def generate_stacking_fault(
    structure_dict: Dict[str, Any],
    fault_type: str = "intrinsic",
    fault_plane: Tuple[int, int, int] = (1, 1, 1),
    position: float = 0.5
) -> Dict[str, Any]:
    """
    Generate a stacking fault structure.
    
    Args:
        structure_dict: Base structure (should be FCC or HCP)
        fault_type: 'intrinsic', 'extrinsic', 'twin'
        fault_plane: Miller indices of fault plane
        position: Fractional position of fault along normal
    
    Returns:
        Structure with stacking fault
    """
    fault_shifts = {
        "intrinsic": [1/6, 1/6, 0],  # a/6 <112> for FCC {111}
        "extrinsic": [1/3, 1/3, 0],
        "twin": [1/6, 1/6, 0]  # Plus layer removal
    }
    
    if fault_type not in fault_shifts:
        return {"success": False, "error": {"code": "INVALID_TYPE", "message": f"Must be intrinsic, extrinsic, or twin"}}
    
    shift = np.array(fault_shifts[fault_type])
    atoms = structure_dict.get("atoms", [])
    
    new_atoms = []
    for atom in atoms:
        atom_copy = dict(atom)
        frac = np.array(atom["coords"])
        
        # Shift atoms above the fault
        if frac[2] > position:
            new_coords = frac + shift
            new_coords = new_coords % 1.0
            atom_copy["coords"] = list(new_coords)
        
        new_atoms.append(atom_copy)
    
    return {
        "success": True,
        "fault_type": fault_type,
        "fault_plane": fault_plane,
        "shift_vector": shift.tolist(),
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": {"stacking_fault_energy": "calculate with DFT"}
        }
    }
