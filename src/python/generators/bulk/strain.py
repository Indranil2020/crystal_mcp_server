"""
bulk/strain.py - Strain Engineering and Martensitic Transformations

Generates strained structures:
- Uniaxial, biaxial, hydrostatic strain
- Bain path (BCC ↔ FCC)
- Martensitic transformations
- Strain tensor applications

Scientific basis:
- Linear elasticity theory
- Bain strain for martensitic transformations
- Epitaxial strain from lattice mismatch
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


def apply_strain(
    structure_dict: Dict[str, Any],
    strain: float,
    strain_type: str = "hydrostatic"
) -> Dict[str, Any]:
    """
    Apply strain to a structure.
    
    Args:
        structure_dict: Base structure
        strain: Strain magnitude (fractional, e.g., 0.02 for 2%)
        strain_type: 'hydrostatic', 'uniaxial_x', 'uniaxial_z', 'biaxial_xy'
    
    Returns:
        Strained structure
    """
    matrix = np.array(structure_dict["lattice"]["matrix"])
    
    if strain_type == "hydrostatic":
        # Equal strain in all directions
        strain_matrix = np.eye(3) * (1 + strain)
    elif strain_type == "uniaxial_x":
        strain_matrix = np.diag([1 + strain, 1, 1])
    elif strain_type == "uniaxial_z":
        strain_matrix = np.diag([1, 1, 1 + strain])
    elif strain_type == "biaxial_xy":
        strain_matrix = np.diag([1 + strain, 1 + strain, 1])
    elif strain_type == "shear_xy":
        strain_matrix = np.array([[1, strain, 0], [0, 1, 0], [0, 0, 1]])
    else:
        return {"success": False, "error": {"code": "INVALID_TYPE", 
                "message": f"Unknown strain type: {strain_type}"}}
    
    new_matrix = matrix @ strain_matrix
    
    new_lattice = {
        "matrix": new_matrix.tolist(),
        "a": np.linalg.norm(new_matrix[0]),
        "b": np.linalg.norm(new_matrix[1]),
        "c": np.linalg.norm(new_matrix[2]),
        "alpha": structure_dict["lattice"].get("alpha", 90),
        "beta": structure_dict["lattice"].get("beta", 90),
        "gamma": structure_dict["lattice"].get("gamma", 90),
        "volume": abs(np.linalg.det(new_matrix))
    }
    
    return {
        "success": True,
        "strain": strain,
        "strain_type": strain_type,
        "volume_change": (abs(np.linalg.det(new_matrix)) / abs(np.linalg.det(matrix)) - 1) * 100,
        "structure": {
            "lattice": new_lattice,
            "atoms": structure_dict.get("atoms", []),
            "metadata": {"strained": True, "strain_percent": strain * 100}
        }
    }


def apply_strain_tensor(
    structure_dict: Dict[str, Any],
    strain_tensor: List[List[float]]
) -> Dict[str, Any]:
    """
    Apply full 3x3 strain tensor to structure.
    
    Args:
        structure_dict: Base structure
        strain_tensor: 3x3 strain tensor (Voigt or full matrix)
    
    Returns:
        Strained structure
    """
    epsilon = np.array(strain_tensor)
    if epsilon.shape != (3, 3):
        return {"success": False, "error": {"code": "INVALID_TENSOR", "message": "Strain tensor must be 3x3"}}
    
    matrix = np.array(structure_dict["lattice"]["matrix"])
    
    # Apply strain: new_lattice = (I + epsilon) @ old_lattice
    deformation = np.eye(3) + epsilon
    new_matrix = matrix @ deformation
    
    new_lattice = {
        "matrix": new_matrix.tolist(),
        "a": np.linalg.norm(new_matrix[0]),
        "b": np.linalg.norm(new_matrix[1]),
        "c": np.linalg.norm(new_matrix[2]),
        "volume": abs(np.linalg.det(new_matrix))
    }
    
    return {
        "success": True,
        "strain_tensor": epsilon.tolist(),
        "structure": {
            "lattice": new_lattice,
            "atoms": structure_dict.get("atoms", []),
            "metadata": {"strain_tensor_applied": True}
        }
    }


def generate_bain_path(
    element: str = "Fe",
    n_steps: int = 10,
    a_bcc: float = 2.87
) -> Dict[str, Any]:
    """
    Generate structures along Bain transformation path (BCC ↔ FCC).
    
    The Bain path connects BCC and FCC through tetragonal distortion:
    BCC (c/a = 1) → tetragonal → FCC (c/a = √2)
    
    Args:
        element: Element symbol
        n_steps: Number of structures along path
        a_bcc: BCC lattice constant
    
    Returns:
        List of structures along Bain path
    """
    # BCC: c/a = 1
    # FCC: c/a = sqrt(2) ≈ 1.414
    
    c_a_bcc = 1.0
    c_a_fcc = np.sqrt(2)
    
    structures = []
    
    for i in range(n_steps):
        fraction = i / (n_steps - 1)
        c_a = c_a_bcc + fraction * (c_a_fcc - c_a_bcc)
        
        # Adjust a to maintain roughly constant volume
        # V_bcc = a^3, V_tetragonal = a^2 * c
        a = a_bcc * (c_a_bcc / c_a) ** (1/3)
        c = a * c_a
        
        lattice = Lattice.tetragonal(a, c)
        
        # BCC has atoms at (0,0,0) and (0.5,0.5,0.5)
        structure = Structure(lattice, [element, element], [[0,0,0], [0.5, 0.5, 0.5]])
        
        structures.append({
            "step": i,
            "fraction": fraction,
            "c_over_a": c_a,
            "a": a,
            "c": c,
            "structure": structure_to_dict(structure)
        })
    
    return {
        "success": True,
        "element": element,
        "n_steps": n_steps,
        "path": "BCC → FCC (Bain)",
        "structures": structures
    }


def generate_epitaxial_strain(
    structure_dict: Dict[str, Any],
    substrate_a: float,
    substrate_b: Optional[float] = None,
    relax_c: bool = True,
    poisson_ratio: float = 0.3
) -> Dict[str, Any]:
    """
    Apply epitaxial strain from substrate.
    
    Args:
        structure_dict: Film structure
        substrate_a: Substrate a lattice constant
        substrate_b: Substrate b (if different from a)
        relax_c: Whether to relax c-axis (Poisson effect)
        poisson_ratio: Poisson's ratio for c relaxation
    
    Returns:
        Epitaxially strained structure
    """
    lattice = structure_dict["lattice"]
    orig_a = lattice["a"]
    orig_b = lattice.get("b", orig_a)
    orig_c = lattice["c"]
    
    if substrate_b is None:
        substrate_b = substrate_a
    
    strain_a = (substrate_a - orig_a) / orig_a
    strain_b = (substrate_b - orig_b) / orig_b
    
    new_a = substrate_a
    new_b = substrate_b
    
    if relax_c:
        # Poisson effect: strain_c = -2 * nu / (1 - nu) * strain_in_plane
        avg_strain_xy = (strain_a + strain_b) / 2
        strain_c = -2 * poisson_ratio / (1 - poisson_ratio) * avg_strain_xy
        new_c = orig_c * (1 + strain_c)
    else:
        new_c = orig_c
    
    matrix = np.array(lattice.get("matrix", np.diag([orig_a, orig_b, orig_c])))
    scale = np.diag([new_a/orig_a, new_b/orig_b, new_c/orig_c])
    new_matrix = matrix @ scale
    
    return {
        "success": True,
        "strain_a_percent": strain_a * 100,
        "strain_b_percent": strain_b * 100,
        "strain_c_percent": (new_c/orig_c - 1) * 100 if relax_c else 0,
        "structure": {
            "lattice": {
                "matrix": new_matrix.tolist(),
                "a": new_a, "b": new_b, "c": new_c,
                "alpha": lattice.get("alpha", 90),
                "beta": lattice.get("beta", 90),
                "gamma": lattice.get("gamma", 90)
            },
            "atoms": structure_dict.get("atoms", []),
            "metadata": {"epitaxial_strain": True, "substrate_a": substrate_a}
        }
    }
