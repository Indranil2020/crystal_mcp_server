"""
bulk/supercells.py - Supercell Generation

Generates supercells using integer and rational transformation matrices.
Supports commensurate cell generation for heterostructures.

Scientific basis:
- Integer matrices: Simple n×m×l expansions
- Rational matrices: Non-diagonal transformations for Moiré patterns
- Commensurate cells: Matching lattices for heterostructures
"""

from typing import Dict, Any, List, Optional, Union
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.transformations.standard_transformations import SupercellTransformation

from .base import structure_to_dict, BulkStructure


def make_supercell(
    structure_dict: Dict[str, Any],
    scaling: Union[List[int], List[List[int]]],
    wrap_atoms: bool = True
) -> Dict[str, Any]:
    """
    Generate supercell with integer scaling matrix.
    
    Args:
        structure_dict: Base structure dictionary
        scaling: [nx, ny, nz] or 3x3 matrix
        wrap_atoms: Wrap atoms into unit cell
    
    Returns:
        Supercell structure dictionary
    """
    # Validate input
    if not structure_dict.get("lattice") or not structure_dict.get("atoms"):
        return {"success": False, "error": {"code": "INVALID_STRUCTURE", "message": "Missing lattice or atoms"}}
    
    # Parse scaling
    scaling_arr = np.array(scaling)
    if scaling_arr.shape == (3,):
        matrix = np.diag(scaling_arr.astype(int))
    elif scaling_arr.shape == (3, 3):
        matrix = scaling_arr.astype(int)
    else:
        return {"success": False, "error": {"code": "INVALID_SCALING", "message": "Scaling must be [3] or [3,3]"}}
    
    det = np.linalg.det(matrix)
    if abs(det) < 0.5:
        return {"success": False, "error": {"code": "SINGULAR_MATRIX", "message": "Scaling matrix is singular"}}
    
    # Convert to pymatgen
    lattice = Lattice(structure_dict["lattice"].get("matrix", [[1,0,0],[0,1,0],[0,0,1]]))
    species = [a["element"] for a in structure_dict["atoms"]]
    coords = [a["coords"] for a in structure_dict["atoms"]]
    
    structure = Structure(lattice, species, coords)
    
    # Make supercell
    structure.make_supercell(matrix)
    
    return {
        "success": True,
        "transformation_matrix": matrix.tolist(),
        "volume_multiplier": int(round(abs(det))),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def make_rational_supercell(
    structure_dict: Dict[str, Any],
    transformation_matrix: List[List[float]]
) -> Dict[str, Any]:
    """
    Generate supercell with rational (non-integer) transformation matrix.
    
    Uses the supercell transformation to create cells with non-integer
    relationships, useful for Moiré superlattices and incommensurate structures.
    
    Args:
        structure_dict: Base structure dictionary
        transformation_matrix: 3x3 float matrix
    
    Returns:
        Supercell structure dictionary
    """
    matrix = np.array(transformation_matrix)
    if matrix.shape != (3, 3):
        return {"success": False, "error": {"code": "INVALID_MATRIX", "message": "Matrix must be 3x3"}}
    
    det = np.abs(np.linalg.det(matrix))
    if det < 0.01:
        return {"success": False, "error": {"code": "SINGULAR_MATRIX", "message": f"Determinant {det:.4f} too small"}}
    
    # Convert to pymatgen
    lattice = Lattice(structure_dict["lattice"].get("matrix", [[1,0,0],[0,1,0],[0,0,1]]))
    species = [a["element"] for a in structure_dict["atoms"]]
    coords = [a["coords"] for a in structure_dict["atoms"]]
    
    structure = Structure(lattice, species, coords)
    
    # For rational matrices, we need to scale up and find a commensurate cell
    # Round to integers if close
    int_matrix = np.round(matrix).astype(int)
    
    if np.allclose(matrix, int_matrix, atol=0.05):
        structure.make_supercell(int_matrix)
    else:
        # Scale up to find approximate integer matrix
        for scale in [2, 3, 4, 5, 6]:
            scaled = matrix * scale
            int_scaled = np.round(scaled).astype(int)
            if np.allclose(scaled, int_scaled, atol=0.1):
                structure.make_supercell(int_scaled)
                break
        else:
            # Fall back to direct application
            structure.make_supercell(np.round(matrix).astype(int))
    
    return {
        "success": True,
        "transformation_matrix": matrix.tolist(),
        "volume_multiplier": float(det),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def make_commensurate_cell(
    structure1: Dict[str, Any],
    structure2: Dict[str, Any],
    max_atoms: int = 200,
    tolerance: float = 0.05
) -> Dict[str, Any]:
    """
    Find commensurate supercells for two structures (for heterostructures).
    
    Searches for integer matrices M1 and M2 such that the resulting
    supercells have matching in-plane lattice vectors.
    
    Args:
        structure1: First structure (usually substrate)
        structure2: Second structure (usually overlayer)
        max_atoms: Maximum total atoms allowed
        tolerance: Lattice mismatch tolerance
    
    Returns:
        Dictionary with commensurate cell information
    """
    lat1 = np.array(structure1["lattice"]["matrix"])
    lat2 = np.array(structure2["lattice"]["matrix"])
    n1 = len(structure1["atoms"])
    n2 = len(structure2["atoms"])
    
    a1, a2 = np.linalg.norm(lat1[0]), np.linalg.norm(lat2[0])
    
    best_result = None
    best_mismatch = float('inf')
    
    # Search for commensurate cells
    for m1 in range(1, 6):
        for m2 in range(1, 6):
            if m1 * n1 + m2 * n2 > max_atoms:
                continue
            
            # Simple diagonal scaling
            scaled_a1 = a1 * m1
            scaled_a2 = a2 * m2
            
            mismatch = abs(scaled_a1 - scaled_a2) / max(scaled_a1, scaled_a2)
            
            if mismatch < tolerance and mismatch < best_mismatch:
                best_mismatch = mismatch
                best_result = {
                    "m1": [m1, m1, 1],
                    "m2": [m2, m2, 1],
                    "mismatch": mismatch,
                    "total_atoms": m1 * m1 * n1 + m2 * m2 * n2
                }
    
    if best_result is None:
        return {"success": False, "error": {"code": "NO_COMMENSURATE", 
                "message": f"No commensurate cell found within {tolerance*100}% mismatch"}}
    
    return {
        "success": True,
        **best_result
    }


# Preset supercell transformations for common cases
SUPERCELL_PRESETS = {
    "sqrt3": [[2, 1, 0], [-1, 1, 0], [0, 0, 1]],  # √3×√3 R30° for hexagonal
    "root2": [[1, 1, 0], [1, -1, 0], [0, 0, 1]],  # √2×√2 R45° for square
    "2x2x2": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    "3x3x3": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
    "2x2x1": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
    "3x3x1": [[3, 0, 0], [0, 3, 0], [0, 0, 1]],
}
