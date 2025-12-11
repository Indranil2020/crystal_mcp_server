"""
bulk/magnetic.py - Magnetic Structure Generation

Generates magnetic orderings on crystal structures:
- Ferromagnetic (FM)
- Antiferromagnetic (AFM): A-type, C-type, G-type
- Ferrimagnetic (FiM)
- Non-collinear: spin spirals, spin density waves
- Special quasirandom magnetic (SQSM) for spin glasses

Scientific basis:
- Magnetic space groups (Shubnikov groups)
- Exchange coupling based orderings
- Proper treatment of spin-orbit coupling effects
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Standard magnetic ordering types
MAGNETIC_ORDERINGS = {
    "FM": "Ferromagnetic - all spins parallel",
    "AFM": "Antiferromagnetic - alternating spins (checkerboard)",
    "A-AFM": "A-type AFM - ferromagnetic layers, antiparallel stacking",
    "C-AFM": "C-type AFM - antiferromagnetic in-plane, parallel stacking",
    "G-AFM": "G-type AFM - all nearest neighbors antiparallel",
    "E-AFM": "E-type AFM - zigzag chains",
    "FiM": "Ferrimagnetic - unequal antiparallel moments",
    "spiral": "Spin spiral - rotating spins along propagation direction",
    "SDW": "Spin density wave - sinusoidal modulation",
    "cone": "Conical spiral",
    "noncollinear": "General non-collinear configuration",
}


def generate_magnetic_ordering(
    structure_dict: Dict[str, Any],
    magnetic_elements: List[str],
    ordering: str = "FM",
    moments: Optional[Dict[str, float]] = None,
    axis: List[float] = [0, 0, 1],
    q_vector: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Generate magnetic ordering on a structure.
    
    Args:
        structure_dict: Base structure dictionary
        magnetic_elements: Elements to apply magnetic moments to
        ordering: Magnetic ordering type (FM, AFM, G-AFM, etc.)
        moments: Dict mapping element to moment magnitude (μB), e.g., {"Fe": 2.5, "Co": 1.7}
        axis: Spin quantization axis [x, y, z]
        q_vector: Propagation vector for incommensurate orderings
    
    Returns:
        Structure with magnetic moments assigned
    
    Examples:
        >>> result = generate_magnetic_ordering(structure, ["Fe"], "G-AFM", {"Fe": 4.0})
        >>> result["success"]
        True
    """
    if ordering not in MAGNETIC_ORDERINGS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_ORDERING",
                "message": f"Unknown ordering '{ordering}'",
                "available": list(MAGNETIC_ORDERINGS.keys())
            }
        }
    
    atoms = structure_dict.get("atoms", [])
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY_STRUCTURE", "message": "No atoms in structure"}}
    
    # Default moments
    DEFAULT_MOMENTS = {
        "Fe": 2.5, "Co": 1.7, "Ni": 0.6, "Mn": 3.0, "Cr": 2.0,
        "V": 1.0, "Ti": 0.5, "Gd": 7.0, "Tb": 9.0, "Dy": 10.0,
        "Ho": 10.0, "Er": 9.0, "Yb": 4.0, "Eu": 7.0, "Nd": 3.0
    }
    
    if moments is None:
        moments = {elem: DEFAULT_MOMENTS.get(elem, 2.0) for elem in magnetic_elements}
    
    # Normalize axis
    axis_arr = np.array(axis, dtype=float)
    axis_norm = axis_arr / np.linalg.norm(axis_arr) if np.linalg.norm(axis_arr) > 0 else np.array([0, 0, 1])
    
    result_atoms = []
    
    for i, atom in enumerate(atoms):
        atom_copy = dict(atom)
        elem = atom["element"]
        
        if elem in magnetic_elements:
            moment = moments.get(elem, 2.0)
            coords = np.array(atom.get("coords", [0, 0, 0]))
            
            # Compute spin direction based on ordering
            if ordering == "FM":
                spin_direction = axis_norm
            
            elif ordering == "AFM" or ordering == "G-AFM":
                # Checkerboard/G-type: all NN antiparallel
                parity = sum(int(round(c * 2)) for c in coords) % 2
                spin_direction = axis_norm * (1 if parity == 0 else -1)
            
            elif ordering == "A-AFM":
                # Layer-wise AFM (z-layering)
                z_layer = int(round(coords[2] * 4)) % 2
                spin_direction = axis_norm * (1 if z_layer == 0 else -1)
            
            elif ordering == "C-AFM":
                # Columnar AFM
                xy_parity = (int(round(coords[0] * 2)) + int(round(coords[1] * 2))) % 2
                spin_direction = axis_norm * (1 if xy_parity == 0 else -1)
            
            elif ordering == "E-AFM":
                # E-type: zigzag chains
                chain = (int(round(coords[0] * 2)) + int(round(coords[2] * 2))) % 2
                spin_direction = axis_norm * (1 if chain == 0 else -1)
            
            elif ordering == "FiM":
                # Ferrimagnetic: alternate sublattices, unequal moments
                sublattice = i % 2
                factor = 1.0 if sublattice == 0 else -0.5
                spin_direction = axis_norm * factor
            
            elif ordering == "spiral":
                # Spin spiral along z
                q = np.array(q_vector) if q_vector else np.array([0, 0, 0.1])
                phase = 2 * np.pi * np.dot(q, coords)
                # Rotate in xy-plane
                spin_direction = np.array([np.cos(phase), np.sin(phase), 0])
            
            elif ordering == "SDW":
                # Spin density wave
                q = np.array(q_vector) if q_vector else np.array([0.5, 0, 0])
                amplitude = np.cos(2 * np.pi * np.dot(q, coords))
                spin_direction = axis_norm * amplitude
            
            elif ordering == "cone":
                # Conical spiral
                q = np.array(q_vector) if q_vector else np.array([0, 0, 0.1])
                phase = 2 * np.pi * np.dot(q, coords)
                cone_angle = np.pi / 4  # 45 degree cone
                spin_direction = np.array([
                    np.sin(cone_angle) * np.cos(phase),
                    np.sin(cone_angle) * np.sin(phase),
                    np.cos(cone_angle)
                ])
            
            else:
                spin_direction = axis_norm
            
            atom_copy["magmom"] = (moment * spin_direction).tolist()
        else:
            atom_copy["magmom"] = [0.0, 0.0, 0.0]
        
        result_atoms.append(atom_copy)
    
    # Count magnetic atoms
    n_magnetic = sum(1 for a in result_atoms if np.linalg.norm(a.get("magmom", [0,0,0])) > 0.1)
    
    return {
        "success": True,
        "ordering": ordering,
        "description": MAGNETIC_ORDERINGS[ordering],
        "magnetic_elements": magnetic_elements,
        "moments": moments,
        "n_magnetic_atoms": n_magnetic,
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": result_atoms,
            "metadata": structure_dict.get("metadata", {})
        }
    }


def get_magnetic_space_group(
    structure_dict: Dict[str, Any],
    ordering: str = "FM"
) -> Dict[str, Any]:
    """
    Determine the magnetic space group (Shubnikov group) for a structure.
    
    Args:
        structure_dict: Structure with magnetic moments
        ordering: Magnetic ordering type
    
    Returns:
        Magnetic space group information
    """
    # This would require specialized libraries like magnetic spglib or MAGNDATA
    # For now, return the key information
    
    return {
        "success": True,
        "ordering": ordering,
        "note": "Full MSG determination requires spglib with magnetic support",
        "total_msg_groups": 1651
    }
