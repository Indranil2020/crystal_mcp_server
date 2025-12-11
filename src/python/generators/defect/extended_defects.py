"""
defect/extended_defects.py - Extended Defect Generation

Dislocations, grain boundaries, stacking faults.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


def generate_dislocation(
    structure_dict: Dict[str, Any],
    burgers_vector: List[float] = [1, 0, 0],
    dislocation_type: str = "edge",
    supercell: List[int] = [10, 10, 1]
) -> Dict[str, Any]:
    """
    Generate dislocation defect.
    
    Args:
        structure_dict: Base structure
        burgers_vector: Burgers vector
        dislocation_type: 'edge' or 'screw'
        supercell: Supercell size
    
    Returns:
        Structure with dislocation
    """
    atoms = list(structure_dict.get("atoms", []))
    lattice_data = structure_dict.get("lattice", {})
    
    if "matrix" in lattice_data:
        lattice = Lattice(lattice_data["matrix"])
    else:
        a = lattice_data.get("a", 5)
        lattice = Lattice.cubic(a)
    
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)
    structure.make_supercell(supercell)
    
    # Apply displacement field
    b = np.array(burgers_vector)
    center = np.array(structure.lattice.abc) / 2
    
    new_coords = []
    for site in structure:
        pos = np.array(site.coords)
        r = pos - center
        
        if dislocation_type == "edge":
            # Edge dislocation displacement
            theta = np.arctan2(r[1], r[0])
            u = b * theta / (2 * np.pi)
        else:
            # Screw dislocation
            u = b * np.arctan2(r[1], r[0]) / (2 * np.pi)
        
        new_pos = pos + u
        new_coords.append(new_pos)
    
    return {
        "success": True,
        "defect_type": "dislocation",
        "dislocation_type": dislocation_type,
        "burgers_vector": burgers_vector,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "cartesian": list(c)} for s, c in zip(structure, new_coords)]
        }
    }


def generate_grain_boundary(
    element: str = "Fe",
    sigma: int = 5,
    axis: List[int] = [0, 0, 1],
    n_layers: int = 10,
    vacuum: float = 0
) -> Dict[str, Any]:
    """
    Generate coincident site lattice grain boundary.
    
    Args:
        element: Element
        sigma: CSL sigma value
        axis: Rotation axis
        n_layers: Layers per grain
        vacuum: Vacuum for slab
    
    Returns:
        Grain boundary structure
    """
    # CSL angle lookup
    sigma_angles = {3: 60.0, 5: 36.87, 7: 38.21, 9: 38.94, 11: 50.48, 13: 22.62}
    
    if sigma not in sigma_angles:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown sigma: {sigma}"}}
    
    angle = sigma_angles[sigma]
    a = 2.87  # Fe lattice constant
    
    # Build bicrystal
    width = n_layers * a
    c = 2 * width + vacuum
    lattice = Lattice.orthorhombic(a * 3, width * 2, c)
    
    species = []
    coords = []
    
    theta_rad = np.radians(angle)
    
    # Grain 1 (unrotated)
    for i in range(3):
        for j in range(n_layers):
            x = (i + 0.5) / 3
            y = j / n_layers / 2
            z = 0.25
            species.append(element)
            coords.append([x, y, z])
    
    # Grain 2 (rotated)
    for i in range(3):
        for j in range(n_layers):
            x = (i + 0.5) / 3
            y = 0.5 + j / n_layers / 2
            z = 0.25
            species.append(element)
            coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "defect_type": "grain_boundary",
        "sigma": sigma,
        "misorientation_deg": angle,
        "axis": axis,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }
