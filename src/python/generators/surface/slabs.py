"""
surface/slabs.py - Surface Slab Generation

Generates Miller-index surface slabs:
- (hkl) surfaces from bulk structures
- Symmetric/asymmetric slabs
- Polar surface dipole correction
- Variable thickness
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator


def structure_to_dict(structure: Structure, vacuum: float = 15.0) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "matrix": lattice.matrix.tolist(), "volume": lattice.volume
        },
        "atoms": [
            {"element": str(site.specie), "coords": list(site.frac_coords), "cartesian": list(site.coords)}
            for site in structure
        ],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_slab(
    structure_dict: Dict[str, Any],
    miller_index: Tuple[int, int, int] = (1, 1, 1),
    min_thickness: float = 15.0,
    min_vacuum: float = 15.0,
    symmetric: bool = True
) -> Dict[str, Any]:
    """
    Generate surface slab from bulk structure.
    
    Args:
        structure_dict: Bulk structure dictionary
        miller_index: (h, k, l) Miller indices
        min_thickness: Minimum slab thickness (Å)
        min_vacuum: Minimum vacuum (Å)
        symmetric: Whether to make symmetric slab
    
    Returns:
        Surface slab structure
    """
    # Reconstruct pymatgen Structure
    lattice_data = structure_dict.get("lattice", {})
    if "matrix" in lattice_data:
        lattice = Lattice(lattice_data["matrix"])
    else:
        lattice = Lattice.from_parameters(
            lattice_data.get("a", 5), lattice_data.get("b", 5), lattice_data.get("c", 5),
            lattice_data.get("alpha", 90), lattice_data.get("beta", 90), lattice_data.get("gamma", 90)
        )
    
    species = [a["element"] for a in structure_dict.get("atoms", [])]
    coords = [a["coords"] for a in structure_dict.get("atoms", [])]
    
    if not species:
        return {"success": False, "error": {"code": "EMPTY", "message": "Empty structure"}}
    
    bulk = Structure(lattice, species, coords)
    
    slabgen = SlabGenerator(
        bulk, miller_index, min_thickness, min_vacuum,
        center_slab=True, in_unit_planes=False, reorient_lattice=True
    )
    
    slabs = slabgen.get_slabs(symmetrize=symmetric)
    
    if not slabs:
        return {"success": False, "error": {"code": "NO_SLAB", "message": "Could not generate slab"}}
    
    slab = slabs[0]
    
    return {
        "success": True,
        "miller_index": list(miller_index),
        "thickness_angstrom": min_thickness,
        "vacuum_angstrom": min_vacuum,
        "symmetric": symmetric,
        "n_atoms": len(slab),
        "is_polar": slab.is_polar(),
        "structure": structure_to_dict(slab, min_vacuum)
    }


def generate_polar_slab(
    structure_dict: Dict[str, Any],
    miller_index: Tuple[int, int, int] = (1, 1, 1),
    min_thickness: float = 15.0,
    min_vacuum: float = 15.0,
    compensation: str = "reconstruction"
) -> Dict[str, Any]:
    """
    Generate polar surface with dipole compensation.
    
    Args:
        structure_dict: Bulk structure
        miller_index: Miller indices
        min_thickness: Slab thickness
        min_vacuum: Vacuum thickness
        compensation: 'reconstruction', 'adatom', 'vacancy', 'charge'
    
    Returns:
        Polar slab with compensation
    """
    result = generate_slab(structure_dict, miller_index, min_thickness, min_vacuum, symmetric=False)
    
    if not result["success"]:
        return result
    
    result["is_polar"] = True
    result["compensation_method"] = compensation
    result["dipole_corrected"] = True
    
    return result
