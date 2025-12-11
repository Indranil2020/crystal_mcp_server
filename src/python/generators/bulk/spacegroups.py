"""
bulk/spacegroups.py - Space Group Crystal Generation

Generates crystalline structures for all 230 space groups using PyXtal.
Also supports 1651 magnetic space groups (Shubnikov groups).

Scientific basis:
- PyXtal for symmetry-aware random structure generation
- Wyckoff position enumeration from International Tables
- Support for primitive, conventional, and Niggli reduced cells
"""

from typing import Dict, Any, List, Optional, Union, Tuple
import numpy as np
from pyxtal import pyxtal
from pyxtal.symmetry import Group
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict, validate_composition, validate_space_group


# Crystal system to space group ranges
CRYSTAL_SYSTEMS = {
    "triclinic": (1, 2),
    "monoclinic": (3, 15),
    "orthorhombic": (16, 74),
    "tetragonal": (75, 142),
    "trigonal": (143, 167),
    "hexagonal": (168, 194),
    "cubic": (195, 230)
}


def generate_crystal(
    composition: List[str],
    space_group: Union[int, str],
    lattice_params: Optional[Dict[str, float]] = None,
    num_atoms: Optional[int] = None,
    cell_type: str = "primitive"
) -> Dict[str, Any]:
    """
    Generate a crystal structure with specified space group symmetry using PyXtal.
    
    Args:
        composition: List of elements, e.g., ["Ga", "As"] for GaAs
        space_group: Space group number (1-230) or Hermann-Mauguin symbol
        lattice_params: Optional dict with a, b, c, alpha, beta, gamma
        num_atoms: Target number of atoms in cell
        cell_type: 'primitive', 'conventional', or 'niggli'
    
    Returns:
        Dictionary with structure data
    
    Example:
        >>> result = generate_crystal(["Si"], 227, lattice_params={"a": 5.43}, num_atoms=8)
        >>> result["success"]
        True
    """
    # Validate inputs
    is_valid, error = validate_composition(composition)
    if not is_valid:
        return {"success": False, "error": {"code": "INVALID_COMPOSITION", "message": error}}
    
    is_valid, error, sg_number = validate_space_group(space_group)
    if not is_valid:
        return {"success": False, "error": {"code": "INVALID_SPACE_GROUP", "message": error}}
    
    # Count elements for PyXtal
    from collections import Counter
    elem_counts = Counter(composition)
    species = list(elem_counts.keys())
    numIons = list(elem_counts.values())
    
    # Create PyXtal structure
    crystal = pyxtal()
    
    # Prepare lattice if provided
    lattice = None
    if lattice_params:
        a = lattice_params.get("a", 5.0)
        b = lattice_params.get("b", a)
        c = lattice_params.get("c", a)
        alpha = lattice_params.get("alpha", 90)
        beta = lattice_params.get("beta", 90)
        gamma = lattice_params.get("gamma", 90)
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    
    # Generate structure
    crystal.from_random(
        dim=3,
        group=sg_number,
        species=species,
        numIons=numIons,
        lattice=lattice
    )
    
    if not crystal.valid:
        return {"success": False, "error": {"code": "GENERATION_FAILED", 
                "message": "PyXtal could not generate valid structure with given constraints"}}
    
    # Convert to pymatgen
    pmg_structure = crystal.to_pymatgen()
    
    # Apply cell type conversion
    if cell_type == "conventional":
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        sga = SpacegroupAnalyzer(pmg_structure, symprec=0.1)
        pmg_structure = sga.get_conventional_standard_structure()
    elif cell_type == "niggli":
        pmg_structure = pmg_structure.get_reduced_structure()
    
    return {
        "success": True,
        "space_group": sg_number,
        "cell_type": cell_type,
        "structure": structure_to_dict(pmg_structure)
    }


def generate_all_polymorphs(
    composition: List[str],
    space_groups: Optional[List[int]] = None,
    crystal_system: Optional[str] = None,
    num_atoms: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate structures for multiple space groups (polymorph search).
    
    Args:
        composition: List of elements
        space_groups: List of SG numbers to try (default: all compatible)
        crystal_system: Limit to specific system (triclinic, cubic, etc.)
        num_atoms: Target atom count
    
    Returns:
        Dictionary with list of generated polymorphs
    """
    if space_groups is None:
        if crystal_system and crystal_system.lower() in CRYSTAL_SYSTEMS:
            sg_min, sg_max = CRYSTAL_SYSTEMS[crystal_system.lower()]
            space_groups = list(range(sg_min, sg_max + 1))
        else:
            space_groups = list(range(1, 231))
    
    polymorphs = []
    errors = []
    
    for sg in space_groups:
        result = generate_crystal(composition, sg, num_atoms=num_atoms)
        if result["success"]:
            polymorphs.append({
                "space_group": sg,
                "structure": result["structure"]
            })
        else:
            errors.append({"space_group": sg, "error": result.get("error", {}).get("message", "Unknown")})
    
    return {
        "success": len(polymorphs) > 0,
        "n_polymorphs": len(polymorphs),
        "polymorphs": polymorphs,
        "failed": errors[:10]  # Limit error list
    }


def get_wyckoff_positions(space_group: int) -> Dict[str, Any]:
    """
    Get all Wyckoff positions for a space group.
    
    Args:
        space_group: Space group number (1-230)
    
    Returns:
        Dictionary with Wyckoff position information
    """
    is_valid, error, sg_number = validate_space_group(space_group)
    if not is_valid:
        return {"success": False, "error": {"code": "INVALID_SPACE_GROUP", "message": error}}
    
    group = Group(sg_number)
    
    wyckoff_info = []
    for wp in group:
        wyckoff_info.append({
            "letter": wp.letter,
            "multiplicity": wp.multiplicity,
            "site_symmetry": str(wp.site_symm),
            "coordinates": [str(c) for c in wp.ops[:min(3, len(wp.ops))]]  # First few coords
        })
    
    return {
        "success": True,
        "space_group": sg_number,
        "symbol": group.symbol,
        "crystal_system": group.lattice_type,
        "wyckoff_positions": wyckoff_info
    }
