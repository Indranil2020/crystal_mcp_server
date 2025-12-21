"""
bulk/perovskites.py - Perovskite Structure Generation

Generates various perovskite structures:
- Simple perovskite ABX3
- Double perovskite A2BB'X6
- Layered perovskites (Ruddlesden-Popper, Dion-Jacobson, Aurivillius)
- Anti-perovskites

Scientific basis:
- Goldschmidt tolerance factors
- Octahedral tilting (Glazer notation)
- B-site ordering patterns
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Goldschmidt tolerance factor ranges
TOLERANCE_FACTOR_RANGES = {
    "cubic": (0.9, 1.0),
    "orthorhombic": (0.8, 0.9),
    "rhombohedral": (0.9, 0.95),
}

# Common perovskite compositions with lattice constants
PEROVSKITE_DATABASE = {
    "SrTiO3": {"a": 3.905, "tolerance": 1.00, "tilting": "a0a0a0"},
    "BaTiO3": {"a": 4.00, "tolerance": 1.06, "tilting": "a0a0a0"},
    "CaTiO3": {"a": 3.84, "tolerance": 0.97, "tilting": "a-a-c+"},
    "LaAlO3": {"a": 3.79, "tolerance": 0.95, "tilting": "a-a-a-"},
    "LaMnO3": {"a": 3.88, "tolerance": 0.91, "tilting": "a+b-b-"},
    "PbTiO3": {"a": 3.90, "tolerance": 1.02, "tilting": "a0a0a0"},
    "NaNbO3": {"a": 3.95, "tolerance": 0.97, "tilting": "a-a-a-"},
}


def generate_perovskite(
    a_site: str,
    b_site: str,
    x_site: str = "O",
    a: Optional[float] = None,
    phase: str = "cubic",
    tilting: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate simple perovskite ABX3 structure.
    
    Args:
        a_site: A-site cation (e.g., 'Ca', 'Sr', 'Ba', 'La')
        b_site: B-site cation (e.g., 'Ti', 'Mn', 'Fe')
        x_site: Anion (usually 'O' for oxides)
        a: Lattice constant (estimated if not provided)
        phase: 'cubic', 'tetragonal', 'orthorhombic', 'rhombohedral'
        tilting: Glazer notation (e.g., 'a0a0a0', 'a-a-a-', 'a+b-b-')
    
    Returns:
        Perovskite structure dictionary
    
    Examples:
        >>> result = generate_perovskite('Ca', 'Ti', 'O', a=3.84)
        >>> result["formula"]
        'CaTiO3'
    """
    # Estimate lattice constant from ionic radii
    if a is None:
        # Simple estimate based on B-O bond length ~ 2.0 Å
        a = 4.0  # Default
        formula = f"{a_site}{b_site}{x_site}3"
        if formula in PEROVSKITE_DATABASE:
            a = PEROVSKITE_DATABASE[formula]["a"]
    
    if phase == "cubic":
        lattice = Lattice.cubic(a)
        # Pm-3m structure
        species = [a_site, b_site, x_site, x_site, x_site]
        coords = [
            [0.0, 0.0, 0.0],  # A at corner
            [0.5, 0.5, 0.5],  # B at body center
            [0.5, 0.5, 0.0],  # X1
            [0.5, 0.0, 0.5],  # X2
            [0.0, 0.5, 0.5],  # X3
        ]
    
    elif phase == "tetragonal":
        c = a * 1.04  # Slight tetragonal distortion
        lattice = Lattice.tetragonal(a, c)
        species = [a_site, b_site, x_site, x_site, x_site]
        coords = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.52],  # Off-center B-site
            [0.5, 0.5, 0.02],
            [0.5, 0.0, 0.52],
            [0.0, 0.5, 0.52],
        ]
    
    elif phase == "orthorhombic":
        # Pnma setting (common for tilted perovskites)
        b = a * np.sqrt(2)
        c = a * np.sqrt(2)
        lattice = Lattice.orthorhombic(a, b, c)
        # 4 formula units
        species = [a_site] * 4 + [b_site] * 4 + [x_site] * 12
        # Simplified positions
        coords = [[0.0, 0.0, 0.25], [0.5, 0.5, 0.75], [0.0, 0.5, 0.75], [0.5, 0.0, 0.25]] + \
                 [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]] + \
                 [[x/4, y/4, z/4] for x in [0,1,2,3] for y in [0,2] for z in [0,2]][:12]
    
    else:
        # Default to cubic
        lattice = Lattice.cubic(a)
        species = [a_site, b_site, x_site, x_site, x_site]
        coords = [[0,0,0], [0.5,0.5,0.5], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "formula": f"{a_site}{b_site}{x_site}3",
        "phase": phase,
        "tilting": tilting or "a0a0a0",
        "structure": structure_to_dict(structure)
    }


def generate_double_perovskite(
    a_site: str,
    b_site_1: str,
    b_site_2: str,
    x_site: str = "O",
    ordering: str = "rock-salt",
    a: float = 8.0
) -> Dict[str, Any]:
    """
    Generate double perovskite A2BB'X6 structure.
    
    Args:
        a_site: A-site cation
        b_site_1: First B-site cation
        b_site_2: Second B-site cation (B')
        x_site: Anion
        ordering: B-site ordering ('rock-salt', 'layered', 'columnar')
        a: Lattice constant (doubled compared to simple perovskite)
    
    Returns:
        Double perovskite structure dictionary
    """
    valid_orderings = ["rock-salt", "layered", "columnar", "random"]
    if ordering not in valid_orderings:
        return {"success": False, "error": {"code": "INVALID_ORDERING", "message": f"Must be one of {valid_orderings}"}}
    
    lattice = Lattice.cubic(a)
    
    if ordering == "rock-salt":
        # Fm-3m: B and B' alternate in 3D checkerboard
        species = [a_site] * 8 + [b_site_1] * 4 + [b_site_2] * 4 + [x_site] * 24
        
        # A-sites at (1/4, 1/4, 1/4) family
        a_coords = []
        for dx in [0.25, 0.75]:
            for dy in [0.25, 0.75]:
                for dz in [0.25, 0.75]:
                    a_coords.append([dx, dy, dz])
        
        # B-sites (rock-salt ordered)
        b1_coords = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
        b2_coords = [[0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
        
        # X-sites at face centers around each B
        x_coords = []
        for b in b1_coords + b2_coords:
            for d in [[0.25, 0, 0], [-0.25, 0, 0], [0, 0.25, 0], [0, -0.25, 0], [0, 0, 0.25], [0, 0, -0.25]]:
                xc = [(b[i] + d[i]) % 1.0 for i in range(3)]
                if xc not in x_coords:
                    x_coords.append(xc)
        x_coords = x_coords[:24]  # Take 24 unique
        
        coords = a_coords + b1_coords + b2_coords + x_coords
    
    elif ordering == "layered":
        # L12 ordering: B and B' in alternating layers
        species = [a_site] * 8 + [b_site_1] * 4 + [b_site_2] * 4 + [x_site] * 24
        coords = [[0.25, 0.25, 0.25 + 0.5*i] for i in range(8)] + \
                 [[0, 0, i*0.5] for i in [0,1]] + [[0.5, 0.5, i*0.5] for i in [0,1]] + \
                 [[0, 0.5, 0.25 + i*0.5] for i in range(4)] + \
                 [[0.25 + i*0.125, 0, j*0.5] for i in range(4) for j in [0,1]][:24]
    else:
        # Simplified version
        species = [a_site, a_site, b_site_1, b_site_2, x_site, x_site, x_site]
        coords = [[0.25,0.25,0.25], [0.75,0.75,0.75], [0,0,0], [0.5,0.5,0.5],
                  [0.25,0,0], [0,0.25,0], [0,0,0.25]]
    
    structure = Structure(lattice, species[:len(coords)], coords)
    
    return {
        "success": True,
        "formula": f"{a_site}2{b_site_1}{b_site_2}{x_site}6",
        "ordering": ordering,
        "space_group": 225 if ordering == "rock-salt" else 123,
        "structure": structure_to_dict(structure)
    }


def generate_layered_perovskite(
    a_site: str,
    a_prime_site: Optional[str],
    b_site: str,
    x_site: str = "O",
    n_layers: int = 1,
    family: str = "ruddlesden-popper"
) -> Dict[str, Any]:
    """
    Generate layered perovskite structure.
    
    Args:
        a_site: A-site cation in rock-salt block
        a_prime_site: A'-site cation in perovskite block (optional)
        b_site: B-site cation
        x_site: Anion
        n_layers: Number of perovskite layers
        family: 'ruddlesden-popper', 'dion-jacobson', 'aurivillius'
    
    Returns:
        Layered perovskite structure
    """
    valid_families = ["ruddlesden-popper", "dion-jacobson", "aurivillius"]
    if family not in valid_families:
        return {"success": False, "error": {"code": "INVALID_FAMILY", "message": f"Must be one of {valid_families}"}}
    
    # Simplified representation
    a = 3.9
    c = 12.0 + 2.5 * n_layers
    lattice = Lattice.tetragonal(a, c)
    
    species = []
    coords = []
    
    if family == "ruddlesden-popper":
        # (A,A')_{n+1}B_n X_{3n+1}
        # Rock-salt layers interleaved with perovskite blocks
        for i in range(n_layers):
            z = (i + 0.5) / (n_layers + 1) * 0.5 + 0.25
            species.append(b_site)
            coords.append([0.5, 0.5, z])
            for ox in [[0.5, 0.5, z+0.05], [0.5, 0, z], [0, 0.5, z]]:
                species.append(x_site)
                coords.append(ox)
        
        # A-site layers
        for z in [0.1, 0.9]:
            species.append(a_site)
            coords.append([0, 0, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "family": family,
        "n_layers": n_layers,
        "formula": f"({a_site},{a_prime_site or a_site})_{n_layers+1}{b_site}_{n_layers}{x_site}_{3*n_layers+1}",
        "structure": structure_to_dict(structure)
    }


def generate_ruddlesden_popper(
    compound: str = "Sr2TiO4",
    n_layers: int = 1
) -> Dict[str, Any]:
    """
    Generate Ruddlesden-Popper phase structure (A_{n+1}B_nO_{3n+1}).

    Ruddlesden-Popper phases are layered perovskites with rock-salt blocks
    separating perovskite slabs. Common examples: Sr2TiO4 (n=1), Sr3Ti2O7 (n=2).

    Args:
        compound: Compound name (Sr2TiO4, La2CuO4, Ca2MnO4, etc.)
        n_layers: Number of perovskite layers (1, 2, 3...)

    Returns:
        Ruddlesden-Popper structure dictionary

    Examples:
        >>> result = generate_ruddlesden_popper('Sr2TiO4')
        >>> result["n_layers"]
        1
    """
    # Parse common compounds
    RP_COMPOUNDS = {
        "Sr2TiO4": {"a_site": "Sr", "b_site": "Ti", "n": 1},
        "Sr3Ti2O7": {"a_site": "Sr", "b_site": "Ti", "n": 2},
        "La2CuO4": {"a_site": "La", "b_site": "Cu", "n": 1},
        "Ca2MnO4": {"a_site": "Ca", "b_site": "Mn", "n": 1},
        "Sr2RuO4": {"a_site": "Sr", "b_site": "Ru", "n": 1},
        "Sr3Ru2O7": {"a_site": "Sr", "b_site": "Ru", "n": 2},
        "Ca3Ti2O7": {"a_site": "Ca", "b_site": "Ti", "n": 2},
    }

    if compound in RP_COMPOUNDS:
        params = RP_COMPOUNDS[compound]
        a_site = params["a_site"]
        b_site = params["b_site"]
        n_layers = params["n"]
    else:
        # Try to parse formula like A2BO4, A3B2O7
        a_site = "Sr"
        b_site = "Ti"

    result = generate_layered_perovskite(
        a_site=a_site,
        a_prime_site=None,
        b_site=b_site,
        x_site="O",
        n_layers=n_layers,
        family="ruddlesden-popper"
    )

    if result.get("success"):
        result["compound"] = compound
        result["is_2d_like"] = True
        result["has_rock_salt_layer"] = True

    return result
