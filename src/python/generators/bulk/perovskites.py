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


def generate_dion_jacobson(
    compound: str = "KLaNb2O7",
    n_layers: int = 2,
    a: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate Dion-Jacobson phase structure (A'[A_{n-1}B_nO_{3n+1}]).

    Dion-Jacobson phases are layered perovskites with a single alkali metal
    layer (typically K, Rb, Cs) separating perovskite slabs. Unlike Ruddlesden-
    Popper phases, they have no rock-salt layer, resulting in different
    stacking sequences.

    Args:
        compound: Compound name (KLaNb2O7, RbLaNb2O7, CsLaNb2O7, etc.)
        n_layers: Number of perovskite layers (typically 2-4)
        a: Lattice constant (optional, estimated if not provided)

    Returns:
        Dion-Jacobson structure dictionary

    Examples:
        >>> result = generate_dion_jacobson('KLaNb2O7')
        >>> result["n_layers"]
        2
    """
    # Database of Dion-Jacobson compounds
    DJ_COMPOUNDS = {
        "KLaNb2O7": {"alkali": "K", "a_site": "La", "b_site": "Nb", "n": 2, "a": 3.88, "c": 10.98},
        "RbLaNb2O7": {"alkali": "Rb", "a_site": "La", "b_site": "Nb", "n": 2, "a": 3.88, "c": 11.43},
        "CsLaNb2O7": {"alkali": "Cs", "a_site": "La", "b_site": "Nb", "n": 2, "a": 3.88, "c": 11.93},
        "KCa2Nb3O10": {"alkali": "K", "a_site": "Ca", "b_site": "Nb", "n": 3, "a": 3.86, "c": 14.58},
        "RbCa2Nb3O10": {"alkali": "Rb", "a_site": "Ca", "b_site": "Nb", "n": 3, "a": 3.86, "c": 15.01},
        "KSr2Nb3O10": {"alkali": "K", "a_site": "Sr", "b_site": "Nb", "n": 3, "a": 3.90, "c": 15.20},
        "KLaTa2O7": {"alkali": "K", "a_site": "La", "b_site": "Ta", "n": 2, "a": 3.88, "c": 10.96},
    }

    if compound in DJ_COMPOUNDS:
        params = DJ_COMPOUNDS[compound]
        alkali = params["alkali"]
        a_site = params["a_site"]
        b_site = params["b_site"]
        n_layers = params["n"]
        a_param = a or params.get("a", 3.88)
        c_param = params.get("c", 10.0 + 3.0 * n_layers)
    else:
        # Default parsing
        alkali = "K"
        a_site = "La"
        b_site = "Nb"
        a_param = a or 3.88
        c_param = 10.0 + 3.0 * n_layers

    lattice = Lattice.tetragonal(a_param, c_param)

    species = []
    coords = []

    # Build Dion-Jacobson structure
    # Single alkali layer at z=0
    species.append(alkali)
    coords.append([0.0, 0.0, 0.0])

    # Perovskite layers
    layer_spacing = 0.8 / n_layers
    for i in range(n_layers):
        z_base = 0.1 + i * layer_spacing

        # B-site (octahedral center)
        species.append(b_site)
        coords.append([0.5, 0.5, z_base])

        # Apical oxygens
        species.append("O")
        coords.append([0.5, 0.5, z_base - 0.06])
        species.append("O")
        coords.append([0.5, 0.5, z_base + 0.06])

        # Equatorial oxygens
        species.append("O")
        coords.append([0.5, 0.0, z_base])
        species.append("O")
        coords.append([0.0, 0.5, z_base])

        # A-site between layers (if n > 1 and not last layer)
        if i < n_layers - 1:
            species.append(a_site)
            coords.append([0.0, 0.0, z_base + layer_spacing / 2])

    structure = Structure(lattice, species, coords)

    formula = f"{alkali}{a_site}_{n_layers-1 if n_layers > 1 else ''}{b_site}_{n_layers}O_{3*n_layers+1}"

    return {
        "success": True,
        "compound": compound,
        "family": "dion-jacobson",
        "n_layers": n_layers,
        "alkali_element": alkali,
        "a_site_element": a_site,
        "b_site_element": b_site,
        "formula": formula,
        "a_angstrom": a_param,
        "c_angstrom": c_param,
        "is_2d_like": True,
        "has_rock_salt_layer": False,
        "ion_exchange_capable": True,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure),
        "description": "Dion-Jacobson layered perovskite with single alkali interlayer"
    }


def generate_aurivillius(
    compound: str = "Bi4Ti3O12",
    n_layers: int = 3,
    a: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate Aurivillius phase structure ([Bi2O2][A_{n-1}B_nO_{3n+1}]).

    Aurivillius phases are layered perovskites with [Bi2O2]^{2+} fluorite-like
    layers separating perovskite slabs. They are important ferroelectric and
    piezoelectric materials.

    Args:
        compound: Compound name (Bi4Ti3O12, SrBi2Ta2O9, Bi2WO6, etc.)
        n_layers: Number of perovskite layers (1-5)
        a: Lattice constant (optional)

    Returns:
        Aurivillius structure dictionary

    Examples:
        >>> result = generate_aurivillius('Bi4Ti3O12')
        >>> result["n_layers"]
        3
    """
    # Database of Aurivillius compounds
    AURIVILLIUS_COMPOUNDS = {
        # n=1 phases
        "Bi2WO6": {"a_site": None, "b_site": "W", "n": 1, "a": 5.46, "c": 16.43, "ferroelectric": True},
        "Bi2MoO6": {"a_site": None, "b_site": "Mo", "n": 1, "a": 5.49, "c": 16.96},
        # n=2 phases
        "SrBi2Ta2O9": {"a_site": "Sr", "b_site": "Ta", "n": 2, "a": 5.53, "c": 25.0, "ferroelectric": True, "fatigue_free": True},
        "SrBi2Nb2O9": {"a_site": "Sr", "b_site": "Nb", "n": 2, "a": 5.51, "c": 25.1, "ferroelectric": True},
        "CaBi2Nb2O9": {"a_site": "Ca", "b_site": "Nb", "n": 2, "a": 5.48, "c": 24.9},
        "BaBi2Nb2O9": {"a_site": "Ba", "b_site": "Nb", "n": 2, "a": 5.56, "c": 25.5},
        "PbBi2Nb2O9": {"a_site": "Pb", "b_site": "Nb", "n": 2, "a": 5.53, "c": 25.2},
        # n=3 phases
        "Bi4Ti3O12": {"a_site": "Bi", "b_site": "Ti", "n": 3, "a": 5.45, "c": 32.8, "ferroelectric": True, "high_Tc": True},
        "SrBi4Ti4O15": {"a_site": "Sr", "b_site": "Ti", "n": 4, "a": 5.44, "c": 41.0, "ferroelectric": True},
        # n=4 phases
        "Sr2Bi4Ti5O18": {"a_site": "Sr", "b_site": "Ti", "n": 5, "a": 5.44, "c": 49.0},
    }

    if compound in AURIVILLIUS_COMPOUNDS:
        params = AURIVILLIUS_COMPOUNDS[compound]
        a_site = params.get("a_site")
        b_site = params["b_site"]
        n_layers = params["n"]
        a_param = a or params.get("a", 5.45)
        c_param = params.get("c", 16.0 + 8.0 * n_layers)
        is_ferroelectric = params.get("ferroelectric", False)
        is_fatigue_free = params.get("fatigue_free", False)
    else:
        a_site = "Sr"
        b_site = "Ti"
        a_param = a or 5.45
        c_param = 16.0 + 8.0 * n_layers
        is_ferroelectric = True
        is_fatigue_free = False

    # Aurivillius phases are typically orthorhombic
    # Simplified as tetragonal for demonstration
    lattice = Lattice.tetragonal(a_param, c_param)

    species = []
    coords = []

    # Calculate proper z-positions based on c-axis length
    # Bi2O2 layer takes about 4 Angstrom, perovskite layer ~4 Angstrom each
    bi2o2_height = 4.0 / c_param  # fractional height of Bi2O2 layer
    perov_height = 4.0 / c_param  # fractional height per perovskite layer

    # Bi2O2 fluorite layer centered at z=0.5
    bi_offset = 1.5 / c_param  # ~1.5 A from center
    species.extend(["Bi", "Bi"])
    coords.append([0.0, 0.0, 0.5 - bi_offset])
    coords.append([0.5, 0.5, 0.5 + bi_offset])

    # O in Bi2O2 layer
    species.extend(["O", "O"])
    coords.append([0.5, 0.0, 0.5])
    coords.append([0.0, 0.5, 0.5])

    # Perovskite layers below and above Bi2O2
    # For n layers, place them symmetrically
    if n_layers >= 1:
        # Calculate starting z for perovskite block (below Bi2O2)
        perov_start = 0.5 - bi2o2_height - n_layers * perov_height

        for i in range(n_layers):
            z_b = perov_start + (i + 0.5) * perov_height
            z_b = z_b % 1.0  # Wrap to [0, 1)

            # B-site
            species.append(b_site)
            coords.append([0.0, 0.0, z_b])

            # Equatorial oxygens
            ox_z = z_b
            species.extend(["O", "O"])
            coords.append([0.5, 0.0, ox_z])
            coords.append([0.0, 0.5, ox_z])

            # Apical oxygen (only between layers and at boundaries)
            if i == 0:  # Bottom apical
                species.append("O")
                coords.append([0.0, 0.0, (z_b - perov_height / 2) % 1.0])

            # Top apical (or between B-sites)
            species.append("O")
            coords.append([0.0, 0.0, (z_b + perov_height / 2) % 1.0])

            # A-site (between perovskite layers)
            if a_site and i < n_layers - 1:
                species.append(a_site)
                coords.append([0.5, 0.5, (z_b + perov_height / 2) % 1.0])

    # Remove any duplicate coordinates
    unique_atoms = []
    unique_coords = []
    for s, c in zip(species, coords):
        is_dup = False
        for uc in unique_coords:
            if all(abs(c[j] - uc[j]) < 0.01 or abs(c[j] - uc[j] - 1) < 0.01 or abs(c[j] - uc[j] + 1) < 0.01 for j in range(3)):
                is_dup = True
                break
        if not is_dup:
            unique_atoms.append(s)
            unique_coords.append(c)

    species = unique_atoms
    coords = unique_coords

    structure = Structure(lattice, species, coords)

    # Build formula
    if n_layers == 1:
        formula = f"Bi2{b_site}O6"
    else:
        formula = f"{a_site if a_site else 'Bi'}Bi2{b_site}{n_layers}O{3*n_layers+3}"

    return {
        "success": True,
        "compound": compound,
        "family": "aurivillius",
        "n_layers": n_layers,
        "formula": formula,
        "a_site_element": a_site,
        "b_site_element": b_site,
        "a_angstrom": a_param,
        "c_angstrom": c_param,
        "is_ferroelectric": is_ferroelectric,
        "is_fatigue_free": is_fatigue_free,
        "has_bi2o2_layer": True,
        "space_group_approx": "I4/mmm or Fmmm (orthorhombic)",
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure),
        "description": "Aurivillius layered ferroelectric with [Bi2O2] interlayer"
    }


def get_layered_perovskite_database() -> Dict[str, Any]:
    """
    Get database of layered perovskite phases.

    Returns database organized by family (Ruddlesden-Popper, Dion-Jacobson,
    Aurivillius) with available compounds and their properties.

    Returns:
        Database dictionary with compounds organized by family
    """
    return {
        "success": True,
        "families": {
            "ruddlesden_popper": {
                "formula": "A_{n+1}B_nO_{3n+1}",
                "interlayer": "Rock-salt AO layer",
                "compounds": ["Sr2TiO4", "Sr3Ti2O7", "La2CuO4", "Ca2MnO4", "Sr2RuO4"],
                "properties": ["superconductivity", "magnetism", "catalysis"]
            },
            "dion_jacobson": {
                "formula": "A'[A_{n-1}B_nO_{3n+1}]",
                "interlayer": "Single alkali layer",
                "compounds": ["KLaNb2O7", "RbLaNb2O7", "CsLaNb2O7", "KCa2Nb3O10"],
                "properties": ["ion_exchange", "photocatalysis", "nanosheet_exfoliation"]
            },
            "aurivillius": {
                "formula": "[Bi2O2][A_{n-1}B_nO_{3n+1}]",
                "interlayer": "Bi2O2 fluorite layer",
                "compounds": ["Bi4Ti3O12", "SrBi2Ta2O9", "Bi2WO6", "SrBi2Nb2O9"],
                "properties": ["ferroelectricity", "piezoelectricity", "fatigue_resistance"]
            }
        },
        "n_layer_range": "Typically n=1 to n=5 perovskite layers"
    }
