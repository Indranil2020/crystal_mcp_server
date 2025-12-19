"""
magnetic/heusler.py - Heusler Alloys

Comprehensive Heusler alloy generation:
- Full Heusler (X2YZ)
- Half Heusler (XYZ)
- Inverse Heusler
- Quaternary Heusler
- Half-metallic ferromagnets
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Full Heusler database (X2YZ, L21 structure, Fm-3m)
FULL_HEUSLER_DATABASE = {
    # Co-based (high Curie temperature)
    "Co2MnSi": {"a": 5.65, "Tc_K": 985, "moment_muB": 5.0, "halfmetal": True},
    "Co2MnGe": {"a": 5.74, "Tc_K": 905, "moment_muB": 5.0, "halfmetal": True},
    "Co2MnGa": {"a": 5.77, "Tc_K": 694, "moment_muB": 4.0, "halfmetal": False},
    "Co2MnSn": {"a": 6.00, "Tc_K": 829, "moment_muB": 5.0, "halfmetal": True},
    "Co2MnAl": {"a": 5.75, "Tc_K": 693, "moment_muB": 4.0, "halfmetal": False},
    "Co2FeSi": {"a": 5.64, "Tc_K": 1100, "moment_muB": 6.0, "halfmetal": True},
    "Co2FeAl": {"a": 5.73, "Tc_K": 1000, "moment_muB": 5.0, "halfmetal": True},
    "Co2FeGe": {"a": 5.75, "Tc_K": 1000, "moment_muB": 5.5, "halfmetal": False},
    "Co2CrAl": {"a": 5.78, "Tc_K": 330, "moment_muB": 3.0, "halfmetal": True},
    
    # Ni-based
    "Ni2MnGa": {"a": 5.82, "Tc_K": 376, "moment_muB": 4.1, "shape_memory": True},
    "Ni2MnIn": {"a": 6.07, "Tc_K": 314, "moment_muB": 4.0, "shape_memory": True},
    "Ni2MnSn": {"a": 6.05, "Tc_K": 340, "moment_muB": 4.0, "shape_memory": True},
    "Ni2FeGa": {"a": 5.76, "Tc_K": 600, "moment_muB": 3.5, "shape_memory": False},
    
    # Cu-based
    "Cu2MnAl": {"a": 5.95, "Tc_K": 630, "moment_muB": 4.0, "halfmetal": False},
    "Cu2MnIn": {"a": 6.20, "Tc_K": 500, "moment_muB": 4.0, "halfmetal": False},
    
    # Fe-based
    "Fe2VAl": {"a": 5.76, "Tc_K": 0, "moment_muB": 0, "thermoelectric": True},
    "Fe2TiSn": {"a": 6.07, "Tc_K": 0, "moment_muB": 0, "thermoelectric": True},
}


# Half Heusler database (XYZ, C1b structure, F-43m)
HALF_HEUSLER_DATABASE = {
    # Semiconducting
    "NiMnSb": {"a": 5.92, "Tc_K": 730, "moment_muB": 4.0, "halfmetal": True, "bandgap_eV": 0},
    "PtMnSb": {"a": 6.21, "Tc_K": 572, "moment_muB": 4.2, "halfmetal": True},
    "CoMnSb": {"a": 5.88, "Tc_K": 478, "moment_muB": 3.0, "halfmetal": True},
    
    # Thermoelectric
    "TiNiSn": {"a": 5.93, "Tc_K": 0, "moment_muB": 0, "bandgap_eV": 0.45, "ZT": 1.0},
    "ZrNiSn": {"a": 6.11, "Tc_K": 0, "moment_muB": 0, "bandgap_eV": 0.50, "ZT": 0.7},
    "HfNiSn": {"a": 6.08, "Tc_K": 0, "moment_muB": 0, "bandgap_eV": 0.52, "ZT": 0.6},
    "TiCoSb": {"a": 5.88, "Tc_K": 0, "moment_muB": 0, "bandgap_eV": 0.95, "ZT": 0.5},
    "ZrCoSb": {"a": 6.07, "Tc_K": 0, "moment_muB": 0, "bandgap_eV": 1.05, "ZT": 0.4},
    
    # Topological
    "LaPtBi": {"a": 6.83, "Tc_K": 0, "moment_muB": 0, "topological": True, "bandgap_eV": 0},
    "LuPtBi": {"a": 6.59, "Tc_K": 0, "moment_muB": 0, "topological": True, "superconducting": True},
    "YPtBi": {"a": 6.65, "Tc_K": 0, "moment_muB": 0, "topological": True},
    "ScPtBi": {"a": 6.50, "Tc_K": 0, "moment_muB": 0, "topological": True},
}


# Quaternary Heusler (XX'YZ)
QUATERNARY_HEUSLER_DATABASE = {
    "CoFeMnSi": {"a": 5.68, "Tc_K": 400, "moment_muB": 4.0, "halfmetal": True},
    "CoFeCrGe": {"a": 5.72, "Tc_K": 350, "moment_muB": 3.0, "halfmetal": False},
    "NiCoMnSn": {"a": 5.95, "Tc_K": 300, "moment_muB": 3.5, "shape_memory": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_full_heusler(
    formula: str = "Co2MnSi",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate full Heusler alloy (X2YZ, L21 structure).
    
    Args:
        formula: Heusler formula from database
        supercell: Supercell dimensions
    
    Returns:
        Heusler structure
    """
    if formula not in FULL_HEUSLER_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_HEUSLER", "message": f"Unknown formula '{formula}'",
                      "available": list(FULL_HEUSLER_DATABASE.keys())}
        }
    
    info = FULL_HEUSLER_DATABASE[formula]
    a = info["a"]
    
    # Parse formula X2YZ
    elements = []
    for i, char in enumerate(formula):
        if char.isupper():
            elem = char
            if i + 1 < len(formula) and formula[i + 1].islower():
                elem += formula[i + 1]
            elements.append(elem)
    
    if len(elements) >= 3:
        X, Y, Z = elements[0], elements[1], elements[2]
    else:
        X, Y, Z = "Co", "Mn", "Si"
    
    lattice = Lattice.cubic(a)
    
    # L21 structure (Fm-3m): X at 8c, Y at 4a, Z at 4b
    # Wyckoff positions
    species = [X, X, X, X, X, X, X, X, Y, Y, Y, Y, Z, Z, Z, Z]
    coords = [
        # X at 8c (1/4, 1/4, 1/4) family
        [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
        [0.75, 0.75, 0.75], [0.25, 0.25, 0.75], [0.25, 0.75, 0.25], [0.75, 0.25, 0.25],
        # Y at 4a (0, 0, 0) family
        [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
        # Z at 4b (1/2, 1/2, 1/2) family
        [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0],
    ]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "formula": formula,
        "heusler_type": "full_L21",
        "spacegroup": "Fm-3m (225)",
        "lattice_constant_A": a,
        "curie_temperature_K": info["Tc_K"],
        "magnetic_moment_muB": info["moment_muB"],
        "is_half_metal": info.get("halfmetal", False),
        "is_shape_memory": info.get("shape_memory", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_half_heusler(
    formula: str = "TiNiSn",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate half Heusler alloy (XYZ, C1b structure).
    
    Args:
        formula: Heusler formula from database
        supercell: Supercell dimensions
    
    Returns:
        Half Heusler structure
    """
    if formula not in HALF_HEUSLER_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_HEUSLER", "message": f"Unknown formula '{formula}'",
                      "available": list(HALF_HEUSLER_DATABASE.keys())}
        }
    
    info = HALF_HEUSLER_DATABASE[formula]
    a = info["a"]
    
    # Parse formula XYZ
    elements = []
    for i, char in enumerate(formula):
        if char.isupper():
            elem = char
            if i + 1 < len(formula) and formula[i + 1].islower():
                elem += formula[i + 1]
            elements.append(elem)
    
    if len(elements) >= 3:
        X, Y, Z = elements[0], elements[1], elements[2]
    else:
        X, Y, Z = "Ti", "Ni", "Sn"
    
    lattice = Lattice.cubic(a)
    
    # C1b structure (F-43m): X at 4a, Y at 4c, Z at 4b
    species = [X, X, X, X, Y, Y, Y, Y, Z, Z, Z, Z]
    coords = [
        # X at 4a (0, 0, 0)
        [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
        # Y at 4c (1/4, 1/4, 1/4)
        [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
        # Z at 4b (1/2, 1/2, 1/2)
        [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0],
    ]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "formula": formula,
        "heusler_type": "half_C1b",
        "spacegroup": "F-43m (216)",
        "lattice_constant_A": a,
        "curie_temperature_K": info.get("Tc_K", 0),
        "magnetic_moment_muB": info.get("moment_muB", 0),
        "bandgap_eV": info.get("bandgap_eV", 0),
        "is_half_metal": info.get("halfmetal", False),
        "is_topological": info.get("topological", False),
        "thermoelectric_ZT": info.get("ZT", 0),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_inverse_heusler(
    formula: str = "Mn2CoAl",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate inverse Heusler alloy (Xa structure).
    
    In inverse Heusler, X occupies 4a and 4d, Y at 4c, Z at 4b.
    """
    # Parse elements from formula
    elements = []
    for i, char in enumerate(formula):
        if char.isupper():
            elem = char
            if i + 1 < len(formula) and formula[i + 1].islower():
                elem += formula[i + 1]
            elements.append(elem)
    
    if len(elements) >= 3:
        X, Y, Z = elements[0], elements[1], elements[2]
    else:
        X, Y, Z = "Mn", "Co", "Al"
    
    a = 5.80  # Approximate
    lattice = Lattice.cubic(a)
    
    # Xa structure
    species = [X, X, X, X, X, X, X, X, Y, Y, Y, Y, Z, Z, Z, Z]
    coords = [
        # X at 4a (0, 0, 0)
        [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
        # X at 4d (3/4, 3/4, 3/4)
        [0.75, 0.75, 0.75], [0.25, 0.25, 0.75], [0.25, 0.75, 0.25], [0.75, 0.25, 0.25],
        # Y at 4c (1/4, 1/4, 1/4)
        [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
        # Z at 4b (1/2, 1/2, 1/2)
        [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0],
    ]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "formula": formula,
        "heusler_type": "inverse_Xa",
        "spacegroup": "F-43m (216)",
        "lattice_constant_A": a,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_quaternary_heusler(
    formula: str = "CoFeMnSi",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate quaternary Heusler alloy (XX'YZ, LiMgPdSn-type).
    """
    if formula not in QUATERNARY_HEUSLER_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_HEUSLER", "message": f"Unknown formula '{formula}'",
                      "available": list(QUATERNARY_HEUSLER_DATABASE.keys())}
        }
    
    info = QUATERNARY_HEUSLER_DATABASE[formula]
    a = info["a"]
    
    # Parse elements
    elements = []
    for i, char in enumerate(formula):
        if char.isupper():
            elem = char
            if i + 1 < len(formula) and formula[i + 1].islower():
                elem += formula[i + 1]
            elements.append(elem)
    
    if len(elements) >= 4:
        X, Xp, Y, Z = elements[0], elements[1], elements[2], elements[3]
    else:
        X, Xp, Y, Z = "Co", "Fe", "Mn", "Si"
    
    lattice = Lattice.cubic(a)
    
    # LiMgPdSn-type structure (F-43m)
    species = [X, X, X, X, Xp, Xp, Xp, Xp, Y, Y, Y, Y, Z, Z, Z, Z]
    coords = [
        # X at 4c
        [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
        # X' at 4d
        [0.75, 0.75, 0.75], [0.25, 0.25, 0.75], [0.25, 0.75, 0.25], [0.75, 0.25, 0.25],
        # Y at 4a
        [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
        # Z at 4b
        [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0],
    ]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "formula": formula,
        "heusler_type": "quaternary",
        "spacegroup": "F-43m (216)",
        "lattice_constant_A": a,
        "curie_temperature_K": info.get("Tc_K", 0),
        "magnetic_moment_muB": info.get("moment_muB", 0),
        "is_half_metal": info.get("halfmetal", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
