"""
magnetic/heusler.py - Heusler Alloy Generation

Full and half Heusler alloys for spintronics.
"""

from typing import Dict, Any, List
from pymatgen.core import Structure, Lattice


HEUSLER_DATABASE = {
    # Half-Metallic Heuslers (spintronic)
    "Co2MnSi": {"a": 5.654, "type": "full", "half_metallic": True, "moment": 5.0, "tc": 985},
    "Co2MnGe": {"a": 5.743, "type": "full", "half_metallic": True, "moment": 5.0, "tc": 905},
    "Co2FeSi": {"a": 5.640, "type": "full", "half_metallic": True, "moment": 6.0, "tc": 1100},
    "Co2FeAl": {"a": 5.730, "type": "full", "half_metallic": True, "moment": 5.0, "tc": 1000},
    
    # Shape Memory Heuslers
    "Ni2MnGa": {"a": 5.825, "type": "full", "half_metallic": False, "moment": 4.2, "tc": 376,
                "shape_memory": True, "Ms": 200},
    "Ni2MnIn": {"a": 6.070, "type": "full", "half_metallic": False, "moment": 4.0, "tc": 320,
                "shape_memory": True},
    
    # Topological Heuslers
    "HgTe": {"a": 6.461, "type": "half", "topological": True, "gap_eV": -0.3},
    
    # Thermoelectric Half-Heuslers
    "ZrNiSn": {"a": 6.115, "type": "half", "thermoelectric": True, "zt": 1.0},
    "HfNiSn": {"a": 6.085, "type": "half", "thermoelectric": True, "zt": 1.0},
    "TiCoSb": {"a": 5.885, "type": "half", "thermoelectric": True, "zt": 0.7},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
    }


def generate_heusler(
    formula: str = "Co2MnSi",
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate Heusler alloy structure.
    
    Args:
        formula: Heusler formula
        size: Supercell size
    
    Returns:
        Heusler structure
    """
    if formula not in HEUSLER_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown Heusler: {formula}",
                "available": list(HEUSLER_DATABASE.keys())}}
    
    info = HEUSLER_DATABASE[formula]
    a = info["a"]
    heusler_type = info["type"]
    
    lattice = Lattice.cubic(a)
    
    if heusler_type == "full":
        # X2YZ full Heusler (L21 structure)
        X = formula[:2]
        Y = formula[3:-2] if len(formula) > 5 else formula[2:-2]
        Z = formula[-2:]
        species = [X, X, Y, Z]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0, 0, 0]]
    else:
        # XYZ half Heusler (C1b structure)
        X = formula[:2] if len(formula) > 4 else formula[0]
        Y = formula[2:-2] if len(formula) > 4 else formula[1:-2]
        Z = formula[-2:]
        species = [X, Y, Z]
        coords = [[0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "formula": formula,
        "type": heusler_type,
        "a_angstrom": a,
        "half_metallic": info.get("half_metallic", False),
        "moment_uB": info.get("moment"),
        "tc_kelvin": info.get("tc"),
        "shape_memory": info.get("shape_memory", False),
        "topological": info.get("topological", False),
        "structure": structure_to_dict(structure)
    }
