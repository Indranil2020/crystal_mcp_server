"""
thermoelectric/materials.py - Thermoelectric Materials

Comprehensive thermoelectric material generation:
- Bi2Te3-based alloys
- PbTe and related chalcogenides
- Half-Heusler thermoelectrics
- Skutterudites
- Clathrates
- Oxide thermoelectrics
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Thermoelectric database with ZT values
THERMOELECTRIC_DATABASE = {
    # Bi2Te3 family (room temperature)
    "Bi2Te3": {"ZT": 1.0, "T_opt_K": 300, "type": "p", "a": 4.38, "c": 30.49,
               "structure": "rhombohedral", "band_gap_eV": 0.15},
    "Bi2Se3": {"ZT": 0.3, "T_opt_K": 300, "type": "n", "a": 4.14, "c": 28.64},
    "Sb2Te3": {"ZT": 0.8, "T_opt_K": 300, "type": "p", "a": 4.26, "c": 30.35},
    "Bi0.5Sb1.5Te3": {"ZT": 1.4, "T_opt_K": 300, "type": "p", "alloy": True,
                      "description": "p-type commercial alloy"},
    "Bi2Te2.7Se0.3": {"ZT": 1.1, "T_opt_K": 300, "type": "n", "alloy": True,
                      "description": "n-type commercial alloy"},
    
    # PbTe family (mid temperature)
    "PbTe": {"ZT": 0.8, "T_opt_K": 700, "type": "both", "a": 6.46,
             "structure": "rocksalt", "band_gap_eV": 0.31},
    "PbSe": {"ZT": 1.0, "T_opt_K": 850, "type": "both", "a": 6.12},
    "PbS": {"ZT": 0.7, "T_opt_K": 900, "type": "n", "a": 5.94},
    "PbTe_Na": {"ZT": 2.0, "T_opt_K": 800, "type": "p", "dopant": "Na",
                "description": "Record p-type PbTe"},
    "PbTe_La": {"ZT": 1.4, "T_opt_K": 800, "type": "n", "dopant": "La"},
    "SnTe": {"ZT": 0.5, "T_opt_K": 700, "type": "p", "a": 6.31, "Pb_free": True},
    "GeTe": {"ZT": 1.3, "T_opt_K": 700, "type": "p", "a": 5.98, "rhombohedral": True},
    
    # Half-Heusler (mid-high temperature)
    "TiNiSn": {"ZT": 1.0, "T_opt_K": 700, "type": "n", "a": 5.93,
               "structure": "half_heusler", "band_gap_eV": 0.45},
    "ZrNiSn": {"ZT": 0.7, "T_opt_K": 900, "type": "n", "a": 6.11},
    "HfNiSn": {"ZT": 0.6, "T_opt_K": 900, "type": "n", "a": 6.08},
    "TiCoSb": {"ZT": 0.5, "T_opt_K": 800, "type": "n", "a": 5.88},
    "ZrCoSb": {"ZT": 0.4, "T_opt_K": 900, "type": "p", "a": 6.07},
    "TiNiSn_Sb": {"ZT": 1.5, "T_opt_K": 700, "type": "n", "dopant": "Sb",
                  "description": "Optimized half-Heusler"},
    
    # Skutterudites (mid temperature)
    "CoSb3": {"ZT": 0.1, "T_opt_K": 700, "type": "n", "a": 9.04,
              "structure": "skutterudite"},
    "CoSb3_filled_Yb": {"ZT": 1.4, "T_opt_K": 800, "type": "n", "filler": "Yb",
                         "description": "Yb-filled skutterudite"},
    "CoSb3_filled_Ba": {"ZT": 1.2, "T_opt_K": 850, "type": "n", "filler": "Ba"},
    "FeSb2Te": {"ZT": 0.8, "T_opt_K": 700, "type": "p"},
    
    # Clathrates
    "Ba8Ga16Ge30": {"ZT": 1.4, "T_opt_K": 900, "type": "n", "a": 10.78,
                    "structure": "clathrate_I"},
    "Ba8Ga16Si30": {"ZT": 0.9, "T_opt_K": 900, "type": "n"},
    
    # Oxides (high temperature)
    "NaCo2O4": {"ZT": 0.8, "T_opt_K": 800, "type": "p", 
                "structure": "layered", "oxide": True},
    "Ca3Co4O9": {"ZT": 0.8, "T_opt_K": 900, "type": "p", "misfit": True},
    "SrTiO3_Nb": {"ZT": 0.4, "T_opt_K": 1000, "type": "n", "dopant": "Nb"},
    "ZnO_Al": {"ZT": 0.3, "T_opt_K": 1000, "type": "n", "dopant": "Al"},
    
    # Silicon-based
    "SiGe": {"ZT": 1.0, "T_opt_K": 1200, "type": "both", 
             "description": "High temp, NASA RTG"},
    "Mg2Si": {"ZT": 0.9, "T_opt_K": 800, "type": "n", "a": 6.35},
    "Mg2Sn": {"ZT": 0.6, "T_opt_K": 700, "type": "n"},
    
    # Novel materials
    "SnSe": {"ZT": 2.6, "T_opt_K": 923, "type": "p", "a": 11.49, "b": 4.15, "c": 4.44,
             "structure": "orthorhombic", "record_ZT": True},
    "Cu2Se": {"ZT": 2.3, "T_opt_K": 1000, "type": "p", "superionic": True},
    "BiCuSeO": {"ZT": 1.4, "T_opt_K": 900, "type": "p", "layered": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_thermoelectric(
    material: str = "Bi2Te3",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate thermoelectric material structure.
    
    Args:
        material: Material from database
        supercell: Supercell dimensions
    
    Returns:
        Thermoelectric structure
    """
    if material not in THERMOELECTRIC_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material",
                      "available": list(THERMOELECTRIC_DATABASE.keys())}
        }
    
    info = THERMOELECTRIC_DATABASE[material]
    
    species = []
    coords = []
    
    structure_type = info.get("structure", "rocksalt")
    
    if structure_type == "rhombohedral" or material.startswith("Bi2"):
        a = info.get("a", 4.38)
        c = info.get("c", 30.49)
        lattice = Lattice.hexagonal(a, c)
        
        # Bi2Te3/Bi2Se3 structure
        if "Te" in material:
            species = ["Bi", "Bi", "Te", "Te", "Te"]
            coords = [[0, 0, 0.4], [0, 0, 0.6], [0, 0, 0], [0, 0, 0.21], [0, 0, 0.79]]
        else:
            species = ["Bi", "Bi", "Se", "Se", "Se"]
            coords = [[0, 0, 0.4], [0, 0, 0.6], [0, 0, 0], [0, 0, 0.21], [0, 0, 0.79]]
            
    elif structure_type == "rocksalt":
        a = info.get("a", 6.46)
        lattice = Lattice.cubic(a)
        
        # Extract elements from formula
        if "PbTe" in material:
            species = ["Pb", "Te"]
        elif "PbSe" in material:
            species = ["Pb", "Se"]
        elif "SnTe" in material:
            species = ["Sn", "Te"]
        else:
            species = ["Pb", "Te"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        
    elif structure_type == "half_heusler":
        a = info.get("a", 6.0)
        lattice = Lattice.cubic(a)
        species = ["Ti", "Ni", "Sn"]
        coords = [[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5]]
        
    elif structure_type == "skutterudite":
        a = info.get("a", 9.04)
        lattice = Lattice.cubic(a)
        # Simplified CoSb3
        species = ["Co", "Co", "Sb", "Sb", "Sb", "Sb", "Sb", "Sb"]
        coords = [
            [0.25, 0.25, 0.25], [0.75, 0.75, 0.75],
            [0, 0.335, 0.158], [0, 0.665, 0.842],
            [0.335, 0.158, 0], [0.665, 0.842, 0],
            [0.158, 0, 0.335], [0.842, 0, 0.665]
        ]
        
    elif structure_type == "clathrate_I":
        a = info.get("a", 10.78)
        lattice = Lattice.cubic(a)
        species = ["Ba", "Ba", "Ga", "Ge"]
        coords = [[0, 0, 0], [0.25, 0.5, 0], [0.18, 0.18, 0.18], [0.12, 0.31, 0]]
        
    elif structure_type == "orthorhombic" or "SnSe" in material:
        a = info.get("a", 11.49)
        b = info.get("b", 4.15)
        c = info.get("c", 4.44)
        lattice = Lattice.orthorhombic(a, b, c)
        species = ["Sn", "Sn", "Se", "Se"]
        coords = [[0.12, 0.25, 0.1], [0.88, 0.75, 0.9], 
                 [0.35, 0.25, 0.48], [0.65, 0.75, 0.52]]
    else:
        a = info.get("a", 5.0)
        lattice = Lattice.cubic(a)
        species = [material[:2]]
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "material": material,
        "ZT": info["ZT"],
        "T_optimal_K": info["T_opt_K"],
        "carrier_type": info["type"],
        "is_record_ZT": info.get("record_ZT", False),
        "is_Pb_free": info.get("Pb_free", False),
        "dopant": info.get("dopant", None),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def get_thermoelectric_database() -> Dict[str, Any]:
    """Get thermoelectrics organized by temperature range."""
    return {
        "success": True,
        "n_materials": len(THERMOELECTRIC_DATABASE),
        "by_temperature": {
            "low_T_300K": [k for k, v in THERMOELECTRIC_DATABASE.items() if v["T_opt_K"] <= 400],
            "mid_T_500_800K": [k for k, v in THERMOELECTRIC_DATABASE.items() if 400 < v["T_opt_K"] <= 800],
            "high_T_900K_plus": [k for k, v in THERMOELECTRIC_DATABASE.items() if v["T_opt_K"] > 800],
        },
        "high_ZT": [k for k, v in THERMOELECTRIC_DATABASE.items() if v["ZT"] >= 1.5]
    }
