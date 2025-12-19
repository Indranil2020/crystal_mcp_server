"""
electronic/superconductors.py - Superconducting Materials

Comprehensive superconductor generation per structure_catalogue.md:
- Conventional (BCS) superconductors
- High-Tc cuprates
- Iron-based superconductors
- Heavy fermion superconductors
- Topological superconductors
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Superconductor database
SUPERCONDUCTOR_DATABASE = {
    # Conventional BCS
    "Nb": {"Tc_K": 9.3, "type": "elemental", "a": 3.30, "structure": "bcc", "gap_meV": 1.55},
    "Pb": {"Tc_K": 7.2, "type": "elemental", "a": 4.95, "structure": "fcc", "gap_meV": 1.35},
    "NbN": {"Tc_K": 16.0, "type": "binary", "a": 4.39, "structure": "rocksalt", "gap_meV": 2.4},
    "NbTi": {"Tc_K": 10.0, "type": "alloy", "a": 3.28, "structure": "bcc"},
    "Nb3Sn": {"Tc_K": 18.3, "type": "A15", "a": 5.29, "structure": "A15", "gap_meV": 3.2},
    "Nb3Ge": {"Tc_K": 23.2, "type": "A15", "a": 5.14, "structure": "A15"},
    "MgB2": {"Tc_K": 39.0, "type": "two_gap", "a": 3.09, "c": 3.52, "structure": "hexagonal",
             "gap1_meV": 7.0, "gap2_meV": 2.3},
    
    # High-Tc Cuprates
    "YBCO": {"Tc_K": 93, "type": "cuprate", "formula": "YBa2Cu3O7", 
             "a": 3.82, "b": 3.89, "c": 11.68, "structure": "orthorhombic",
             "CuO_planes": 2, "gap_meV": 30},
    "Bi2212": {"Tc_K": 85, "type": "cuprate", "formula": "Bi2Sr2CaCu2O8",
               "a": 5.41, "c": 30.89, "structure": "tetragonal", "CuO_planes": 2},
    "Bi2223": {"Tc_K": 110, "type": "cuprate", "formula": "Bi2Sr2Ca2Cu3O10",
               "CuO_planes": 3},
    "Tl2201": {"Tc_K": 90, "type": "cuprate", "formula": "Tl2Ba2CuO6",
               "CuO_planes": 1},
    "Hg1201": {"Tc_K": 94, "type": "cuprate", "formula": "HgBa2CuO4",
               "CuO_planes": 1},
    "Hg1223": {"Tc_K": 135, "type": "cuprate", "formula": "HgBa2Ca2Cu3O8",
               "CuO_planes": 3, "record_Tc": True},
    "LSCO": {"Tc_K": 40, "type": "cuprate", "formula": "La1.85Sr0.15CuO4",
             "a": 3.78, "c": 13.23, "structure": "tetragonal", "CuO_planes": 1},
    
    # Iron-based
    "LaFeAsO_F": {"Tc_K": 26, "type": "iron_1111", "formula": "LaFeAsO0.9F0.1",
                  "a": 4.03, "c": 8.74, "structure": "tetragonal"},
    "SmFeAsO_F": {"Tc_K": 55, "type": "iron_1111", "formula": "SmFeAsO0.85F0.15"},
    "BaFe2As2_K": {"Tc_K": 38, "type": "iron_122", "formula": "Ba0.6K0.4Fe2As2",
                   "a": 3.94, "c": 13.30},
    "FeSe": {"Tc_K": 8, "type": "iron_11", "a": 3.77, "c": 5.52,
             "monolayer_Tc": 65},  # Enhanced in thin films
    "FeTe0.5Se0.5": {"Tc_K": 14, "type": "iron_11"},
    "LiFeAs": {"Tc_K": 18, "type": "iron_111", "a": 3.79, "c": 6.36},
    
    # Heavy fermion
    "CeCoIn5": {"Tc_K": 2.3, "type": "heavy_fermion", "a": 4.62, "c": 7.51,
                "d_wave": True},
    "UPt3": {"Tc_K": 0.5, "type": "heavy_fermion", "unconventional": True},
    "CeCu2Si2": {"Tc_K": 0.6, "type": "heavy_fermion", "first_HF": True},
    
    # Topological
    "Cu_Bi2Se3": {"Tc_K": 3.8, "type": "topological", "formula": "CuxBi2Se3",
                  "p_wave": True, "Majorana": True},
    "Sn_TaS2": {"Tc_K": 3.0, "type": "topological", "formula": "SnxTaS2"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_superconductor(
    material: str = "YBCO",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate superconductor structure.
    
    Args:
        material: Superconductor from database
        supercell: Supercell dimensions
    
    Returns:
        Superconductor structure
    """
    if material not in SUPERCONDUCTOR_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material '{material}'",
                      "available": list(SUPERCONDUCTOR_DATABASE.keys())}
        }
    
    info = SUPERCONDUCTOR_DATABASE[material]
    
    species = []
    coords = []
    
    # Generate structure based on type
    if info.get("structure") == "bcc":
        a = info["a"]
        lattice = Lattice.cubic(a)
        species = ["Nb", "Nb"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        
    elif info.get("structure") == "fcc":
        a = info["a"]
        lattice = Lattice.cubic(a)
        species = ["Pb"] * 4
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        
    elif info.get("structure") == "A15":
        a = info["a"]
        lattice = Lattice.cubic(a)
        # A15 structure (Nb3Sn type)
        species = ["Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "Sn", "Sn"]
        coords = [
            [0.25, 0, 0.5], [0.75, 0, 0.5], [0.5, 0.25, 0], [0.5, 0.75, 0],
            [0, 0.5, 0.25], [0, 0.5, 0.75], [0, 0, 0], [0.5, 0.5, 0.5]
        ]
        
    elif info.get("structure") == "hexagonal" and "MgB2" in material:
        a = info["a"]
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        species = ["Mg", "B", "B"]
        coords = [[0, 0, 0], [1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]
        
    elif info["type"] == "cuprate":
        # YBCO structure (simplified)
        a = info.get("a", 3.85)
        b = info.get("b", a)
        c = info.get("c", 11.68)
        lattice = Lattice.orthorhombic(a, b, c)
        
        # Layer structure: CuO chains, BaO, CuO2 planes, Y, CuO2, BaO, CuO chains
        species = ["Y", "Ba", "Ba", "Cu", "Cu", "Cu", "O", "O", "O", "O", "O", "O", "O"]
        coords = [
            [0.5, 0.5, 0.5],  # Y
            [0.5, 0.5, 0.18], [0.5, 0.5, 0.82],  # Ba
            [0, 0, 0], [0, 0, 0.35], [0, 0, 0.65],  # Cu
            [0, 0.5, 0], [0.5, 0, 0],  # Chain O
            [0, 0.5, 0.38], [0.5, 0, 0.38],  # Plane O
            [0, 0.5, 0.62], [0.5, 0, 0.62],  # Plane O
            [0, 0, 0.16]  # Apical O
        ]
        
    elif info["type"] in ["iron_1111", "iron_122", "iron_11"]:
        a = info.get("a", 4.0)
        c = info.get("c", 8.5)
        lattice = Lattice.tetragonal(a, c)
        
        if "FeSe" in material or info["type"] == "iron_11":
            species = ["Fe", "Fe", "Se", "Se"]
            coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.25], [0.5, 0, 0.75]]
        else:
            # 122 structure
            species = ["Ba", "Fe", "Fe", "As", "As"]
            coords = [[0, 0, 0], [0.5, 0, 0.25], [0.5, 0, 0.75], 
                     [0, 0.5, 0.35], [0, 0.5, 0.65]]
    else:
        a = info.get("a", 4.0)
        lattice = Lattice.cubic(a)
        species = [material[:2]]
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "material": material,
        "Tc_K": info["Tc_K"],
        "type": info["type"],
        "gap_meV": info.get("gap_meV", 0),
        "CuO_planes": info.get("CuO_planes", 0),
        "is_unconventional": info.get("unconventional", False) or info["type"] in ["cuprate", "heavy_fermion"],
        "has_Majorana": info.get("Majorana", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_cuprate_superconductor(
    base: str = "YBCO",
    doping: float = 0.15,
    n_CuO_planes: int = 2,
    supercell: List[int] = [2, 2, 1]
) -> Dict[str, Any]:
    """
    Generate cuprate superconductor with doping control.
    
    Args:
        base: Base cuprate type
        doping: Hole doping level
        n_CuO_planes: Number of CuO2 planes
        supercell: Supercell dimensions
    
    Returns:
        Cuprate structure
    """
    result = generate_superconductor(base, supercell)
    
    if result["success"]:
        result["doping"] = doping
        result["optimal_doping"] = 0.16
        result["underdoped"] = doping < 0.16
        result["overdoped"] = doping > 0.16
        
        # Estimate Tc based on doping (dome shape)
        Tc_max = SUPERCONDUCTOR_DATABASE[base]["Tc_K"]
        result["estimated_Tc_K"] = round(Tc_max * (1 - 82.6 * (doping - 0.16)**2), 1)
    
    return result


def get_superconductor_database() -> Dict[str, Any]:
    """Get all superconductors organized by type."""
    return {
        "success": True,
        "n_materials": len(SUPERCONDUCTOR_DATABASE),
        "by_type": {
            "conventional": [k for k, v in SUPERCONDUCTOR_DATABASE.items() 
                           if v["type"] in ["elemental", "binary", "alloy", "A15", "two_gap"]],
            "cuprate": [k for k, v in SUPERCONDUCTOR_DATABASE.items() if v["type"] == "cuprate"],
            "iron_based": [k for k, v in SUPERCONDUCTOR_DATABASE.items() if "iron" in v["type"]],
            "heavy_fermion": [k for k, v in SUPERCONDUCTOR_DATABASE.items() if v["type"] == "heavy_fermion"],
            "topological": [k for k, v in SUPERCONDUCTOR_DATABASE.items() if v["type"] == "topological"],
        }
    }
