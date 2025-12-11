"""
electronic/superconductors.py - Superconductor Materials

Conventional BCS, cuprates, iron-based, heavy fermion superconductors.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


SUPERCONDUCTOR_DATABASE = {
    # Conventional BCS
    "Nb": {"a": 3.301, "structure": "bcc", "tc_kelvin": 9.2, "type": "conventional"},
    "V": {"a": 3.024, "structure": "bcc", "tc_kelvin": 5.4, "type": "conventional"},
    "NbN": {"a": 4.392, "structure": "rocksalt", "tc_kelvin": 16.0, "type": "conventional"},
    "MgB2": {"a": 3.086, "c": 3.524, "structure": "hexagonal", "tc_kelvin": 39.0, "type": "conventional"},
    "Nb3Sn": {"a": 5.289, "structure": "A15", "tc_kelvin": 18.3, "type": "conventional"},
    
    # Cuprate High-Tc
    "YBCO": {"a": 3.817, "b": 3.884, "c": 11.68, "structure": "orthorhombic", 
             "tc_kelvin": 93, "type": "cuprate", "formula": "YBa2Cu3O7"},
    "LSCO": {"a": 3.78, "c": 13.2, "structure": "tetragonal",
             "tc_kelvin": 38, "type": "cuprate", "formula": "La1.85Sr0.15CuO4"},
    "Bi2212": {"a": 5.41, "b": 5.41, "c": 30.89, "structure": "tetragonal",
               "tc_kelvin": 85, "type": "cuprate", "formula": "Bi2Sr2CaCu2O8"},
    "Hg1223": {"a": 3.85, "c": 15.85, "structure": "tetragonal",
               "tc_kelvin": 135, "type": "cuprate", "formula": "HgBa2Ca2Cu3O8"},
    
    # Iron-based
    "LaFeAsO": {"a": 4.032, "c": 8.741, "structure": "tetragonal", 
                "tc_kelvin": 26, "type": "iron-based", "family": "1111"},
    "FeSe": {"a": 3.765, "c": 5.518, "structure": "tetragonal",
             "tc_kelvin": 8, "type": "iron-based", "family": "11"},
    "BaFe2As2": {"a": 3.963, "c": 13.017, "structure": "tetragonal",
                 "tc_kelvin": 38, "type": "iron-based", "family": "122"},
    
    # Heavy Fermion
    "CeCoIn5": {"a": 4.613, "c": 7.551, "structure": "tetragonal",
                "tc_kelvin": 2.3, "type": "heavy-fermion"},
    "UPt3": {"a": 5.764, "c": 4.884, "structure": "hexagonal",
             "tc_kelvin": 0.5, "type": "heavy-fermion"},
    
    # Room-temperature candidates
    "H3S": {"a": 3.08, "structure": "cubic", "tc_kelvin": 203, 
            "type": "hydride", "pressure_GPa": 155},
    "LaH10": {"a": 5.12, "structure": "fcc", "tc_kelvin": 250,
              "type": "hydride", "pressure_GPa": 170},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_superconductor(
    material: str = "MgB2",
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate superconductor structure.
    
    Args:
        material: Superconductor name
        size: Supercell size
    
    Returns:
        Superconductor structure with Tc
    """
    if material not in SUPERCONDUCTOR_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {material}",
                "available": list(SUPERCONDUCTOR_DATABASE.keys())}}
    
    info = SUPERCONDUCTOR_DATABASE[material]
    a = info["a"]
    struct_type = info["structure"]
    
    if struct_type == "bcc":
        lattice = Lattice.cubic(a)
        species = [material, material]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    
    elif struct_type == "hexagonal":
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        if material == "MgB2":
            species = ["Mg", "B", "B"]
            coords = [[0, 0, 0], [1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]
        else:
            species = ["U", "Pt", "Pt", "Pt"]
            coords = [[0, 0, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0]]
    
    elif struct_type == "A15":
        lattice = Lattice.cubic(a)
        species = ["Nb"] * 6 + ["Sn"] * 2
        coords = [[0.25, 0, 0.5], [0.75, 0, 0.5], [0.5, 0.25, 0], [0.5, 0.75, 0],
                  [0, 0.5, 0.25], [0, 0.5, 0.75], [0, 0, 0], [0.5, 0.5, 0.5]]
    
    elif struct_type == "rocksalt":
        lattice = Lattice.cubic(a)
        species = ["Nb", "Nb", "Nb", "Nb", "N", "N", "N", "N"]
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
    
    elif struct_type == "tetragonal":
        c = info["c"]
        lattice = Lattice.tetragonal(a, c)
        if material == "FeSe":
            species = ["Fe", "Fe", "Se", "Se"]
            coords = [[0, 0.5, 0], [0.5, 0, 0], [0, 0, 0.23], [0.5, 0.5, 0.77]]
        else:
            species = ["Ba", "Fe", "Fe", "As", "As"]
            coords = [[0, 0, 0], [0.5, 0, 0.25], [0, 0.5, 0.75], [0, 0, 0.35], [0.5, 0.5, 0.65]]
    
    elif struct_type == "orthorhombic":
        b = info["b"]
        c = info["c"]
        lattice = Lattice.orthorhombic(a, b, c)
        # Simplified YBCO
        species = ["Y", "Ba", "Ba", "Cu", "Cu", "Cu", "O", "O", "O", "O", "O", "O", "O"]
        coords = [[0.5, 0.5, 0.5], [0.5, 0.5, 0.18], [0.5, 0.5, 0.82],
                  [0, 0, 0], [0, 0, 0.35], [0, 0, 0.65]] + [[0.5, 0, i/7] for i in range(7)]
    
    else:
        lattice = Lattice.cubic(a)
        if material == "H3S":
            species = ["S", "H", "H", "H"]
            coords = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
        else:
            species = ["La"] + ["H"] * 10
            coords = [[0, 0, 0]] + [[i/10, i/10, i/10] for i in range(1, 11)]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": material,
        "formula": info.get("formula", material),
        "tc_kelvin": info["tc_kelvin"],
        "type": info["type"],
        "pressure_GPa": info.get("pressure_GPa"),
        "structure": structure_to_dict(structure)
    }
