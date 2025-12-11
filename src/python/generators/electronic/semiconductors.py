"""
electronic/semiconductors.py - Semiconductor Materials

Comprehensive semiconductor structure generation with accurate lattice parameters.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


SEMICONDUCTOR_DATABASE = {
    # Group IV
    "Si": {"a": 5.431, "structure": "diamond", "bandgap": 1.12, "type": "indirect", "group": "IV"},
    "Ge": {"a": 5.658, "structure": "diamond", "bandgap": 0.66, "type": "indirect", "group": "IV"},
    "C-diamond": {"a": 3.567, "structure": "diamond", "bandgap": 5.47, "type": "indirect", "group": "IV"},
    "SiC-3C": {"a": 4.359, "structure": "zincblende", "bandgap": 2.36, "type": "indirect", "group": "IV-IV"},
    "SiC-4H": {"a": 3.073, "c": 10.053, "structure": "wurtzite", "bandgap": 3.23, "type": "indirect", "group": "IV-IV"},
    "SiC-6H": {"a": 3.081, "c": 15.117, "structure": "wurtzite", "bandgap": 3.00, "type": "indirect", "group": "IV-IV"},
    
    # III-V Zincblende
    "GaAs": {"a": 5.653, "structure": "zincblende", "bandgap": 1.42, "type": "direct", "group": "III-V"},
    "GaP": {"a": 5.451, "structure": "zincblende", "bandgap": 2.26, "type": "indirect", "group": "III-V"},
    "GaSb": {"a": 6.096, "structure": "zincblende", "bandgap": 0.73, "type": "direct", "group": "III-V"},
    "InAs": {"a": 6.058, "structure": "zincblende", "bandgap": 0.36, "type": "direct", "group": "III-V"},
    "InP": {"a": 5.869, "structure": "zincblende", "bandgap": 1.35, "type": "direct", "group": "III-V"},
    "InSb": {"a": 6.479, "structure": "zincblende", "bandgap": 0.17, "type": "direct", "group": "III-V"},
    "AlAs": {"a": 5.661, "structure": "zincblende", "bandgap": 2.16, "type": "indirect", "group": "III-V"},
    "AlP": {"a": 5.467, "structure": "zincblende", "bandgap": 2.45, "type": "indirect", "group": "III-V"},
    "AlSb": {"a": 6.136, "structure": "zincblende", "bandgap": 1.62, "type": "indirect", "group": "III-V"},
    
    # III-V Wurtzite (Nitrides)
    "GaN": {"a": 3.189, "c": 5.185, "structure": "wurtzite", "bandgap": 3.44, "type": "direct", "group": "III-N"},
    "AlN": {"a": 3.112, "c": 4.982, "structure": "wurtzite", "bandgap": 6.2, "type": "direct", "group": "III-N"},
    "InN": {"a": 3.545, "c": 5.703, "structure": "wurtzite", "bandgap": 0.7, "type": "direct", "group": "III-N"},
    
    # II-VI
    "ZnO": {"a": 3.250, "c": 5.207, "structure": "wurtzite", "bandgap": 3.37, "type": "direct", "group": "II-VI"},
    "ZnS": {"a": 5.409, "structure": "zincblende", "bandgap": 3.68, "type": "direct", "group": "II-VI"},
    "ZnSe": {"a": 5.668, "structure": "zincblende", "bandgap": 2.70, "type": "direct", "group": "II-VI"},
    "ZnTe": {"a": 6.101, "structure": "zincblende", "bandgap": 2.26, "type": "direct", "group": "II-VI"},
    "CdS": {"a": 5.832, "structure": "zincblende", "bandgap": 2.42, "type": "direct", "group": "II-VI"},
    "CdSe": {"a": 6.052, "structure": "zincblende", "bandgap": 1.74, "type": "direct", "group": "II-VI"},
    "CdTe": {"a": 6.482, "structure": "zincblende", "bandgap": 1.49, "type": "direct", "group": "II-VI"},
    "HgTe": {"a": 6.461, "structure": "zincblende", "bandgap": -0.30, "type": "inverted", "group": "II-VI"},
    
    # Oxides
    "Cu2O": {"a": 4.269, "structure": "cubic", "bandgap": 2.17, "type": "direct", "group": "oxide"},
    "TiO2-rutile": {"a": 4.594, "c": 2.959, "structure": "tetragonal", "bandgap": 3.0, "type": "indirect", "group": "oxide"},
    "TiO2-anatase": {"a": 3.785, "c": 9.514, "structure": "tetragonal", "bandgap": 3.2, "type": "indirect", "group": "oxide"},
    "SnO2": {"a": 4.737, "c": 3.186, "structure": "tetragonal", "bandgap": 3.6, "type": "direct", "group": "oxide"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_semiconductor(
    material: str = "Si",
    size: List[int] = [1, 1, 1],
    doping: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    """
    Generate semiconductor structure.
    
    Args:
        material: Semiconductor name
        size: Supercell size
        doping: Doping configuration {"element": concentration}
    
    Returns:
        Semiconductor structure with properties
    """
    if material not in SEMICONDUCTOR_DATABASE:
        return {
            "success": False,
            "error": {"code": "UNKNOWN", "message": f"Unknown semiconductor: {material}",
                      "available": list(SEMICONDUCTOR_DATABASE.keys())}
        }
    
    info = SEMICONDUCTOR_DATABASE[material]
    a = info["a"]
    struct_type = info["structure"]
    
    if struct_type == "diamond":
        lattice = Lattice.cubic(a)
        element = material.split("-")[0] if "-" in material else material
        species = [element] * 8
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
    
    elif struct_type == "zincblende":
        lattice = Lattice.cubic(a)
        # Parse elements from formula
        if len(material) == 4:  # e.g., GaAs
            e1, e2 = material[:2], material[2:]
        else:
            e1, e2 = material[0], material[1:]
        species = [e1, e1, e1, e1, e2, e2, e2, e2]
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
    
    elif struct_type == "wurtzite":
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        if material.startswith("GaN") or material.startswith("AlN") or material.startswith("InN"):
            e1, e2 = material[:2], "N"
        elif material == "ZnO":
            e1, e2 = "Zn", "O"
        else:
            e1 = material[:2] if len(material) > 2 else material[0]
            e2 = material[2:] if len(material) > 2 else material[1:]
        species = [e1, e1, e2, e2]
        u = 0.375
        coords = [[1/3, 2/3, 0], [2/3, 1/3, 0.5], [1/3, 2/3, u], [2/3, 1/3, 0.5 + u]]
    
    elif struct_type == "tetragonal":
        c = info["c"]
        lattice = Lattice.tetragonal(a, c)
        if "TiO2" in material:
            species = ["Ti", "Ti", "O", "O", "O", "O"]
            coords = [[0, 0, 0], [0.5, 0.5, 0.5],
                      [0.305, 0.305, 0], [0.695, 0.695, 0], [0.195, 0.805, 0.5], [0.805, 0.195, 0.5]]
        else:
            species = ["Sn", "Sn", "O", "O", "O", "O"]
            coords = [[0, 0, 0], [0.5, 0.5, 0.5],
                      [0.307, 0.307, 0], [0.693, 0.693, 0], [0.193, 0.807, 0.5], [0.807, 0.193, 0.5]]
    
    else:
        lattice = Lattice.cubic(a)
        species = ["Cu", "Cu", "O"]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Apply doping
    if doping:
        for dopant, conc in doping.items():
            n_dopant = max(1, int(len(structure) * conc))
            for i in range(n_dopant):
                if i < len(structure):
                    structure[i] = dopant
    
    return {
        "success": True,
        "material": material,
        "structure_type": struct_type,
        "bandgap_eV": info["bandgap"],
        "gap_type": info["type"],
        "group": info["group"],
        "a_angstrom": a,
        "c_angstrom": info.get("c"),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
