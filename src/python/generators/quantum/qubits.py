"""
quantum/qubits.py - Qubit Host Materials

Color centers and qubit materials: NV centers, SiV, defects.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


QUBIT_DATABASE = {
    # Diamond Color Centers
    "NV-diamond": {"host": "diamond", "a": 3.567, "defect": "NV", "spin": 1,
                   "zpl_nm": 637, "t2_us": 1800, "description": "Nitrogen-vacancy center"},
    "SiV-diamond": {"host": "diamond", "a": 3.567, "defect": "SiV", "spin": 0.5,
                    "zpl_nm": 737, "t2_us": 13, "description": "Silicon-vacancy center"},
    "GeV-diamond": {"host": "diamond", "a": 3.567, "defect": "GeV", "spin": 0.5,
                    "zpl_nm": 602, "t2_us": 20, "description": "Germanium-vacancy center"},
    "SnV-diamond": {"host": "diamond", "a": 3.567, "defect": "SnV", "spin": 0.5,
                    "zpl_nm": 620, "t2_us": 100, "description": "Tin-vacancy center"},
    
    # SiC Color Centers
    "VSi-4H-SiC": {"host": "4H-SiC", "a": 3.073, "c": 10.053, "defect": "VSi", "spin": 1.5,
                   "zpl_nm": 862, "t2_us": 160, "description": "Silicon vacancy in 4H-SiC"},
    "divacancy-4H-SiC": {"host": "4H-SiC", "a": 3.073, "c": 10.053, "defect": "VSiVC", "spin": 1,
                         "zpl_nm": 1078, "t2_us": 1000, "description": "Divacancy in SiC"},
    
    # hBN Defects
    "VB-hBN": {"host": "hBN", "a": 2.50, "c": 6.66, "defect": "VB", "spin": 1,
               "zpl_nm": 850, "description": "Boron vacancy in hBN"},
    "VNCB-hBN": {"host": "hBN", "a": 2.50, "c": 6.66, "defect": "VNCB", "spin": 0.5,
                 "zpl_nm": 623, "description": "Complex defect in hBN"},
    
    # Rare Earth Qubits
    "Er-Si": {"host": "Si", "a": 5.431, "defect": "Er", "spin": 0.5,
              "zpl_nm": 1536, "description": "Erbium in silicon (telecom)"},
    "Eu-YSO": {"host": "YSO", "a": 10.41, "b": 6.72, "c": 12.49, "defect": "Eu",
               "spin": 0, "t2_us": 6000, "description": "Europium in Y2SiO5"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    return {
        "lattice": {"a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                    "matrix": structure.lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
    }


def generate_qubit_host(
    qubit_type: str = "NV-diamond",
    size: List[int] = [3, 3, 3],
    include_defect: bool = True
) -> Dict[str, Any]:
    """
    Generate qubit host structure with color center.
    
    Args:
        qubit_type: Qubit/defect type
        size: Supercell size
        include_defect: Whether to include the defect
    
    Returns:
        Qubit host structure
    """
    if qubit_type not in QUBIT_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {qubit_type}",
                "available": list(QUBIT_DATABASE.keys())}}
    
    info = QUBIT_DATABASE[qubit_type]
    a = info["a"]
    host = info["host"]
    
    if host == "diamond":
        lattice = Lattice.cubic(a)
        species = ["C"] * 8
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
    
    elif "SiC" in host:
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        species = ["Si", "Si", "C", "C"]
        coords = [[0, 0, 0], [1/3, 2/3, 0.5], [0, 0, 0.125], [1/3, 2/3, 0.625]]
    
    elif host == "hBN":
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        species = ["B", "N"]
        coords = [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]
    
    else:  # Si or YSO
        lattice = Lattice.cubic(a)
        species = ["Si"] * 8
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
    
    structure = Structure(lattice, species, coords)
    structure.make_supercell(size)
    
    defect_position = None
    if include_defect:
        center = len(structure) // 2
        defect = info["defect"]
        
        if defect == "NV":
            # Remove C, add N at adjacent site
            structure.remove_sites([center])
            structure[center - 1] = "N"
            defect_position = "center"
        elif defect == "SiV":
            structure.remove_sites([center])
            structure.append("Si", [0.5, 0.5, 0.5])
            defect_position = "center"
        elif "V" in defect:
            structure.remove_sites([center])
            defect_position = "center"
    
    return {
        "success": True,
        "qubit_type": qubit_type,
        "host": host,
        "defect": info["defect"],
        "spin": info["spin"],
        "zpl_nm": info.get("zpl_nm"),
        "coherence_time_us": info.get("t2_us"),
        "description": info["description"],
        "defect_included": include_defect,
        "structure": structure_to_dict(structure)
    }
