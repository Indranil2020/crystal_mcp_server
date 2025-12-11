"""
high_pressure/phases.py - High Pressure Phases

Generates high-pressure structure phases:
- Post-perovskite MgSiO3
- Polymeric nitrogen (cg-N, LP-N)
- Super-ionic phases (ice-X)
- Superhydrides (LaH10, YH9)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# High-pressure phase database
HP_PHASE_DATABASE = {
    # Post-perovskite phases
    "MgSiO3-ppv": {
        "formula": "MgSiO3", "phase": "post-perovskite", "pressure_GPa": 125,
        "a": 2.456, "b": 8.042, "c": 6.093, "structure": "CaIrO3-type"
    },
    "CaIrO3-ppv": {
        "formula": "CaIrO3", "phase": "post-perovskite", "pressure_GPa": 0,
        "a": 3.145, "b": 9.865, "c": 7.297, "structure": "prototype"
    },
    
    # High-pressure CaCO3
    "CaCO3-aragonite": {
        "formula": "CaCO3", "phase": "aragonite", "pressure_GPa": 0.3,
        "a": 4.961, "b": 7.967, "c": 5.741, "structure": "orthorhombic"
    },
    "CaCO3-post-aragonite": {
        "formula": "CaCO3", "phase": "post-aragonite", "pressure_GPa": 40,
        "a": 4.21, "b": 4.21, "c": 4.21, "structure": "cubic"
    },
    
    # Black phosphorus phases
    "P-black": {
        "formula": "P", "phase": "black", "pressure_GPa": 0,
        "a": 3.314, "b": 4.376, "c": 10.478, "structure": "orthorhombic"
    },
    "P-blue": {
        "formula": "P", "phase": "blue", "pressure_GPa": 5,
        "a": 3.28, "b": 3.28, "c": 5.24, "structure": "hexagonal"
    },
    "P-simple-cubic": {
        "formula": "P", "phase": "simple-cubic", "pressure_GPa": 10,
        "a": 2.38, "structure": "cubic"
    },
}


# Superhydride database
SUPERHYDRIDE_DATABASE = {
    "LaH10": {
        "metal": "La", "n_H": 10, "structure": "Fm-3m", 
        "a": 5.12, "Tc": 250, "pressure_GPa": 170,
        "description": "Room-temperature superconductor candidate"
    },
    "YH9": {
        "metal": "Y", "n_H": 9, "structure": "P63/mmc",
        "a": 3.65, "c": 5.40, "Tc": 243, "pressure_GPa": 201,
        "description": "YBHO structure"
    },
    "YH6": {
        "metal": "Y", "n_H": 6, "structure": "Im-3m",
        "a": 3.54, "Tc": 220, "pressure_GPa": 166,
        "description": "Clathrate H cage"
    },
    "CaH6": {
        "metal": "Ca", "n_H": 6, "structure": "Im-3m",
        "a": 3.43, "Tc": 215, "pressure_GPa": 172,
        "description": "Sodalite-like cage"
    },
    "H3S": {
        "metal": "S", "n_H": 3, "structure": "Im-3m",
        "a": 3.08, "Tc": 203, "pressure_GPa": 155,
        "description": "First high-Tc superhydride"
    },
    "ThH10": {
        "metal": "Th", "n_H": 10, "structure": "Fm-3m",
        "a": 5.28, "Tc": 161, "pressure_GPa": 175,
        "description": "Heavy element superhydride"
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_high_pressure_phase(
    phase: str = "MgSiO3-ppv",
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate high-pressure phase structure.
    
    Args:
        phase: Phase name from database
        supercell: Supercell size
    
    Returns:
        High-pressure phase structure
    """
    if phase not in HP_PHASE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_PHASE", "message": f"Unknown phase '{phase}'",
                      "available": list(HP_PHASE_DATABASE.keys())}
        }
    
    info = HP_PHASE_DATABASE[phase]
    
    a = info["a"]
    b = info.get("b", a)
    c = info.get("c", a)
    
    if info["structure"] == "cubic":
        lattice = Lattice.cubic(a)
    elif info["structure"] == "hexagonal":
        lattice = Lattice.hexagonal(a, c)
    else:
        lattice = Lattice.orthorhombic(a, b, c)
    
    # Generate simplified structure
    formula = info["formula"]
    
    if "MgSiO3" in formula:
        # Post-perovskite structure
        species = ["Mg", "Si", "O", "O", "O"]
        coords = [[0, 0, 0.25], [0.5, 0.5, 0.5], [0, 0.5, 0.25], [0.5, 0, 0.75], [0.25, 0.25, 0]]
    elif "CaCO3" in formula:
        species = ["Ca", "C", "O", "O", "O"]
        coords = [[0, 0, 0], [0.25, 0.25, 0.25], [0.3, 0.3, 0.1], [0.1, 0.3, 0.3], [0.3, 0.1, 0.3]]
    elif formula == "P":
        species = ["P"] * 4
        coords = [[0, 0, 0], [0, 0.5, 0.25], [0.5, 0, 0.5], [0.5, 0.5, 0.75]]
    else:
        species = ["X"]
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "phase": phase,
        "formula": formula,
        "phase_name": info["phase"],
        "pressure_GPa": info["pressure_GPa"],
        "structure_type": info["structure"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_polymeric_nitrogen(
    phase: str = "cg-N",
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate polymeric nitrogen structure.
    
    Args:
        phase: 'cg-N' (cubic gauche) or 'LP-N' (layered polymeric)
        supercell: Supercell size
    
    Returns:
        Polymeric nitrogen structure
    """
    nitrogen_phases = {
        "cg-N": {"a": 3.45, "structure": "cubic", "pressure_GPa": 110, "coordination": 3},
        "LP-N": {"a": 4.07, "c": 4.30, "structure": "hexagonal", "pressure_GPa": 120, "coordination": 2},
        "bp-N": {"a": 3.80, "structure": "cubic", "pressure_GPa": 150, "coordination": 4},
    }
    
    if phase not in nitrogen_phases:
        return {"success": False, "error": {"code": "INVALID_PHASE", "message": f"Unknown phase '{phase}'"}}
    
    info = nitrogen_phases[phase]
    
    if info["structure"] == "cubic":
        lattice = Lattice.cubic(info["a"])
        # Cubic gauche structure (I213)
        species = ["N"] * 8
        coords = [
            [0.069, 0.069, 0.069], [0.431, 0.931, 0.569],
            [0.931, 0.569, 0.431], [0.569, 0.431, 0.931],
            [0.931, 0.431, 0.069], [0.069, 0.569, 0.931],
            [0.431, 0.069, 0.569], [0.569, 0.931, 0.431]
        ]
    else:
        lattice = Lattice.hexagonal(info["a"], info["c"])
        species = ["N"] * 4
        coords = [[0, 0, 0], [0.33, 0.67, 0.25], [0.67, 0.33, 0.5], [0, 0, 0.75]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "phase": phase,
        "pressure_GPa": info["pressure_GPa"],
        "coordination": info["coordination"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_super_ionic_phase(
    compound: str = "H2O",
    phase: str = "ice-XVIII",
    temperature_K: float = 5000
) -> Dict[str, Any]:
    """
    Generate super-ionic phase structure.
    
    Args:
        compound: Compound (H2O, NaCl-H2O, Li3N-LiOH)
        phase: Phase name
        temperature_K: Approximate temperature
    
    Returns:
        Super-ionic structure
    """
    super_ionic_phases = {
        "ice-X": {"a": 2.54, "pressure_GPa": 70, "mobile_ion": "H"},
        "ice-XVIII": {"a": 2.62, "pressure_GPa": 200, "mobile_ion": "H"},
        "NaCl-H2O": {"a": 5.64, "pressure_GPa": 50, "mobile_ion": "H"},
    }
    
    if phase not in super_ionic_phases:
        return {"success": False, "error": {"code": "INVALID_PHASE", "message": f"Unknown phase"}}
    
    info = super_ionic_phases[phase]
    
    lattice = Lattice.cubic(info["a"])
    
    if "ice" in phase:
        # Ice structure with random H positions
        species = ["O"] * 2
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        # Add delocalized H (represented as multiple fractional positions)
        for _ in range(4):
            np.random.seed(42)
            h_pos = np.random.uniform(0, 1, 3).tolist()
            species.append("H")
            coords.append(h_pos)
    else:
        species = ["Na", "Cl", "O", "H", "H"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.3, 0.3, 0.3], [0.2, 0.2, 0.2]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "compound": compound,
        "phase": phase,
        "pressure_GPa": info["pressure_GPa"],
        "temperature_K": temperature_K,
        "mobile_ion": info["mobile_ion"],
        "is_superionic": True,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_superhydride(
    compound: str = "LaH10",
    pressure_GPa: Optional[float] = None,
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate superhydride structure.
    
    Args:
        compound: Superhydride formula
        pressure_GPa: Override pressure
        supercell: Supercell size
    
    Returns:
        Superhydride structure
    """
    if compound not in SUPERHYDRIDE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_COMPOUND", "message": f"Unknown superhydride '{compound}'",
                      "available": list(SUPERHYDRIDE_DATABASE.keys())}
        }
    
    info = SUPERHYDRIDE_DATABASE[compound]
    
    a = info["a"]
    c = info.get("c", a)
    
    if info["structure"] == "P63/mmc":
        lattice = Lattice.hexagonal(a, c)
    else:
        lattice = Lattice.cubic(a)
    
    metal = info["metal"]
    n_H = info["n_H"]
    
    # Generate H cage structure
    species = [metal]
    coords = [[0, 0, 0]]
    
    # Add H atoms forming cage
    if n_H == 10:
        # LaH10 - H32 cage
        for i in range(n_H):
            theta = 2 * np.pi * i / n_H
            phi = np.pi * ((i % 3) + 1) / 4
            r = 0.3
            x = 0.5 + r * np.sin(phi) * np.cos(theta)
            y = 0.5 + r * np.sin(phi) * np.sin(theta)
            z = 0.5 + r * np.cos(phi)
            species.append("H")
            coords.append([x % 1, y % 1, z % 1])
    elif n_H == 6:
        # Clathrate H6 cage
        cage_positions = [
            [0.25, 0.25, 0], [0.75, 0.25, 0], [0.25, 0.75, 0],
            [0.25, 0, 0.25], [0.75, 0, 0.25], [0, 0.25, 0.25]
        ]
        for pos in cage_positions:
            species.append("H")
            coords.append(pos)
    else:
        # Generic
        for i in range(n_H):
            angle = 2 * np.pi * i / n_H
            species.append("H")
            coords.append([0.5 + 0.2 * np.cos(angle), 0.5 + 0.2 * np.sin(angle), 0.5])
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "compound": compound,
        "description": info["description"],
        "metal": metal,
        "n_hydrogen": n_H,
        "space_group": info["structure"],
        "superconducting_Tc_K": info["Tc"],
        "synthesis_pressure_GPa": pressure_GPa or info["pressure_GPa"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
