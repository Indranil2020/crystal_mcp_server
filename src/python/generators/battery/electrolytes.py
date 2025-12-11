"""
battery/electrolytes.py - Solid State Electrolytes

Comprehensive solid electrolyte generation:
- Garnet (LLZO)
- Sulfide (LGPS, Li6PS5Cl)
- NASICON (LATP, LAGP)
- Perovskite (LLTO)
- Polymer electrolytes
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Comprehensive electrolyte database
ELECTROLYTE_DATABASE = {
    # Garnets
    "LLZO": {
        "formula": "Li7La3Zr2O12", "structure": "garnet", "spacegroup": "Ia-3d",
        "a": 12.96, "conductivity": 1e-4, "activation_eV": 0.3,
        "stability_window": [0, 5.5], "elements": {"Li": 7, "La": 3, "Zr": 2, "O": 12}
    },
    "LLZO_Ta": {
        "formula": "Li6.5La3Zr1.5Ta0.5O12", "structure": "garnet", "spacegroup": "Ia-3d",
        "a": 12.94, "conductivity": 5e-4, "activation_eV": 0.28,
        "stability_window": [0, 5.5], "elements": {"Li": 6.5, "La": 3, "Zr": 1.5, "Ta": 0.5, "O": 12}
    },
    "LLZO_Al": {
        "formula": "Li6.25Al0.25La3Zr2O12", "structure": "garnet", "spacegroup": "Ia-3d",
        "a": 12.93, "conductivity": 3e-4, "activation_eV": 0.29,
        "stability_window": [0, 5.5], "elements": {"Li": 6.25, "Al": 0.25, "La": 3, "Zr": 2, "O": 12}
    },
    
    # Sulfides
    "LGPS": {
        "formula": "Li10GeP2S12", "structure": "LGPS", "spacegroup": "P42/nmc",
        "a": 8.72, "c": 12.63, "conductivity": 1.2e-2, "activation_eV": 0.22,
        "stability_window": [1.5, 2.5], "elements": {"Li": 10, "Ge": 1, "P": 2, "S": 12}
    },
    "LGPS_Si": {
        "formula": "Li10SiP2S12", "structure": "LGPS", "spacegroup": "P42/nmc",
        "a": 8.68, "c": 12.55, "conductivity": 2.3e-3, "activation_eV": 0.24,
        "stability_window": [1.5, 2.5], "elements": {"Li": 10, "Si": 1, "P": 2, "S": 12}
    },
    "Li6PS5Cl": {
        "formula": "Li6PS5Cl", "structure": "argyrodite", "spacegroup": "F-43m",
        "a": 9.86, "conductivity": 3e-3, "activation_eV": 0.25,
        "stability_window": [1.5, 2.8], "elements": {"Li": 6, "P": 1, "S": 5, "Cl": 1}
    },
    "Li6PS5Br": {
        "formula": "Li6PS5Br", "structure": "argyrodite", "spacegroup": "F-43m",
        "a": 9.98, "conductivity": 1e-3, "activation_eV": 0.27,
        "stability_window": [1.5, 2.8], "elements": {"Li": 6, "P": 1, "S": 5, "Br": 1}
    },
    "Li3PS4": {
        "formula": "Li3PS4", "structure": "orthorhombic", "spacegroup": "Pnma",
        "a": 12.82, "b": 8.13, "c": 6.11, "conductivity": 2e-4, "activation_eV": 0.35,
        "stability_window": [1.5, 2.5], "elements": {"Li": 3, "P": 1, "S": 4}
    },
    
    # NASICON
    "LATP": {
        "formula": "Li1.3Al0.3Ti1.7(PO4)3", "structure": "NASICON", "spacegroup": "R-3c",
        "a": 8.51, "c": 20.88, "conductivity": 1e-3, "activation_eV": 0.32,
        "stability_window": [2.0, 4.5], "elements": {"Li": 1.3, "Al": 0.3, "Ti": 1.7, "P": 3, "O": 12}
    },
    "LAGP": {
        "formula": "Li1.5Al0.5Ge1.5(PO4)3", "structure": "NASICON", "spacegroup": "R-3c",
        "a": 8.25, "c": 20.67, "conductivity": 4e-4, "activation_eV": 0.33,
        "stability_window": [2.0, 5.0], "elements": {"Li": 1.5, "Al": 0.5, "Ge": 1.5, "P": 3, "O": 12}
    },
    "LiZr2(PO4)3": {
        "formula": "LiZr2(PO4)3", "structure": "NASICON", "spacegroup": "R-3c",
        "a": 8.86, "c": 22.24, "conductivity": 2e-6, "activation_eV": 0.45,
        "stability_window": [0, 4.5], "elements": {"Li": 1, "Zr": 2, "P": 3, "O": 12}
    },
    
    # Perovskite
    "LLTO": {
        "formula": "Li0.33La0.56TiO3", "structure": "perovskite", "spacegroup": "P4/mmm",
        "a": 3.87, "c": 7.75, "conductivity": 1e-3, "activation_eV": 0.34,
        "stability_window": [2.5, 4.2], "elements": {"Li": 0.33, "La": 0.56, "Ti": 1, "O": 3}
    },
    
    # Thiophosphates
    "Li7P3S11": {
        "formula": "Li7P3S11", "structure": "triclinic", "spacegroup": "P-1",
        "a": 12.50, "b": 6.03, "c": 12.53, "conductivity": 1.7e-2, "activation_eV": 0.18,
        "stability_window": [1.5, 2.1], "elements": {"Li": 7, "P": 3, "S": 11}
    },
    
    # Sodium electrolytes
    "Na3PS4": {
        "formula": "Na3PS4", "structure": "cubic", "spacegroup": "I-43m",
        "a": 6.99, "conductivity": 2e-4, "activation_eV": 0.30,
        "stability_window": [1.0, 2.5], "elements": {"Na": 3, "P": 1, "S": 4}, "ion": "Na"
    },
    "NASICON_Na": {
        "formula": "Na3Zr2Si2PO12", "structure": "NASICON", "spacegroup": "C2/c",
        "a": 15.63, "c": 9.22, "conductivity": 1e-3, "activation_eV": 0.29,
        "stability_window": [0, 5.0], "elements": {"Na": 3, "Zr": 2, "Si": 2, "P": 1, "O": 12}, "ion": "Na"
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


def generate_electrolyte(
    electrolyte_type: str = "LLZO",
    supercell: List[int] = [1, 1, 1],
    li_site_disorder: bool = False
) -> Dict[str, Any]:
    """
    Generate solid electrolyte structure.
    
    Args:
        electrolyte_type: Electrolyte type from database
        supercell: Supercell dimensions
        li_site_disorder: Apply Li site disorder
    
    Returns:
        Electrolyte structure
    """
    if electrolyte_type not in ELECTROLYTE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ELECTROLYTE", "message": f"Unknown electrolyte '{electrolyte_type}'",
                      "available": list(ELECTROLYTE_DATABASE.keys())}
        }
    
    info = ELECTROLYTE_DATABASE[electrolyte_type]
    structure_type = info["structure"]
    
    if structure_type == "garnet":
        a = info["a"]
        lattice = Lattice.cubic(a)
        
        # Simplified garnet structure
        species = ["Li"] * 7 + ["La"] * 3 + ["Zr"] * 2 + ["O"] * 12
        
        # Approximate positions
        coords = []
        # Li sites (tetrahedral and octahedral)
        for i in range(7):
            coords.append([i * 0.1, 0.5 + i * 0.05, 0.25 + i * 0.1])
        # La sites
        for i in range(3):
            coords.append([0.125 + i * 0.25, 0, 0.25])
        # Zr sites
        coords.extend([[0, 0, 0], [0.5, 0.5, 0.5]])
        # O sites
        for i in range(12):
            coords.append([np.sin(i * np.pi/6) * 0.3 + 0.5, np.cos(i * np.pi/6) * 0.3 + 0.5, i / 12])
        
    elif structure_type == "LGPS":
        a = info["a"]
        c = info["c"]
        lattice = Lattice.tetragonal(a, c)
        
        # Simplified LGPS
        species = ["Li"] * 10 + ["Ge"] + ["P"] * 2 + ["S"] * 12
        coords = []
        
        for i in range(10):
            coords.append([0.1 * i, 0.5, i * 0.1])
        coords.append([0, 0, 0])  # Ge
        coords.extend([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])  # P
        for i in range(12):
            coords.append([i / 12, i / 12, 0.5 + i * 0.04])
        
    elif structure_type == "argyrodite":
        a = info["a"]
        lattice = Lattice.cubic(a)
        
        species = ["Li"] * 6 + ["P"] + ["S"] * 5 + ["Cl"]
        coords = [
            [0.12, 0, 0], [0.25, 0.25, 0.25], [0.38, 0, 0],
            [0, 0.12, 0], [0, 0.25, 0.25], [0, 0.38, 0],
            [0, 0, 0],
            [0.12, 0.12, 0.12], [0.38, 0.38, 0.12], [0.12, 0.38, 0.38],
            [0.38, 0.12, 0.38], [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.75]
        ]
        
    elif structure_type == "NASICON":
        a = info["a"]
        c = info.get("c", a * 2.5)
        lattice = Lattice.rhombohedral(a, 60)
        
        species = ["Li", "Al", "Ti", "Ti", "P", "P", "P"] + ["O"] * 12
        coords = [[0.1 * i, 0.1 * i, i * 0.05] for i in range(len(species))]
        
    elif structure_type == "perovskite":
        a = info["a"]
        c = info.get("c", a * 2)
        lattice = Lattice.tetragonal(a, c)
        
        # A-site deficient perovskite
        species = ["Li", "La", "La", "Ti", "O", "O", "O"]
        coords = [
            [0, 0, 0], [0, 0, 0.5], [0.5, 0.5, 0.25],
            [0.5, 0.5, 0], [0.5, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0.5]
        ]
    
    else:
        a = info.get("a", 10)
        lattice = Lattice.cubic(a)
        species = ["Li", "P", "S", "S"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    # Apply Li disorder if requested
    if li_site_disorder:
        li_indices = [i for i, s in enumerate(structure) if str(s.specie) == "Li"]
        np.random.seed(42)
        for idx in li_indices:
            site = structure[idx]
            # Small random displacement
            new_coords = site.frac_coords + np.random.uniform(-0.05, 0.05, 3)
            new_coords = new_coords % 1
            structure[idx] = "Li", new_coords
    
    stability_low, stability_high = info["stability_window"]
    
    return {
        "success": True,
        "electrolyte_type": electrolyte_type,
        "formula": info["formula"],
        "structure_type": structure_type,
        "spacegroup": info["spacegroup"],
        "ionic_conductivity_S_cm": info["conductivity"],
        "activation_energy_eV": info["activation_eV"],
        "stability_window_V": info["stability_window"],
        "mobile_ion": info.get("ion", "Li"),
        "li_site_disorder": li_site_disorder,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_cathode_electrolyte_interface(
    cathode: str = "NCM811",
    electrolyte: str = "LLZO",
    interface_orientation: str = "(001)"
) -> Dict[str, Any]:
    """
    Generate cathode-electrolyte interface.
    
    Args:
        cathode: Cathode type
        electrolyte: Electrolyte type
        interface_orientation: Interface plane
    
    Returns:
        Interface structure
    """
    # Import cathode generator
    from .cathodes import generate_cathode, CATHODE_DATABASE
    
    if cathode not in CATHODE_DATABASE:
        return {"success": False, "error": {"code": "INVALID_CATHODE", "message": f"Unknown cathode"}}
    
    if electrolyte not in ELECTROLYTE_DATABASE:
        return {"success": False, "error": {"code": "INVALID_ELECTROLYTE", "message": f"Unknown electrolyte"}}
    
    # Generate both structures
    cathode_struct = generate_cathode(cathode, supercell=[3, 3, 2])
    electrolyte_struct = generate_electrolyte(electrolyte, supercell=[2, 2, 2])
    
    if not cathode_struct["success"] or not electrolyte_struct["success"]:
        return {"success": False, "error": {"code": "GENERATION_FAILED", "message": "Failed to generate structures"}}
    
    # Combine atoms (simplified - stack vertically)
    cathode_atoms = cathode_struct["structure"]["atoms"]
    electrolyte_atoms = electrolyte_struct["structure"]["atoms"]
    
    combined_atoms = list(cathode_atoms)
    
    # Offset electrolyte atoms
    for atom in electrolyte_atoms:
        new_coords = [atom["coords"][0], atom["coords"][1], atom["coords"][2] + 1.1]
        combined_atoms.append({**atom, "coords": new_coords, "region": "electrolyte"})
    
    return {
        "success": True,
        "cathode": cathode,
        "electrolyte": electrolyte,
        "interface_orientation": interface_orientation,
        "n_cathode_atoms": len(cathode_atoms),
        "n_electrolyte_atoms": len(electrolyte_atoms),
        "n_total_atoms": len(combined_atoms),
        "structure": {"atoms": combined_atoms}
    }


def compare_electrolytes(
    electrolyte_list: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Compare properties of multiple electrolytes.
    
    Args:
        electrolyte_list: List of electrolytes to compare
    
    Returns:
        Comparison table
    """
    if electrolyte_list is None:
        electrolyte_list = list(ELECTROLYTE_DATABASE.keys())
    
    comparison = []
    for elec in electrolyte_list:
        if elec in ELECTROLYTE_DATABASE:
            info = ELECTROLYTE_DATABASE[elec]
            comparison.append({
                "electrolyte": elec,
                "formula": info["formula"],
                "structure": info["structure"],
                "conductivity_S_cm": info["conductivity"],
                "activation_eV": info["activation_eV"],
                "stability_low_V": info["stability_window"][0],
                "stability_high_V": info["stability_window"][1],
                "ion": info.get("ion", "Li")
            })
    
    # Sort by conductivity
    comparison.sort(key=lambda x: x["conductivity_S_cm"], reverse=True)
    
    return {
        "success": True,
        "n_electrolytes": len(comparison),
        "comparison": comparison,
        "highest_conductivity": comparison[0]["electrolyte"] if comparison else None
    }
