"""
adsorption/ionic_liquids.py - Ionic Liquid Interfaces

Comprehensive ionic liquid interface generation:
- Room temperature ionic liquids (RTILs)
- Electrode-ionic liquid interfaces
- Confined ionic liquids
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Ionic liquid database
IONIC_LIQUID_DATABASE = {
    # Imidazolium-based
    "EMIM_BF4": {
        "cation": "EMIM", "anion": "BF4",
        "cation_formula": "C6H11N2+", "anion_formula": "BF4-",
        "cation_size_A": 5.5, "anion_size_A": 2.8,
        "density_gcm3": 1.29, "viscosity_mPas": 43
    },
    "BMIM_PF6": {
        "cation": "BMIM", "anion": "PF6",
        "cation_formula": "C8H15N2+", "anion_formula": "PF6-",
        "cation_size_A": 6.5, "anion_size_A": 3.2,
        "density_gcm3": 1.38, "viscosity_mPas": 312
    },
    "BMIM_BF4": {
        "cation": "BMIM", "anion": "BF4",
        "cation_formula": "C8H15N2+", "anion_formula": "BF4-",
        "cation_size_A": 6.5, "anion_size_A": 2.8,
        "density_gcm3": 1.21, "viscosity_mPas": 154
    },
    "BMIM_TFSI": {
        "cation": "BMIM", "anion": "TFSI",
        "cation_formula": "C8H15N2+", "anion_formula": "C2F6NO4S2-",
        "cation_size_A": 6.5, "anion_size_A": 5.5,
        "density_gcm3": 1.44, "viscosity_mPas": 52,
        "hydrophobic": True
    },
    "EMIM_TFSI": {
        "cation": "EMIM", "anion": "TFSI",
        "cation_formula": "C6H11N2+", "anion_formula": "C2F6NO4S2-",
        "cation_size_A": 5.5, "anion_size_A": 5.5,
        "density_gcm3": 1.52, "viscosity_mPas": 34
    },
    
    # Pyrrolidinium-based
    "PYR14_TFSI": {
        "cation": "PYR14", "anion": "TFSI",
        "cation_formula": "C9H20N+", "anion_formula": "C2F6NO4S2-",
        "cation_size_A": 6.0, "anion_size_A": 5.5,
        "density_gcm3": 1.41, "electrochemical_window_V": 5.5
    },
    "PYR13_FSI": {
        "cation": "PYR13", "anion": "FSI",
        "cation_formula": "C8H18N+", "anion_formula": "F2NO4S2-",
        "cation_size_A": 5.8, "anion_size_A": 4.5,
        "density_gcm3": 1.34
    },
    
    # Ammonium-based
    "N1114_TFSI": {
        "cation": "N1114", "anion": "TFSI",
        "cation_formula": "C7H18N+", "anion_formula": "C2F6NO4S2-",
        "cation_size_A": 5.5, "anion_size_A": 5.5,
        "electrochemical_window_V": 6.0
    },
}


# Electrode surface database
ELECTRODE_SURFACES = {
    "Au_111": {"element": "Au", "miller": [1, 1, 1], "a": 2.88, "workfunction_eV": 5.3},
    "Pt_111": {"element": "Pt", "miller": [1, 1, 1], "a": 2.77, "workfunction_eV": 5.65},
    "graphite_0001": {"element": "C", "miller": [0, 0, 0, 1], "a": 2.46, "workfunction_eV": 4.6},
    "graphene": {"element": "C", "miller": [0, 0, 0, 1], "a": 2.46, "workfunction_eV": 4.5},
    "LiFePO4_010": {"formula": "LiFePO4", "miller": [0, 1, 0], "battery": True},
    "Li_metal": {"element": "Li", "a": 3.50, "alkali": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_ionic_liquid_bulk(
    ionic_liquid: str = "BMIM_BF4",
    box_size_A: float = 30.0,
    n_ion_pairs: int = 20
) -> Dict[str, Any]:
    """
    Generate bulk ionic liquid structure.
    
    Args:
        ionic_liquid: Ionic liquid from database
        box_size_A: Box size in Angstrom
        n_ion_pairs: Number of ion pairs
    
    Returns:
        Ionic liquid structure
    """
    if ionic_liquid not in IONIC_LIQUID_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_IL", "message": f"Unknown ionic liquid",
                      "available": list(IONIC_LIQUID_DATABASE.keys())}
        }
    
    info = IONIC_LIQUID_DATABASE[ionic_liquid]
    
    lattice = Lattice.cubic(box_size_A)
    
    species = []
    coords = []
    
    cation_size = info["cation_size_A"]
    anion_size = info["anion_size_A"]
    
    np.random.seed(42)
    
    # Place ion pairs randomly
    for i in range(n_ion_pairs):
        # Cation (represented as single atom)
        while True:
            x = np.random.random()
            y = np.random.random()
            z = np.random.random()
            
            # Check minimum distance from existing atoms
            too_close = False
            for existing in coords:
                dist = np.sqrt(sum((x - existing[0])**2 + (y - existing[1])**2 + (z - existing[2])**2))
                if dist * box_size_A < 3.0:
                    too_close = True
                    break
            
            if not too_close:
                break
        
        species.append("N")  # Represent cation center
        coords.append([x, y, z])
        
        # Anion nearby
        offset = 4.0 / box_size_A  # ~4 Angstrom away
        anion_x = (x + offset * (np.random.random() - 0.5) * 2) % 1
        anion_y = (y + offset * (np.random.random() - 0.5) * 2) % 1
        anion_z = (z + offset * (np.random.random() - 0.5) * 2) % 1
        
        if "BF4" in ionic_liquid:
            species.append("B")
        elif "PF6" in ionic_liquid:
            species.append("P")
        elif "TFSI" in ionic_liquid or "FSI" in ionic_liquid:
            species.append("S")
        else:
            species.append("X")
        coords.append([anion_x, anion_y, anion_z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "ionic_liquid": ionic_liquid,
        "cation": info["cation"],
        "anion": info["anion"],
        "n_ion_pairs": n_ion_pairs,
        "box_size_A": box_size_A,
        "density_gcm3": info["density_gcm3"],
        "viscosity_mPas": info.get("viscosity_mPas", 0),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_electrode_il_interface(
    electrode: str = "Au_111",
    ionic_liquid: str = "BMIM_BF4",
    electrode_layers: int = 4,
    il_thickness_A: float = 30.0,
    n_ion_pairs: int = 15
) -> Dict[str, Any]:
    """
    Generate electrode-ionic liquid interface.
    
    Args:
        electrode: Electrode surface
        ionic_liquid: Ionic liquid
        electrode_layers: Number of electrode layers
        il_thickness_A: IL region thickness
        n_ion_pairs: Number of ion pairs
    
    Returns:
        Interface structure
    """
    if electrode not in ELECTRODE_SURFACES:
        return {"success": False, "error": {"code": "INVALID_ELECTRODE", "message": f"Unknown electrode"}}
    
    if ionic_liquid not in IONIC_LIQUID_DATABASE:
        return {"success": False, "error": {"code": "INVALID_IL", "message": f"Unknown IL"}}
    
    elec_info = ELECTRODE_SURFACES[electrode]
    il_info = IONIC_LIQUID_DATABASE[ionic_liquid]
    
    a = elec_info.get("a", 4.0)
    surface_size = a * 5
    
    electrode_height = electrode_layers * 2.5
    total_c = electrode_height + il_thickness_A + 10
    
    lattice = Lattice.orthorhombic(surface_size, surface_size * np.sqrt(3), total_c)
    
    species = []
    coords = []
    
    # Build electrode slab
    elem = elec_info.get("element", electrode.split("_")[0])
    for layer in range(electrode_layers):
        z = (layer * 2.5) / total_c
        for i in range(5):
            for j in range(5):
                x = (i + 0.5 * (layer % 2)) / 5
                y = (j + (layer % 3) * 0.33) / 5
                species.append(elem)
                coords.append([x % 1, y % 1, z])
    
    # Add ionic liquid
    il_start_z = electrode_height / total_c
    
    np.random.seed(42)
    
    for i in range(n_ion_pairs):
        # Cation
        z = il_start_z + (np.random.random() * il_thickness_A) / total_c
        x = np.random.random()
        y = np.random.random()
        species.append("N")
        coords.append([x, y, z])
        
        # Anion
        z_anion = il_start_z + (np.random.random() * il_thickness_A) / total_c
        species.append("B" if "BF4" in ionic_liquid else "S")
        coords.append([np.random.random(), np.random.random(), z_anion])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "electrode": electrode,
        "ionic_liquid": ionic_liquid,
        "electrode_layers": electrode_layers,
        "il_thickness_A": il_thickness_A,
        "n_ion_pairs": n_ion_pairs,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_confined_il(
    ionic_liquid: str = "EMIM_BF4",
    confinement: str = "nanotube",
    pore_diameter_A: float = 20.0
) -> Dict[str, Any]:
    """
    Generate confined ionic liquid.
    
    Args:
        ionic_liquid: Ionic liquid
        confinement: 'nanotube', 'graphene_slit', 'micropore'
        pore_diameter_A: Pore diameter
    
    Returns:
        Confined IL structure
    """
    if ionic_liquid not in IONIC_LIQUID_DATABASE:
        return {"success": False, "error": {"code": "INVALID_IL", "message": f"Unknown IL"}}
    
    il_info = IONIC_LIQUID_DATABASE[ionic_liquid]
    
    length = 50.0
    
    if confinement == "nanotube":
        lattice = Lattice.orthorhombic(pore_diameter_A + 10, pore_diameter_A + 10, length)
    else:
        lattice = Lattice.orthorhombic(30, 30, pore_diameter_A + 20)
    
    species = []
    coords = []
    
    # Add confining wall atoms
    if confinement == "nanotube":
        r = pore_diameter_A / 2 + 3.4  # Wall radius
        for z_idx in range(int(length / 2)):
            for theta_idx in range(12):
                theta = theta_idx * 2 * np.pi / 12
                x = 0.5 + (r * np.cos(theta)) / (pore_diameter_A + 10)
                y = 0.5 + (r * np.sin(theta)) / (pore_diameter_A + 10)
                z = z_idx / (length / 2)
                species.append("C")
                coords.append([x % 1, y % 1, z])
    
    # Add IL inside
    n_pairs = int((pore_diameter_A / 5) ** 2)
    np.random.seed(42)
    
    for i in range(n_pairs):
        if confinement == "nanotube":
            r_ion = np.random.random() * (pore_diameter_A / 2 - 3) / (pore_diameter_A + 10)
            theta = np.random.random() * 2 * np.pi
            x = 0.5 + r_ion * np.cos(theta)
            y = 0.5 + r_ion * np.sin(theta)
            z = 0.1 + np.random.random() * 0.8
        else:
            x = 0.2 + np.random.random() * 0.6
            y = 0.2 + np.random.random() * 0.6
            z = 0.3 + np.random.random() * 0.4
        
        species.append("N")
        coords.append([x, y, z])
        species.append("B" if "BF4" in ionic_liquid else "S")
        coords.append([x + 0.05, y + 0.05, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "ionic_liquid": ionic_liquid,
        "confinement": confinement,
        "pore_diameter_A": pore_diameter_A,
        "n_ion_pairs": n_pairs,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
