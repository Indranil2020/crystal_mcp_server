"""
electronic/semiconductors.py - Semiconductor Heterostructures

Comprehensive semiconductor generation per structure_catalogue.md Category 8:
- Type-I, Type-II, Type-III band alignment
- Quantum wells, superlattices
- 2DEG interfaces
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Semiconductor database
SEMICONDUCTOR_DATABASE = {
    # III-V semiconductors
    "GaAs": {"a": 5.65, "Eg_eV": 1.42, "type": "direct", "me": 0.067, "mh": 0.47, "eps": 13.1},
    "AlAs": {"a": 5.66, "Eg_eV": 2.16, "type": "indirect", "me": 0.15, "mh": 0.79, "eps": 10.1},
    "InAs": {"a": 6.06, "Eg_eV": 0.35, "type": "direct", "me": 0.026, "mh": 0.40, "eps": 15.1},
    "InP": {"a": 5.87, "Eg_eV": 1.35, "type": "direct", "me": 0.077, "mh": 0.60, "eps": 12.5},
    "GaP": {"a": 5.45, "Eg_eV": 2.26, "type": "indirect", "me": 0.13, "mh": 0.67, "eps": 11.1},
    "GaN": {"a": 3.19, "c": 5.19, "Eg_eV": 3.4, "type": "direct", "wurtzite": True},
    "AlN": {"a": 3.11, "c": 4.98, "Eg_eV": 6.2, "type": "direct", "wurtzite": True},
    "InN": {"a": 3.54, "c": 5.70, "Eg_eV": 0.7, "type": "direct", "wurtzite": True},
    "GaSb": {"a": 6.10, "Eg_eV": 0.73, "type": "direct", "me": 0.041, "eps": 15.7},
    "InSb": {"a": 6.48, "Eg_eV": 0.17, "type": "direct", "me": 0.014, "eps": 16.8},
    "AlSb": {"a": 6.14, "Eg_eV": 1.58, "type": "indirect"},
    
    # II-VI semiconductors
    "ZnS": {"a": 5.41, "Eg_eV": 3.68, "type": "direct"},
    "ZnSe": {"a": 5.67, "Eg_eV": 2.70, "type": "direct"},
    "ZnTe": {"a": 6.10, "Eg_eV": 2.26, "type": "direct"},
    "CdS": {"a": 5.82, "Eg_eV": 2.42, "type": "direct"},
    "CdSe": {"a": 6.05, "Eg_eV": 1.74, "type": "direct"},
    "CdTe": {"a": 6.48, "Eg_eV": 1.44, "type": "direct"},
    
    # Group IV
    "Si": {"a": 5.43, "Eg_eV": 1.12, "type": "indirect", "me": 1.08, "mh": 0.56, "eps": 11.7},
    "Ge": {"a": 5.66, "Eg_eV": 0.66, "type": "indirect", "me": 0.55, "mh": 0.37, "eps": 16.0},
    "SiGe": {"a": 5.54, "Eg_eV": 0.90, "type": "indirect", "alloy": True},
}


# Common heterostructure systems
HETEROSTRUCTURE_SYSTEMS = {
    "GaAs/AlGaAs": {"material1": "GaAs", "material2": "AlAs", "x_range": [0, 0.45],
                   "type_I": True, "2DEG": True, "HEMT": True},
    "InGaAs/InP": {"material1": "InAs", "material2": "InP", "x_range": [0.53, 0.53],
                  "lattice_matched": True, "type_I": True},
    "InAs/GaSb": {"material1": "InAs", "material2": "GaSb", 
                 "type_III": True, "broken_gap": True, "topological": True},
    "GaN/AlGaN": {"material1": "GaN", "material2": "AlN", "x_range": [0, 0.3],
                 "2DEG": True, "piezoelectric": True},
    "InGaN/GaN": {"material1": "InN", "material2": "GaN", "x_range": [0, 0.3],
                 "LED": True, "quantum_well": True},
    "SiGe/Si": {"material1": "Ge", "material2": "Si", "x_range": [0, 0.3],
               "type_II": True, "strained": True},
    "CdTe/HgTe": {"material1": "CdTe", "material2": "HgTe",
                 "topological": True, "inverted_gap": True},
}


# Band offsets (eV) for common interfaces
BAND_OFFSETS = {
    "GaAs/AlAs": {"VBO": 0.53, "CBO": 0.89, "alignment": "Type_I"},
    "InAs/GaSb": {"VBO": -0.56, "CBO": 0.90, "alignment": "Type_III"},
    "GaAs/InAs": {"VBO": 0.17, "CBO": 0.68, "alignment": "Type_I"},
    "ZnSe/ZnTe": {"VBO": 0.64, "CBO": -0.20, "alignment": "Type_II"},
    "CdTe/HgTe": {"VBO": 0.35, "CBO": -0.25, "alignment": "Type_III"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_quantum_well(
    well_material: str = "GaAs",
    barrier_material: str = "AlAs",
    well_thickness_nm: float = 10.0,
    barrier_thickness_nm: float = 20.0,
    n_wells: int = 1,
    alloy_fraction: float = 0.3
) -> Dict[str, Any]:
    """
    Generate quantum well heterostructure.
    
    Args:
        well_material: Well layer material
        barrier_material: Barrier layer material
        well_thickness_nm: Well thickness in nm
        barrier_thickness_nm: Barrier thickness in nm
        n_wells: Number of wells
        alloy_fraction: Alloy composition if barrier is alloy
    
    Returns:
        Quantum well structure
    """
    if well_material not in SEMICONDUCTOR_DATABASE:
        return {"success": False, "error": {"code": "INVALID_MATERIAL", "message": f"Unknown well material",
                "available": list(SEMICONDUCTOR_DATABASE.keys())}}
    
    if barrier_material not in SEMICONDUCTOR_DATABASE:
        return {"success": False, "error": {"code": "INVALID_MATERIAL", "message": f"Unknown barrier material"}}
    
    well_info = SEMICONDUCTOR_DATABASE[well_material]
    barrier_info = SEMICONDUCTOR_DATABASE[barrier_material]
    
    a_well = well_info["a"]
    a_barrier = barrier_info["a"]
    a_avg = (a_well + a_barrier) / 2
    
    well_thick = well_thickness_nm * 10  # Angstrom
    barrier_thick = barrier_thickness_nm * 10
    total_c = (well_thick + barrier_thick) * n_wells + barrier_thick
    
    lattice = Lattice.tetragonal(a_avg, total_c)
    
    species = []
    coords = []
    
    z = 0
    for period in range(n_wells + 1):
        # Barrier layer
        n_barrier_layers = int(barrier_thick / (a_barrier / 2))
        for layer in range(n_barrier_layers):
            z_frac = (z + layer * a_barrier / 2) / total_c
            if z_frac < 1:
                for i in range(2):
                    for j in range(2):
                        x = (i + 0.25 * (layer % 2)) / 2
                        y = (j + 0.25 * (layer % 2)) / 2
                        
                        # Alloy: mix elements
                        if alloy_fraction > 0 and period < n_wells:
                            elem = barrier_material[:2] if np.random.random() > alloy_fraction else well_material[:2]
                        else:
                            elem = barrier_material[:2]
                        species.append(elem)
                        coords.append([x % 1, y % 1, z_frac])
                        
                        # Anion
                        species.append(barrier_material[-2:] if len(barrier_material) > 2 else barrier_material[-1])
                        coords.append([(x + 0.25) % 1, (y + 0.25) % 1, z_frac])
        
        z += barrier_thick
        
        # Well layer
        if period < n_wells:
            n_well_layers = int(well_thick / (a_well / 2))
            for layer in range(n_well_layers):
                z_frac = (z + layer * a_well / 2) / total_c
                if z_frac < 1:
                    for i in range(2):
                        for j in range(2):
                            x = (i + 0.25 * (layer % 2)) / 2
                            y = (j + 0.25 * (layer % 2)) / 2
                            
                            species.append(well_material[:2])
                            coords.append([x % 1, y % 1, z_frac])
                            
                            species.append(well_material[-2:] if len(well_material) > 2 else well_material[-1])
                            coords.append([(x + 0.25) % 1, (y + 0.25) % 1, z_frac])
            
            z += well_thick
    
    structure = Structure(lattice, species, coords)
    
    # Calculate confinement energy (simple infinite well)
    me = well_info.get("me", 0.1)
    h = 4.136e-15  # eV*s
    m0 = 9.109e-31  # kg
    L = well_thickness_nm * 1e-9
    
    E1 = (h**2 / (8 * me * m0 * L**2)) * 1e20 * 6.242e18  # First level in eV
    
    return {
        "success": True,
        "well_material": well_material,
        "barrier_material": barrier_material,
        "well_thickness_nm": well_thickness_nm,
        "barrier_thickness_nm": barrier_thickness_nm,
        "n_wells": n_wells,
        "alloy_fraction": alloy_fraction,
        "confinement_energy_meV": round(E1 * 1000, 1),
        "well_bandgap_eV": well_info["Eg_eV"],
        "barrier_bandgap_eV": barrier_info["Eg_eV"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_superlattice(
    material1: str = "GaAs",
    material2: str = "AlAs",
    period1_nm: float = 5.0,
    period2_nm: float = 5.0,
    n_periods: int = 5
) -> Dict[str, Any]:
    """
    Generate superlattice structure.
    
    Args:
        material1: First layer material
        material2: Second layer material
        period1_nm: First layer thickness
        period2_nm: Second layer thickness
        n_periods: Number of periods
    
    Returns:
        Superlattice structure
    """
    result = generate_quantum_well(
        material1, material2,
        well_thickness_nm=period1_nm,
        barrier_thickness_nm=period2_nm,
        n_wells=n_periods
    )
    
    if result["success"]:
        result["structure_type"] = "superlattice"
        result["period_nm"] = period1_nm + period2_nm
        result["total_periods"] = n_periods
        
        # Miniband formation
        result["miniband_formation"] = n_periods >= 3
    
    return result


def generate_2deg_interface(
    system: str = "GaAs/AlGaAs",
    spacer_nm: float = 5.0,
    doping_nm: float = 10.0,
    doping_density_cm2: float = 1e12
) -> Dict[str, Any]:
    """
    Generate 2DEG heterostructure interface.
    
    Args:
        system: Heterostructure system
        spacer_nm: Undoped spacer thickness
        doping_nm: Doped region thickness
        doping_density_cm2: 2D doping density
    
    Returns:
        2DEG interface structure
    """
    if system not in HETEROSTRUCTURE_SYSTEMS:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", "message": f"Unknown system",
                      "available": list(HETEROSTRUCTURE_SYSTEMS.keys())}
        }
    
    sys_info = HETEROSTRUCTURE_SYSTEMS[system]
    mat1 = sys_info["material1"]
    mat2 = sys_info["material2"]
    
    # Generate base structure
    result = generate_quantum_well(mat1, mat2, 50.0, 30.0, 1)
    
    if result["success"]:
        result["system"] = system
        result["has_2DEG"] = True
        result["spacer_nm"] = spacer_nm
        result["doping_nm"] = doping_nm
        result["doping_density_cm2"] = doping_density_cm2
        result["is_HEMT"] = sys_info.get("HEMT", False)
        result["is_piezoelectric"] = sys_info.get("piezoelectric", False)
        result["is_topological"] = sys_info.get("topological", False)
    
    return result
