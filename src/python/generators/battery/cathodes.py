"""
battery/cathodes.py - Battery Cathode Materials

Comprehensive cathode material generation:
- Layered oxides (NCM, NCA, LCO)
- Olivines (LFP, LMP)
- Spinels (LMO)
- Disordered rock-salt (DRX)
- Lithiation/delithiation states
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Comprehensive cathode database
CATHODE_DATABASE = {
    # Layered oxides
    "LCO": {
        "formula": "LiCoO2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.816, "c": 14.05, "voltage": 3.9, "capacity_mAh_g": 274,
        "elements": {"Li": 1, "Co": 1, "O": 2}
    },
    "NCM111": {
        "formula": "LiNi1/3Co1/3Mn1/3O2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.867, "c": 14.23, "voltage": 3.7, "capacity_mAh_g": 280,
        "elements": {"Li": 1, "Ni": 0.33, "Co": 0.33, "Mn": 0.33, "O": 2}
    },
    "NCM523": {
        "formula": "LiNi0.5Co0.2Mn0.3O2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.871, "c": 14.26, "voltage": 3.7, "capacity_mAh_g": 275,
        "elements": {"Li": 1, "Ni": 0.5, "Co": 0.2, "Mn": 0.3, "O": 2}
    },
    "NCM622": {
        "formula": "LiNi0.6Co0.2Mn0.2O2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.875, "c": 14.28, "voltage": 3.65, "capacity_mAh_g": 285,
        "elements": {"Li": 1, "Ni": 0.6, "Co": 0.2, "Mn": 0.2, "O": 2}
    },
    "NCM811": {
        "formula": "LiNi0.8Co0.1Mn0.1O2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.881, "c": 14.31, "voltage": 3.6, "capacity_mAh_g": 290,
        "elements": {"Li": 1, "Ni": 0.8, "Co": 0.1, "Mn": 0.1, "O": 2}
    },
    "NCA": {
        "formula": "LiNi0.8Co0.15Al0.05O2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.865, "c": 14.18, "voltage": 3.7, "capacity_mAh_g": 200,
        "elements": {"Li": 1, "Ni": 0.8, "Co": 0.15, "Al": 0.05, "O": 2}
    },
    "LNO": {
        "formula": "LiNiO2", "structure": "layered", "spacegroup": "R-3m",
        "a": 2.879, "c": 14.19, "voltage": 3.8, "capacity_mAh_g": 275,
        "elements": {"Li": 1, "Ni": 1, "O": 2}
    },
    
    # Olivines
    "LFP": {
        "formula": "LiFePO4", "structure": "olivine", "spacegroup": "Pnma",
        "a": 10.334, "b": 6.008, "c": 4.694, "voltage": 3.4, "capacity_mAh_g": 170,
        "elements": {"Li": 1, "Fe": 1, "P": 1, "O": 4}
    },
    "LMP": {
        "formula": "LiMnPO4", "structure": "olivine", "spacegroup": "Pnma",
        "a": 10.45, "b": 6.10, "c": 4.75, "voltage": 4.1, "capacity_mAh_g": 171,
        "elements": {"Li": 1, "Mn": 1, "P": 1, "O": 4}
    },
    "LFMP": {
        "formula": "LiFe0.5Mn0.5PO4", "structure": "olivine", "spacegroup": "Pnma",
        "a": 10.39, "b": 6.05, "c": 4.72, "voltage": 3.75, "capacity_mAh_g": 170,
        "elements": {"Li": 1, "Fe": 0.5, "Mn": 0.5, "P": 1, "O": 4}
    },
    "LCP": {
        "formula": "LiCoPO4", "structure": "olivine", "spacegroup": "Pnma",
        "a": 10.20, "b": 5.92, "c": 4.69, "voltage": 4.8, "capacity_mAh_g": 167,
        "elements": {"Li": 1, "Co": 1, "P": 1, "O": 4}
    },
    
    # Spinels
    "LMO": {
        "formula": "LiMn2O4", "structure": "spinel", "spacegroup": "Fd-3m",
        "a": 8.24, "voltage": 4.0, "capacity_mAh_g": 148,
        "elements": {"Li": 1, "Mn": 2, "O": 4}
    },
    "LNMO": {
        "formula": "LiNi0.5Mn1.5O4", "structure": "spinel", "spacegroup": "P4332",
        "a": 8.17, "voltage": 4.7, "capacity_mAh_g": 147,
        "elements": {"Li": 1, "Ni": 0.5, "Mn": 1.5, "O": 4}
    },
    
    # Li-rich
    "Li-rich_NCM": {
        "formula": "Li1.2Ni0.2Mn0.6O2", "structure": "layered", "spacegroup": "C2/m",
        "a": 4.93, "b": 8.54, "c": 5.02, "voltage": 3.5, "capacity_mAh_g": 300,
        "elements": {"Li": 1.2, "Ni": 0.2, "Mn": 0.6, "O": 2}
    },
    
    # Disordered rock-salt
    "DRX_LiMnO2F": {
        "formula": "Li1.2Mn0.6Ti0.2O1.8F0.2", "structure": "DRX", "spacegroup": "Fm-3m",
        "a": 4.20, "voltage": 3.3, "capacity_mAh_g": 250,
        "elements": {"Li": 1.2, "Mn": 0.6, "Ti": 0.2, "O": 1.8, "F": 0.2}
    },
    
    # Sodium cathodes
    "NaCoO2": {
        "formula": "NaCoO2", "structure": "layered", "spacegroup": "P63/mmc",
        "a": 2.889, "c": 10.81, "voltage": 3.0, "capacity_mAh_g": 235,
        "elements": {"Na": 1, "Co": 1, "O": 2}, "ion": "Na"
    },
    "Na-NCM": {
        "formula": "NaNi1/3Co1/3Mn1/3O2", "structure": "layered", "spacegroup": "P63/mmc",
        "a": 2.96, "c": 11.2, "voltage": 3.2, "capacity_mAh_g": 200,
        "elements": {"Na": 1, "Ni": 0.33, "Co": 0.33, "Mn": 0.33, "O": 2}, "ion": "Na"
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


def generate_cathode(
    cathode_type: str = "NCM111",
    lithiation: float = 1.0,
    supercell: List[int] = [2, 2, 1]
) -> Dict[str, Any]:
    """
    Generate cathode structure.
    
    Args:
        cathode_type: Cathode type from database
        lithiation: Li content (0.0 to 1.0)
        supercell: Supercell dimensions
    
    Returns:
        Cathode structure
    """
    if cathode_type not in CATHODE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_CATHODE", "message": f"Unknown cathode '{cathode_type}'",
                      "available": list(CATHODE_DATABASE.keys())}
        }
    
    info = CATHODE_DATABASE[cathode_type]
    structure_type = info["structure"]
    
    if structure_type == "layered":
        a = info["a"]
        c = info["c"]
        lattice = Lattice.hexagonal(a, c)
        
        # R-3m structure
        elements = info["elements"]
        species = []
        coords = []
        
        # Li layer
        if lithiation > 0:
            species.append("Li")
            coords.append([0, 0, 0.5])
        
        # TM layer
        for elem, frac in elements.items():
            if elem not in ["Li", "O", "Na"]:
                species.append(elem)
                coords.append([0, 0, 0])
                break
        
        # O layers
        species.extend(["O", "O"])
        coords.extend([[0, 0, 0.24], [0, 0, 0.76]])
        
    elif structure_type == "olivine":
        a, b, c = info["a"], info["b"], info["c"]
        lattice = Lattice.orthorhombic(a, b, c)
        
        # Pnma structure (simplified)
        species = []
        coords = []
        
        if lithiation > 0:
            species.append("Li")
            coords.append([0, 0, 0])
        
        species.append("Fe" if "Fe" in info["elements"] else "Mn")
        coords.append([0.28, 0.25, 0.97])
        
        species.append("P")
        coords.append([0.09, 0.25, 0.42])
        
        species.extend(["O"] * 4)
        coords.extend([
            [0.1, 0.25, 0.74],
            [0.45, 0.25, 0.21],
            [0.16, 0.04, 0.29],
            [0.16, 0.46, 0.29]
        ])
        
    elif structure_type == "spinel":
        a = info["a"]
        lattice = Lattice.cubic(a)
        
        # Fd-3m spinel
        species = ["Li", "Mn", "Mn", "O", "O", "O", "O"]
        coords = [
            [0.125, 0.125, 0.125],
            [0.5, 0.5, 0.5],
            [0.5, 0, 0],
            [0.26, 0.26, 0.26],
            [0.74, 0.74, 0.74],
            [0.26, 0.26, 0.74],
            [0.74, 0.26, 0.26]
        ]
    
    else:  # DRX or generic
        a = info["a"]
        lattice = Lattice.cubic(a)
        species = ["Li", "Mn", "O", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
    
    structure = Structure(lattice, species, coords)
    
    nx, ny, nz = supercell
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    # Apply lithiation (remove Li atoms)
    if lithiation < 1.0:
        n_li = sum(1 for s in structure if str(s.specie) == "Li")
        n_remove = int(n_li * (1 - lithiation))
        
        if n_remove > 0:
            li_indices = [i for i, s in enumerate(structure) if str(s.specie) == "Li"]
            np.random.seed(42)
            remove_indices = np.random.choice(li_indices, n_remove, replace=False)
            structure.remove_sites(sorted(remove_indices, reverse=True))
    
    # Calculate practical capacity
    practical_capacity = info["capacity_mAh_g"] * lithiation
    
    return {
        "success": True,
        "cathode_type": cathode_type,
        "formula": info["formula"],
        "structure_type": structure_type,
        "spacegroup": info["spacegroup"],
        "voltage_V": info["voltage"],
        "theoretical_capacity_mAh_g": info["capacity_mAh_g"],
        "lithiation": lithiation,
        "practical_capacity_mAh_g": round(practical_capacity, 1),
        "ion": info.get("ion", "Li"),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_lithiated_series(
    cathode_type: str = "LCO",
    n_points: int = 5
) -> Dict[str, Any]:
    """
    Generate series of lithiated structures.
    
    Args:
        cathode_type: Cathode type
        n_points: Number of lithiation states
    
    Returns:
        Series of structures from x=0 to x=1
    """
    lithiations = np.linspace(0, 1, n_points)
    structures = []
    
    for x in lithiations:
        result = generate_cathode(cathode_type, lithiation=x)
        if result["success"]:
            structures.append({
                "lithiation": round(x, 2),
                "voltage_V": result["voltage_V"],
                "capacity_mAh_g": result["practical_capacity_mAh_g"],
                "n_atoms": result["n_atoms"]
            })
    
    return {
        "success": True,
        "cathode_type": cathode_type,
        "n_structures": len(structures),
        "series": structures
    }


def generate_surface_coated_cathode(
    cathode_type: str = "NCM811",
    coating: str = "Al2O3",
    coating_thickness_nm: float = 2.0,
    supercell: List[int] = [3, 3, 2]
) -> Dict[str, Any]:
    """
    Generate surface-coated cathode.
    
    Args:
        cathode_type: Base cathode
        coating: Coating material (Al2O3, ZrO2, TiO2, LiNbO3)
        coating_thickness_nm: Coating thickness
        supercell: Supercell size
    
    Returns:
        Coated cathode structure
    """
    coatings = {
        "Al2O3": {"elements": ["Al", "O"], "a": 4.76, "protection": "HF_scavenger"},
        "ZrO2": {"elements": ["Zr", "O"], "a": 5.15, "protection": "structural"},
        "TiO2": {"elements": ["Ti", "O"], "a": 4.59, "protection": "electronic"},
        "LiNbO3": {"elements": ["Li", "Nb", "O"], "a": 5.15, "protection": "ionic_conductor"},
        "Li3PO4": {"elements": ["Li", "P", "O"], "a": 6.11, "protection": "ionic_conductor"},
    }
    
    if coating not in coatings:
        return {"success": False, "error": {"code": "INVALID_COATING", "message": f"Unknown coating"}}
    
    # Generate base cathode
    base = generate_cathode(cathode_type, supercell=supercell)
    
    if not base["success"]:
        return base
    
    atoms = base["structure"]["atoms"]
    coating_info = coatings[coating]
    
    # Add coating layer on top surface
    base_lattice = base["structure"]["lattice"]
    c = base_lattice.get("c", 10)
    
    thickness = coating_thickness_nm * 10  # Convert to Angstrom
    n_coating_atoms = int(thickness * 5)  # Approximate
    
    np.random.seed(42)
    for _ in range(n_coating_atoms):
        elem = np.random.choice(coating_info["elements"])
        x = np.random.uniform(0, 1)
        y = np.random.uniform(0, 1)
        z = 1.0 + np.random.uniform(0.05, 0.05 + thickness / c)
        
        atoms.append({
            "element": elem,
            "coords": [x, y, z],
            "coating": True
        })
    
    return {
        "success": True,
        "cathode_type": cathode_type,
        "coating": coating,
        "coating_thickness_nm": coating_thickness_nm,
        "protection_mechanism": coating_info["protection"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": base_lattice}
    }
