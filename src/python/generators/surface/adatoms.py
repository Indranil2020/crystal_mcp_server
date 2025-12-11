"""
surface/adatoms.py - Adatom and Advacancy Superlattices

Generates ordered adatom/advacancy patterns on surfaces:
- O on Ru(0001) - p(2×2), (2×1)
- N on Fe(100) - c(2×2)
- S on Ni(111) - (√3×√3)R30°
- Various coverage patterns
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Adatom pattern database
ADATOM_PATTERNS = {
    "p(2x2)": {
        "supercell": [2, 2, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 0.25,
        "description": "1/4 ML coverage, square pattern"
    },
    "p(1x1)": {
        "supercell": [1, 1, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 1.0,
        "description": "Full monolayer"
    },
    "c(2x2)": {
        "supercell": [2, 2, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 0.25,
        "description": "Centered 2×2, 1/4 ML"
    },
    "(2x1)": {
        "supercell": [2, 1, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 0.5,
        "description": "1/2 ML row pattern"
    },
    "(sqrt3xsqrt3)R30": {
        "supercell": [2, 2, 1],  # Approximation
        "adatom_positions": [[0.33, 0.33], [0.67, 0.67]],
        "coverage": 0.33,
        "description": "(√3×√3)R30° reconstruction"
    },
    "p(3x1)": {
        "supercell": [3, 1, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 0.33,
        "description": "1/3 ML stripe pattern"
    },
    "p(4x4)": {
        "supercell": [4, 4, 1],
        "adatom_positions": [[0.5, 0.5]],
        "coverage": 0.0625,
        "description": "Isolated adatom, 1/16 ML"
    },
    "honeycomb": {
        "supercell": [2, 2, 1],
        "adatom_positions": [[0.33, 0.67], [0.67, 0.33]],
        "coverage": 0.5,
        "description": "Honeycomb adatom pattern"
    },
    "kagome": {
        "supercell": [2, 2, 1],
        "adatom_positions": [[0.25, 0.5], [0.5, 0.25], [0.75, 0.75]],
        "coverage": 0.75,
        "description": "Kagome adatom lattice"
    },
}

# Common adatom systems
ADATOM_SYSTEMS = {
    "O-Ru(0001)": {"substrate": "Ru", "adatom": "O", "facet": "0001", "height": 1.2, "a": 2.70},
    "O-Pt(111)": {"substrate": "Pt", "adatom": "O", "facet": "111", "height": 1.3, "a": 2.77},
    "N-Fe(100)": {"substrate": "Fe", "adatom": "N", "facet": "100", "height": 0.3, "a": 2.87},
    "S-Ni(111)": {"substrate": "Ni", "adatom": "S", "facet": "111", "height": 1.5, "a": 2.49},
    "S-Au(111)": {"substrate": "Au", "adatom": "S", "facet": "111", "height": 1.8, "a": 2.88},
    "Cl-Ag(111)": {"substrate": "Ag", "adatom": "Cl", "facet": "111", "height": 2.0, "a": 2.89},
    "C-Ni(111)": {"substrate": "Ni", "adatom": "C", "facet": "111", "height": 1.0, "a": 2.49},
    "H-Pd(111)": {"substrate": "Pd", "adatom": "H", "facet": "111", "height": 0.9, "a": 2.75},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords), "cartesian": list(s.coords)}
                  for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_adatom_superlattice(
    system: str = "O-Ru(0001)",
    pattern: str = "p(2x2)",
    n_layers: int = 4,
    vacuum: float = 15.0,
    site_type: str = "fcc"
) -> Dict[str, Any]:
    """
    Generate adatom superlattice on metal surface.
    
    Args:
        system: Pre-defined system (O-Ru(0001), N-Fe(100), etc.)
        pattern: Adatom pattern (p(2x2), c(2x2), (√3×√3)R30°, etc.)
        n_layers: Number of substrate layers
        vacuum: Vacuum thickness
        site_type: Adsorption site (fcc, hcp, top, bridge)
    
    Returns:
        Adatom superlattice structure
    """
    if system not in ADATOM_SYSTEMS:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", "message": f"Unknown system '{system}'",
                      "available": list(ADATOM_SYSTEMS.keys())}
        }
    
    if pattern not in ADATOM_PATTERNS:
        return {
            "success": False,
            "error": {"code": "INVALID_PATTERN", "message": f"Unknown pattern '{pattern}'",
                      "available": list(ADATOM_PATTERNS.keys())}
        }
    
    sys_info = ADATOM_SYSTEMS[system]
    pat_info = ADATOM_PATTERNS[pattern]
    
    a = sys_info["a"]
    substrate = sys_info["substrate"]
    adatom = sys_info["adatom"]
    ad_height = sys_info["height"]
    
    supercell = pat_info["supercell"]
    adatom_frac = pat_info["adatom_positions"]
    
    # Build substrate slab
    a_super = a * supercell[0]
    b_super = a * supercell[1]
    z_sep = a / np.sqrt(3)  # FCC/HCP interlayer
    
    c = n_layers * z_sep + ad_height + vacuum
    
    if "111" in sys_info["facet"] or "0001" in sys_info["facet"]:
        lattice = Lattice.hexagonal(a_super, c)
    else:
        lattice = Lattice.orthorhombic(a_super, b_super, c)
    
    species = []
    coords = []
    
    # Substrate atoms
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c
        
        # Shift for ABC stacking
        shift_x = (layer % 3) / 3 / supercell[0]
        shift_y = (layer % 3) / 3 / supercell[1]
        
        for i in range(supercell[0]):
            for j in range(supercell[1]):
                x = (i + 0.5 + shift_x * supercell[0]) / supercell[0]
                y = (j + 0.5 + shift_y * supercell[1]) / supercell[1]
                
                species.append(substrate)
                coords.append([x % 1, y % 1, z_layer])
    
    # Adatom layer
    z_adatom = (n_layers * z_sep + ad_height + vacuum / 2) / c
    
    # Site position adjustments
    site_offsets = {
        "fcc": [0, 0],
        "hcp": [0.17, 0.17],
        "top": [0.5 / supercell[0], 0.5 / supercell[1]],
        "bridge": [0.25 / supercell[0], 0.25 / supercell[1]],
    }
    offset = site_offsets.get(site_type, [0, 0])
    
    for pos in adatom_frac:
        x = (pos[0] + offset[0]) % 1
        y = (pos[1] + offset[1]) % 1
        species.append(adatom)
        coords.append([x, y, z_adatom])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "system": system,
        "pattern": pattern,
        "coverage_ML": pat_info["coverage"],
        "site_type": site_type,
        "n_substrate_atoms": n_layers * supercell[0] * supercell[1],
        "n_adatoms": len(adatom_frac),
        "adatom_height_angstrom": ad_height,
        "structure": structure_to_dict(structure)
    }


def generate_advacancy_superlattice(
    element: str = "Cu",
    facet: str = "111",
    pattern: str = "p(2x2)",
    n_layers: int = 4,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate ordered vacancy pattern on surface.
    
    Args:
        element: Surface element
        facet: Surface facet
        pattern: Vacancy pattern
        n_layers: Number of layers
        vacuum: Vacuum thickness
    
    Returns:
        Surface with ordered vacancies
    """
    if pattern not in ADATOM_PATTERNS:
        return {"success": False, "error": {"code": "INVALID_PATTERN", "message": f"Unknown: {pattern}"}}
    
    pat_info = ADATOM_PATTERNS[pattern]
    supercell = pat_info["supercell"]
    vacancy_frac = pat_info["adatom_positions"]
    
    a_elem = {"Cu": 2.56, "Ag": 2.89, "Au": 2.88, "Pt": 2.77, "Pd": 2.75, "Ni": 2.49}.get(element, 2.5)
    
    a_super = a_elem * supercell[0]
    z_sep = a_elem / np.sqrt(3)
    c = n_layers * z_sep + vacuum
    
    lattice = Lattice.hexagonal(a_super, c)
    
    species = []
    coords = []
    
    # Build slab and remove vacancies from top layer
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c
        is_top_layer = (layer == n_layers - 1)
        
        shift_x = (layer % 3) / 3 / supercell[0]
        shift_y = (layer % 3) / 3 / supercell[1]
        
        for i in range(supercell[0]):
            for j in range(supercell[1]):
                x = (i + 0.5 + shift_x * supercell[0]) / supercell[0]
                y = (j + 0.5 + shift_y * supercell[1]) / supercell[1]
                
                # Check if this is a vacancy position in top layer
                is_vacancy = False
                if is_top_layer:
                    for vac_pos in vacancy_frac:
                        if abs(x - vac_pos[0]) < 0.1 and abs(y - vac_pos[1]) < 0.1:
                            is_vacancy = True
                            break
                
                if not is_vacancy:
                    species.append(element)
                    coords.append([x % 1, y % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "facet": facet,
        "pattern": pattern,
        "vacancy_concentration": pat_info["coverage"],
        "n_vacancies": len(vacancy_frac),
        "structure": structure_to_dict(structure)
    }
