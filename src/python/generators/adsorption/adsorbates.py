"""
adsorption/adsorbates.py - Adsorbate Molecules on Surfaces

Comprehensive adsorbate generation per structure_catalogue.md Category 6:
(i) Small molecule adsorbates (H, O, N, C, CO, NO, NO2, CO2, H2O, NH3, CH4, C2H4, etc.)
(ii) Atomic coverage patterns (1/4 ML, 1/2 ML, full ML)
(iii) Co-adsorption configurations
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice, Molecule


# Adsorbate molecule database
ADSORBATE_DATABASE = {
    # Atomic adsorbates
    "H": {"atoms": ["H"], "coords": [[0, 0, 0]], "adsorption_height_A": 1.5, "binding": "covalent"},
    "O": {"atoms": ["O"], "coords": [[0, 0, 0]], "adsorption_height_A": 1.3, "binding": "covalent"},
    "N": {"atoms": ["N"], "coords": [[0, 0, 0]], "adsorption_height_A": 1.4, "binding": "covalent"},
    "C": {"atoms": ["C"], "coords": [[0, 0, 0]], "adsorption_height_A": 1.5, "binding": "covalent"},
    "S": {"atoms": ["S"], "coords": [[0, 0, 0]], "adsorption_height_A": 1.8, "binding": "covalent"},
    "Cl": {"atoms": ["Cl"], "coords": [[0, 0, 0]], "adsorption_height_A": 2.0, "binding": "ionic"},
    
    # Diatomic molecules
    "CO": {"atoms": ["C", "O"], "coords": [[0, 0, 0], [0, 0, 1.13]], 
           "adsorption_height_A": 1.8, "binding": "atop", "orientation": "vertical"},
    "NO": {"atoms": ["N", "O"], "coords": [[0, 0, 0], [0, 0, 1.15]], 
           "adsorption_height_A": 1.7, "binding": "bridge", "orientation": "tilted"},
    "O2": {"atoms": ["O", "O"], "coords": [[0, 0, 0], [0, 0, 1.21]], 
           "adsorption_height_A": 1.5, "binding": "bridge", "orientation": "horizontal"},
    "N2": {"atoms": ["N", "N"], "coords": [[0, 0, 0], [0, 0, 1.10]], 
           "adsorption_height_A": 2.0, "binding": "atop", "orientation": "vertical"},
    "H2": {"atoms": ["H", "H"], "coords": [[0, 0, 0], [0, 0, 0.74]], 
           "adsorption_height_A": 2.5, "binding": "physisorbed"},
    
    # Triatomic molecules
    "H2O": {"atoms": ["O", "H", "H"], 
            "coords": [[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]], 
            "adsorption_height_A": 2.5, "binding": "atop_O", "orientation": "flat"},
    "CO2": {"atoms": ["C", "O", "O"], 
            "coords": [[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]], 
            "adsorption_height_A": 3.0, "binding": "physisorbed", "orientation": "horizontal"},
    "NO2": {"atoms": ["N", "O", "O"], 
            "coords": [[0, 0, 0], [1.19, 0.46, 0], [-1.19, 0.46, 0]], 
            "adsorption_height_A": 2.0, "binding": "O_down"},
    "H2S": {"atoms": ["S", "H", "H"], 
            "coords": [[0, 0, 0], [0.96, 0.59, 0], [-0.96, 0.59, 0]], 
            "adsorption_height_A": 2.2, "binding": "atop_S"},
    "SO2": {"atoms": ["S", "O", "O"], 
            "coords": [[0, 0, 0], [1.43, 0.52, 0], [-1.43, 0.52, 0]], 
            "adsorption_height_A": 2.5, "binding": "bridge"},
    
    # Larger molecules
    "NH3": {"atoms": ["N", "H", "H", "H"], 
            "coords": [[0, 0, 0], [0.94, 0.38, 0], [-0.47, 0.38, 0.81], [-0.47, 0.38, -0.81]], 
            "adsorption_height_A": 2.0, "binding": "atop_N", "orientation": "N_down"},
    "CH4": {"atoms": ["C", "H", "H", "H", "H"], 
            "coords": [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], 
                      [0.63, -0.63, -0.63], [-0.63, 0.63, -0.63]], 
            "adsorption_height_A": 3.5, "binding": "physisorbed"},
    "C2H4": {"atoms": ["C", "C", "H", "H", "H", "H"], 
             "coords": [[0, 0, 0.67], [0, 0, -0.67], [0, 0.92, 1.24], [0, -0.92, 1.24],
                       [0, 0.92, -1.24], [0, -0.92, -1.24]], 
             "adsorption_height_A": 2.0, "binding": "di-sigma", "orientation": "flat"},
    "C2H2": {"atoms": ["C", "C", "H", "H"], 
             "coords": [[0, 0, 0.60], [0, 0, -0.60], [0, 0, 1.67], [0, 0, -1.67]], 
             "adsorption_height_A": 2.0, "binding": "bridge"},
    "CH3OH": {"atoms": ["C", "O", "H", "H", "H", "H"], 
              "coords": [[0, 0, 0], [1.43, 0, 0], [-0.38, 1.03, 0], [-0.38, -0.51, 0.89],
                        [-0.38, -0.51, -0.89], [1.80, 0, 0.87]], 
              "adsorption_height_A": 2.5, "binding": "O_down"},
    "HCOOH": {"atoms": ["C", "O", "O", "H", "H"], 
              "coords": [[0, 0, 0], [1.20, 0.25, 0], [-0.30, 1.20, 0], [-0.80, -0.75, 0], [1.40, 1.15, 0]], 
              "adsorption_height_A": 2.5, "binding": "bidentate"},
    
    # Reaction intermediates
    "OH": {"atoms": ["O", "H"], "coords": [[0, 0, 0], [0, 0, 0.97]], 
           "adsorption_height_A": 1.5, "binding": "bridge"},
    "OOH": {"atoms": ["O", "O", "H"], 
            "coords": [[0, 0, 0], [1.30, 0.50, 0], [1.80, 1.30, 0]], 
            "adsorption_height_A": 1.8, "binding": "atop"},
    "CHO": {"atoms": ["C", "H", "O"], 
            "coords": [[0, 0, 0], [-0.60, 0.90, 0], [1.20, 0, 0]], 
            "adsorption_height_A": 1.8, "binding": "atop_C"},
    "COOH": {"atoms": ["C", "O", "O", "H"], 
             "coords": [[0, 0, 0], [1.20, 0.25, 0], [-0.30, 1.20, 0], [1.40, 1.15, 0]], 
             "adsorption_height_A": 2.0, "binding": "bidentate"},
}


# Coverage patterns database
COVERAGE_PATTERNS = {
    "1x1": {"coverage_ML": 1.0, "pattern": [[0, 0]], "symmetry": "p1"},
    "2x2_quarter": {"coverage_ML": 0.25, "pattern": [[0, 0]], "symmetry": "p4mm"},
    "2x2_half": {"coverage_ML": 0.50, "pattern": [[0, 0], [0.5, 0.5]], "symmetry": "p2mm"},
    "sqrt3_sqrt3": {"coverage_ML": 0.33, "pattern": [[0, 0]], "symmetry": "p6mm"},
    "c2x2": {"coverage_ML": 0.50, "pattern": [[0, 0]], "symmetry": "c2mm"},
    "2x1": {"coverage_ML": 0.50, "pattern": [[0, 0]], "symmetry": "p2mm"},
    "3x3_R30": {"coverage_ML": 0.11, "pattern": [[0, 0]], "symmetry": "p6mm"},
    "honeycomb": {"coverage_ML": 0.33, "pattern": [[0.33, 0.33], [0.67, 0.67]], "symmetry": "p6mm"},
    "kagome": {"coverage_ML": 0.75, "pattern": [[0.25, 0], [0.75, 0], [0.5, 0.5]], "symmetry": "p6mm"},
}


# Adsorption site types
ADSORPTION_SITES = {
    "atop": {"coordination": 1, "offset": [0, 0]},
    "bridge": {"coordination": 2, "offset": [0.25, 0]},
    "fcc_hollow": {"coordination": 3, "offset": [0.33, 0.33]},
    "hcp_hollow": {"coordination": 3, "offset": [0.67, 0.33]},
    "4fold_hollow": {"coordination": 4, "offset": [0.5, 0.5]},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_adsorbate_on_surface(
    surface_structure: Dict[str, Any] = None,
    adsorbate: str = "CO",
    site: str = "atop",
    coverage: str = "2x2_quarter",
    height_offset: float = 0.0,
    rotation_deg: float = 0.0,
    **kwargs  # Accept structure, site_type, height as aliases
) -> Dict[str, Any]:
    """
    Place adsorbate molecule on surface.

    Args:
        surface_structure: Surface slab structure dict
        adsorbate: Adsorbate molecule from database
        site: Adsorption site type
        coverage: Coverage pattern
        height_offset: Additional height offset (Å)
        rotation_deg: Rotation around surface normal
        **kwargs: Accept aliases (structure, site_type, height)

    Returns:
        Surface with adsorbate
    """
    # Handle parameter aliases
    if surface_structure is None:
        surface_structure = kwargs.get('structure')
    if surface_structure is None:
        return {"success": False, "error": {"code": "MISSING_STRUCTURE",
                "message": "Must provide 'surface_structure' or 'structure'"}}

    # Handle pymatgen Structure objects by converting to dict
    from pymatgen.core import Structure
    if isinstance(surface_structure, Structure):
        lattice = surface_structure.lattice
        surface_structure = {
            "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                        "matrix": lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in surface_structure]
        }

    # Handle site_type alias
    if 'site_type' in kwargs:
        site_raw = kwargs['site_type'].lower()
        # Map common names
        site_map = {'ontop': 'atop', 'on_top': 'atop', 'top': 'atop'}
        site = site_map.get(site_raw, site_raw)

    # Handle height alias
    if 'height' in kwargs and height_offset == 0.0:
        height_offset = kwargs['height']

    if adsorbate not in ADSORBATE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ADSORBATE", "message": f"Unknown adsorbate '{adsorbate}'",
                      "available": list(ADSORBATE_DATABASE.keys())}
        }
    
    if site not in ADSORPTION_SITES:
        return {
            "success": False,
            "error": {"code": "INVALID_SITE", "message": f"Unknown site '{site}'",
                      "available": list(ADSORPTION_SITES.keys())}
        }
    
    if coverage not in COVERAGE_PATTERNS:
        return {
            "success": False,
            "error": {"code": "INVALID_COVERAGE", "message": f"Unknown coverage '{coverage}'",
                      "available": list(COVERAGE_PATTERNS.keys())}
        }
    
    ads_info = ADSORBATE_DATABASE[adsorbate]
    site_info = ADSORPTION_SITES[site]
    cov_info = COVERAGE_PATTERNS[coverage]
    
    # Get surface atoms
    surface_atoms = surface_structure.get("atoms", [])
    lattice_info = surface_structure.get("lattice", {})
    
    a = lattice_info.get("a", 10)
    b = lattice_info.get("b", 10)
    c = lattice_info.get("c", 30)
    
    # Find surface z position
    z_coords = [atom["coords"][2] for atom in surface_atoms]
    z_max = max(z_coords) if z_coords else 0.5
    
    # Calculate adsorption height in fractional coords
    ads_height = (ads_info["adsorption_height_A"] + height_offset) / c
    z_ads = z_max + ads_height
    
    # Apply rotation
    rotation_rad = np.radians(rotation_deg)
    
    # Add adsorbate at each coverage position
    new_atoms = list(surface_atoms)
    n_adsorbates = 0
    
    for pattern_pos in cov_info["pattern"]:
        # Calculate position
        x_base = pattern_pos[0] + site_info["offset"][0]
        y_base = pattern_pos[1] + site_info["offset"][1]
        
        # Add all atoms of adsorbate
        for i, (elem, coord) in enumerate(zip(ads_info["atoms"], ads_info["coords"])):
            # Apply rotation around z-axis
            x_rot = coord[0] * np.cos(rotation_rad) - coord[1] * np.sin(rotation_rad)
            y_rot = coord[0] * np.sin(rotation_rad) + coord[1] * np.cos(rotation_rad)
            
            x = x_base + x_rot / a
            y = y_base + y_rot / b
            z = z_ads + coord[2] / c
            
            new_atoms.append({
                "element": elem,
                "coords": [x % 1, y % 1, z],
                "is_adsorbate": True
            })
        
        n_adsorbates += 1
    
    return {
        "success": True,
        "adsorbate": adsorbate,
        "n_adsorbate_atoms": len(ads_info["atoms"]),
        "adsorption_site": site,
        "coverage_pattern": coverage,
        "coverage_ML": cov_info["coverage_ML"],
        "n_adsorbate_molecules": n_adsorbates,
        "binding_mode": ads_info.get("binding", "unknown"),
        "orientation": ads_info.get("orientation", "unknown"),
        "adsorption_height_A": ads_info["adsorption_height_A"] + height_offset,
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms, "lattice": lattice_info}
    }


def generate_coadsorption(
    surface_structure: Dict[str, Any],
    adsorbates: List[str] = ["CO", "O"],
    sites: List[str] = ["atop", "fcc_hollow"],
    positions: List[Tuple[float, float]] = [(0.25, 0.25), (0.75, 0.75)]
) -> Dict[str, Any]:
    """
    Generate co-adsorption configuration.
    
    Args:
        surface_structure: Surface slab structure
        adsorbates: List of adsorbate molecules
        sites: Adsorption sites for each
        positions: (x, y) positions for each
    
    Returns:
        Surface with multiple adsorbates
    """
    if len(adsorbates) != len(sites) or len(adsorbates) != len(positions):
        return {
            "success": False,
            "error": {"code": "LENGTH_MISMATCH", "message": "adsorbates, sites, positions must have same length"}
        }
    
    surface_atoms = list(surface_structure.get("atoms", []))
    lattice_info = surface_structure.get("lattice", {})
    
    c = lattice_info.get("c", 30)
    z_coords = [atom["coords"][2] for atom in surface_atoms]
    z_max = max(z_coords) if z_coords else 0.5
    
    new_atoms = list(surface_atoms)
    
    for adsorbate, site, (x_pos, y_pos) in zip(adsorbates, sites, positions):
        if adsorbate not in ADSORBATE_DATABASE:
            continue
        
        ads_info = ADSORBATE_DATABASE[adsorbate]
        site_info = ADSORPTION_SITES.get(site, {"offset": [0, 0]})
        
        ads_height = ads_info["adsorption_height_A"] / c
        z_ads = z_max + ads_height
        
        x_base = x_pos + site_info["offset"][0]
        y_base = y_pos + site_info["offset"][1]
        
        for elem, coord in zip(ads_info["atoms"], ads_info["coords"]):
            a = lattice_info.get("a", 10)
            b = lattice_info.get("b", 10)
            
            x = x_base + coord[0] / a
            y = y_base + coord[1] / b
            z = z_ads + coord[2] / c
            
            new_atoms.append({
                "element": elem,
                "coords": [x % 1, y % 1, z],
                "is_adsorbate": True,
                "adsorbate_type": adsorbate
            })
    
    return {
        "success": True,
        "adsorbates": adsorbates,
        "sites": sites,
        "positions": positions,
        "is_coadsorption": True,
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms, "lattice": lattice_info}
    }


def get_adsorbate_molecule(adsorbate: str) -> Dict[str, Any]:
    """Get adsorbate molecule information."""
    if adsorbate not in ADSORBATE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ADSORBATE", "message": f"Unknown adsorbate",
                      "available": list(ADSORBATE_DATABASE.keys())}
        }
    
    info = ADSORBATE_DATABASE[adsorbate]
    return {
        "success": True,
        "adsorbate": adsorbate,
        "atoms": info["atoms"],
        "coords": info["coords"],
        "n_atoms": len(info["atoms"]),
        "adsorption_height_A": info["adsorption_height_A"],
        "preferred_binding": info.get("binding", "unknown"),
        "orientation": info.get("orientation", "unknown")
    }
