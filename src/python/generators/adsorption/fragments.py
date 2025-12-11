"""
adsorption/fragments.py - Dissociated Adsorbate Fragments

Generates dissociated fragment adsorbates:
- OH*, O*, H* radicals
- CHx*, NHx*, SH* fragments
- CxHyOz* intermediates
- Coverage control from 0.01 to 2 ML
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Fragment database with adsorption heights
FRAGMENT_DATABASE = {
    # Atomic adsorbates
    "O*": {"elements": ["O"], "height": 1.2, "charge": -1},
    "H*": {"elements": ["H"], "height": 1.0, "charge": 0},
    "N*": {"elements": ["N"], "height": 1.3, "charge": 0},
    "C*": {"elements": ["C"], "height": 1.4, "charge": 0},
    "S*": {"elements": ["S"], "height": 1.6, "charge": 0},
    
    # Hydroxyl and related
    "OH*": {"elements": ["O", "H"], "positions": [[0, 0, 0], [0, 0, 0.97]], "height": 1.2},
    "OOH*": {"elements": ["O", "O", "H"], "positions": [[0, 0, 0], [0, 0, 1.4], [0.8, 0, 1.8]], "height": 1.3},
    
    # CHx fragments
    "CH*": {"elements": ["C", "H"], "positions": [[0, 0, 0], [0, 0, 1.09]], "height": 1.3},
    "CH2*": {"elements": ["C", "H", "H"], "positions": [[0, 0, 0], [0.94, 0, 0.5], [-0.94, 0, 0.5]], "height": 1.3},
    "CH3*": {"elements": ["C", "H", "H", "H"], 
             "positions": [[0, 0, 0], [0.95, 0.55, 0], [-0.95, 0.55, 0], [0, -1.09, 0]], "height": 1.5},
    
    # NHx fragments
    "NH*": {"elements": ["N", "H"], "positions": [[0, 0, 0], [0, 0, 1.01]], "height": 1.3},
    "NH2*": {"elements": ["N", "H", "H"], "positions": [[0, 0, 0], [0.8, 0.6, 0], [-0.8, 0.6, 0]], "height": 1.3},
    
    # Carbon-oxygen fragments
    "CO*": {"elements": ["C", "O"], "positions": [[0, 0, 0], [0, 0, 1.14]], "height": 1.5, "orientation": "vertical"},
    "COH*": {"elements": ["C", "O", "H"], "positions": [[0, 0, 0], [0, 0, 1.3], [0, 0.8, 1.7]], "height": 1.5},
    "CHO*": {"elements": ["C", "H", "O"], "positions": [[0, 0, 0], [0.9, 0, 0.5], [0, 0, 1.2]], "height": 1.4},
    "COOH*": {"elements": ["C", "O", "O", "H"], 
              "positions": [[0, 0, 0], [0, 1.0, 0.6], [0, -0.6, 1.2], [0, -0.6, 2.1]], "height": 1.5},
    
    # Formate and related
    "HCOO*": {"elements": ["H", "C", "O", "O"], 
              "positions": [[0, 0, 0], [0, 0, 1.1], [0.9, 0, 1.6], [-0.9, 0, 1.6]], "height": 1.3},
    
    # SH fragment
    "SH*": {"elements": ["S", "H"], "positions": [[0, 0, 0], [0, 0, 1.34]], "height": 1.6},
}


def generate_fragment_adsorbate(
    surface: Dict[str, Any],
    fragment: str = "OH*",
    site: str = "top",
    coverage: float = 0.25,
    site_index: Optional[int] = None
) -> Dict[str, Any]:
    """
    Add dissociated fragment to surface.
    
    Args:
        surface: Surface structure dict
        fragment: Fragment type from database
        site: Adsorption site (top, bridge, hollow, fcc, hcp)
        coverage: Coverage in ML (0.01-2.0)
        site_index: Specific site index (for precise placement)
    
    Returns:
        Surface with fragment adsorbate
    """
    if fragment not in FRAGMENT_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_FRAGMENT", "message": f"Unknown fragment '{fragment}'",
                      "available": list(FRAGMENT_DATABASE.keys())}
        }
    
    frag_info = FRAGMENT_DATABASE[fragment]
    elements = frag_info["elements"]
    height = frag_info["height"]
    
    # Get positions (default is single atom at origin)
    positions = frag_info.get("positions", [[0, 0, 0]])
    
    # Get surface atoms
    surface_atoms = surface.get("structure", {}).get("atoms", [])
    if not surface_atoms:
        return {"success": False, "error": {"code": "NO_SURFACE", "message": "Invalid surface"}}
    
    # Find top surface layer
    z_coords = [a.get("cartesian", a.get("coords", [0,0,0]))[2] for a in surface_atoms]
    z_max = max(z_coords)
    
    # Find surface sites
    surface_layer = [a for a in surface_atoms if abs(a.get("cartesian", a.get("coords", [0,0,0]))[2] - z_max) < 0.5]
    
    # Calculate number of fragments based on coverage
    n_surface = len(surface_layer)
    n_fragments = max(1, int(n_surface * coverage))
    
    # Select sites
    if site_index is not None:
        selected_sites = [surface_layer[site_index % len(surface_layer)]]
    else:
        np.random.seed(42)
        indices = np.random.choice(len(surface_layer), min(n_fragments, len(surface_layer)), replace=False)
        selected_sites = [surface_layer[i] for i in indices]
    
    # Create new atoms list
    new_atoms = list(surface_atoms)
    
    for surf_atom in selected_sites:
        base_pos = surf_atom.get("cartesian", surf_atom.get("coords", [0, 0, 0]))
        
        # Site-specific offsets
        site_offsets = {
            "top": [0, 0],
            "bridge": [0.7, 0],
            "hollow": [0.7, 0.7],
            "fcc": [0.7, 0.4],
            "hcp": [0.4, 0.7],
        }
        offset = site_offsets.get(site, [0, 0])
        
        # Add fragment atoms
        for elem, pos in zip(elements, positions):
            new_pos = [
                base_pos[0] + pos[0] + offset[0],
                base_pos[1] + pos[1] + offset[1],
                z_max + height + pos[2]
            ]
            new_atoms.append({"element": elem, "cartesian": new_pos, "adsorbate": True, "fragment": fragment})
    
    return {
        "success": True,
        "fragment": fragment,
        "site": site,
        "coverage_ML": coverage,
        "n_fragments": len(selected_sites),
        "fragment_composition": elements,
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms}
    }


def generate_coverage_pattern(
    surface: Dict[str, Any],
    adsorbate: str = "O*",
    pattern: str = "p(2x2)",
    height: float = 1.2
) -> Dict[str, Any]:
    """
    Generate ordered adsorbate coverage pattern.
    
    Args:
        surface: Surface structure
        adsorbate: Adsorbate species
        pattern: Coverage pattern (p(2x2), c(4x2), (√3×√3)R30°, etc.)
        height: Adsorbate height
    
    Returns:
        Surface with ordered adsorbate pattern
    """
    COVERAGE_PATTERNS = {
        "p(1x1)": {"coverage": 1.0, "positions": [[0.5, 0.5]]},
        "p(2x2)": {"coverage": 0.25, "positions": [[0.5, 0.5]]},
        "p(2x1)": {"coverage": 0.5, "positions": [[0.5, 0.5]]},
        "c(2x2)": {"coverage": 0.5, "positions": [[0, 0.5], [0.5, 0]]},
        "c(4x2)": {"coverage": 0.25, "positions": [[0.25, 0.5], [0.75, 0]]},
        "(sqrt3xsqrt3)R30": {"coverage": 0.33, "positions": [[1/3, 1/3], [2/3, 2/3]]},
        "p(3x3)": {"coverage": 0.11, "positions": [[0.5, 0.5]]},
        "p(3x1)": {"coverage": 0.33, "positions": [[0.5, 0.5]]},
        "p(4x4)": {"coverage": 0.0625, "positions": [[0.5, 0.5]]},
    }
    
    if pattern not in COVERAGE_PATTERNS:
        return {
            "success": False,
            "error": {"code": "INVALID_PATTERN", "message": f"Unknown pattern '{pattern}'",
                      "available": list(COVERAGE_PATTERNS.keys())}
        }
    
    pat_info = COVERAGE_PATTERNS[pattern]
    
    # Get surface info
    surface_atoms = surface.get("structure", {}).get("atoms", [])
    z_coords = [a.get("cartesian", [0,0,0])[2] for a in surface_atoms]
    z_max = max(z_coords) if z_coords else 0
    
    # Approximate surface dimensions
    x_coords = [a.get("cartesian", [0,0,0])[0] for a in surface_atoms]
    y_coords = [a.get("cartesian", [0,0,0])[1] for a in surface_atoms]
    x_span = max(x_coords) - min(x_coords) if x_coords else 10
    y_span = max(y_coords) - min(y_coords) if y_coords else 10
    
    # Get adsorbate elements
    if adsorbate in FRAGMENT_DATABASE:
        elements = FRAGMENT_DATABASE[adsorbate]["elements"]
    else:
        elements = [adsorbate.replace("*", "")]
    
    new_atoms = list(surface_atoms)
    
    # Add adsorbates at pattern positions
    for frac_pos in pat_info["positions"]:
        x = min(x_coords) + frac_pos[0] * x_span if x_coords else frac_pos[0] * 10
        y = min(y_coords) + frac_pos[1] * y_span if y_coords else frac_pos[1] * 10
        z = z_max + height
        
        for elem in elements:
            new_atoms.append({"element": elem, "cartesian": [x, y, z], "adsorbate": True})
    
    return {
        "success": True,
        "adsorbate": adsorbate,
        "pattern": pattern,
        "coverage_ML": pat_info["coverage"],
        "n_adsorbates": len(pat_info["positions"]),
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms}
    }
