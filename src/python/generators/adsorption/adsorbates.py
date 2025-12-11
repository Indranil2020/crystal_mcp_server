"""
adsorption/adsorbates.py - Adsorbate Addition

Adds molecular adsorbates to surface structures.
"""

from typing import Dict, Any, List, Optional
import numpy as np


ADSORPTION_SITES = {
    "top": {"description": "On top of surface atom"},
    "bridge": {"description": "Between two surface atoms"},
    "hollow": {"description": "Central position in surface unit"},
    "fcc": {"description": "FCC hollow site"},
    "hcp": {"description": "HCP hollow site"},
}


def add_adsorbate(
    surface_dict: Dict[str, Any],
    adsorbate: str = "CO",
    site: str = "top",
    height: float = 2.0,
    position: List[float] = [0.5, 0.5]
) -> Dict[str, Any]:
    """
    Add adsorbate to surface structure.
    
    Args:
        surface_dict: Surface structure dictionary
        adsorbate: Adsorbate molecule
        site: Adsorption site type
        height: Adsorbate height above surface
        position: (x, y) fractional position
    
    Returns:
        Surface with adsorbate
    """
    # Get surface atoms
    atoms = list(surface_dict.get("atoms", []))
    lattice = surface_dict.get("lattice", {})
    c = lattice.get("c", 20.0)
    
    # Find top z position
    z_max = max(a["coords"][2] for a in atoms) if atoms else 0.5
    z_ads = z_max + height / c
    
    # Add adsorbate atoms
    adsorbate_data = {
        "CO": [("C", 0), ("O", 1.13)],
        "H2O": [("O", 0), ("H", 0.96), ("H", 0.96)],
        "O": [("O", 0)],
        "H": [("H", 0)],
        "N": [("N", 0)],
        "OH": [("O", 0), ("H", 0.97)],
    }
    
    mol = adsorbate_data.get(adsorbate, [("X", 0)])
    
    for elem, dz in mol:
        atoms.append({
            "element": elem,
            "coords": [position[0], position[1], z_ads + dz / c]
        })
    
    result = dict(surface_dict)
    result["atoms"] = atoms
    result["adsorbate"] = adsorbate
    result["adsorption_site"] = site
    
    return {"success": True, "structure": result}
