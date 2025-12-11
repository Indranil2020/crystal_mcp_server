"""
molecule/cages.py - Cage and Fullerene Generation

Generates cage molecular structures:
- Fullerenes (C60, C70, C80, ...)
- Endohedral fullerenes (M@C60)
- Metal-organic cages
- Cucurbiturils
"""

from typing import Dict, Any, List, Optional
import numpy as np


FULLERENE_DATABASE = {
    "C60": {"n_atoms": 60, "diameter": 7.1, "symmetry": "Ih"},
    "C70": {"n_atoms": 70, "diameter": 7.8, "symmetry": "D5h"},
    "C80": {"n_atoms": 80, "diameter": 8.2, "symmetry": "D5d"},
    "C84": {"n_atoms": 84, "diameter": 8.4, "symmetry": "D2d"},
    "C540": {"n_atoms": 540, "diameter": 14.0, "symmetry": "Ih"},
}


def generate_fullerene(
    formula: str = "C60",
    endohedral: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate fullerene structure.
    
    Args:
        formula: Fullerene formula (C60, C70, etc.)
        endohedral: Endohedral guest atom (e.g., "La" for La@C60)
    
    Returns:
        Fullerene structure
    """
    if formula not in FULLERENE_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown fullerene: {formula}"}}
    
    info = FULLERENE_DATABASE[formula]
    n_atoms = info["n_atoms"]
    radius = info["diameter"] / 2
    
    # Generate approximate positions on sphere
    positions = []
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    
    # Fibonacci sphere for approximately uniform distribution
    for i in range(n_atoms):
        y = 1 - (i / (n_atoms - 1)) * 2  # -1 to 1
        r = np.sqrt(1 - y * y)
        theta = 2 * np.pi * i / phi
        
        x = np.cos(theta) * r * radius
        z = np.sin(theta) * r * radius
        y_scaled = y * radius
        
        positions.append([x, y_scaled, z])
    
    atoms = [{"element": "C", "cartesian": list(p)} for p in positions]
    
    # Add endohedral guest
    if endohedral:
        atoms.append({"element": endohedral, "cartesian": [0, 0, 0]})
    
    name = f"{endohedral}@{formula}" if endohedral else formula
    
    return {
        "success": True,
        "formula": name,
        "n_carbon": n_atoms,
        "diameter_angstrom": info["diameter"],
        "symmetry": info["symmetry"],
        "endohedral": endohedral,
        "structure": {"atoms": atoms}
    }


def generate_cage(
    cage_type: str = "M12L24",
    metal: str = "Pd"
) -> Dict[str, Any]:
    """
    Generate metal-organic cage structure.
    
    Args:
        cage_type: Cage topology (M12L24, M24L48, Fe4L6)
        metal: Metal element
    
    Returns:
        Cage structure
    """
    cage_types = {
        "M12L24": {"n_metal": 12, "n_ligand": 24, "radius": 15.0},
        "M24L48": {"n_metal": 24, "n_ligand": 48, "radius": 20.0},
        "Fe4L6": {"n_metal": 4, "n_ligand": 6, "radius": 8.0},
        "Pd6L4": {"n_metal": 6, "n_ligand": 4, "radius": 10.0},
    }
    
    if cage_type not in cage_types:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown cage type"}}
    
    info = cage_types[cage_type]
    n_metal = info["n_metal"]
    radius = info["radius"]
    
    # Generate metal positions on sphere
    positions = []
    phi = (1 + np.sqrt(5)) / 2
    
    for i in range(n_metal):
        y = 1 - (i / (n_metal - 1)) * 2 if n_metal > 1 else 0
        r = np.sqrt(max(0, 1 - y * y))
        theta = 2 * np.pi * i / phi
        
        x = np.cos(theta) * r * radius
        z = np.sin(theta) * r * radius
        y_scaled = y * radius
        
        positions.append([x, y_scaled, z])
    
    atoms = [{"element": metal, "cartesian": list(p)} for p in positions]
    
    return {
        "success": True,
        "cage_type": cage_type,
        "metal": metal,
        "n_metal": n_metal,
        "n_ligand": info["n_ligand"],
        "cavity_diameter_angstrom": radius * 2 * 0.8,  # Approximate
        "structure": {"atoms": atoms}
    }
