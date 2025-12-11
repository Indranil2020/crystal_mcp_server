"""
adsorption/ionic_liquids.py - Ionic Liquid Monolayers

Generates ionic liquid adlayers:
- [BMIM][PF6] on graphite/graphene
- [EMIM][BF4] on various surfaces
- Various cation/anion combinations
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Ionic liquid database
IL_DATABASE = {
    # Imidazolium-based ILs
    "[BMIM][PF6]": {
        "cation": "BMIM", "anion": "PF6",
        "cation_size": 6.5, "anion_size": 3.2,
        "description": "1-Butyl-3-methylimidazolium hexafluorophosphate"
    },
    "[BMIM][BF4]": {
        "cation": "BMIM", "anion": "BF4",
        "cation_size": 6.5, "anion_size": 2.8,
        "description": "1-Butyl-3-methylimidazolium tetrafluoroborate"
    },
    "[BMIM][NTf2]": {
        "cation": "BMIM", "anion": "NTf2",
        "cation_size": 6.5, "anion_size": 5.5,
        "description": "1-Butyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide"
    },
    "[EMIM][BF4]": {
        "cation": "EMIM", "anion": "BF4",
        "cation_size": 5.5, "anion_size": 2.8,
        "description": "1-Ethyl-3-methylimidazolium tetrafluoroborate"
    },
    "[EMIM][OTf]": {
        "cation": "EMIM", "anion": "OTf",
        "cation_size": 5.5, "anion_size": 4.0,
        "description": "1-Ethyl-3-methylimidazolium trifluoromethanesulfonate"
    },
    # Pyrrolidinium-based
    "[P14][NTf2]": {
        "cation": "P14", "anion": "NTf2",
        "cation_size": 6.0, "anion_size": 5.5,
        "description": "1-Butyl-1-methylpyrrolidinium bis(trifluoromethylsulfonyl)imide"
    },
}


def generate_imidazolium_cation(center: List[float], rotation: float = 0) -> List[Dict[str, Any]]:
    """Generate simplified imidazolium cation atoms."""
    atoms = []
    r = 1.35  # Ring radius
    
    # Imidazole ring (5 atoms: 2N, 3C)
    ring_elements = ["N", "C", "N", "C", "C"]
    for i, elem in enumerate(ring_elements):
        angle = 2 * np.pi * i / 5 + rotation
        x = center[0] + r * np.cos(angle)
        y = center[1] + r * np.sin(angle)
        atoms.append({"element": elem, "cartesian": [x, y, center[2]]})
    
    # Methyl group
    atoms.append({"element": "C", "cartesian": [center[0] - 2.5, center[1], center[2]]})
    
    # Butyl chain (4 C atoms)
    for i in range(4):
        atoms.append({"element": "C", "cartesian": [center[0] + 2.5 + i * 1.5, center[1], center[2]]})
    
    return atoms


def generate_anion(anion_type: str, center: List[float]) -> List[Dict[str, Any]]:
    """Generate simplified anion atoms."""
    atoms = []
    
    if anion_type == "PF6":
        atoms.append({"element": "P", "cartesian": center})
        for dx, dy, dz in [[1.6,0,0], [-1.6,0,0], [0,1.6,0], [0,-1.6,0], [0,0,1.6], [0,0,-1.6]]:
            atoms.append({"element": "F", "cartesian": [center[0]+dx, center[1]+dy, center[2]+dz]})
    elif anion_type == "BF4":
        atoms.append({"element": "B", "cartesian": center})
        for i in range(4):
            angle = 2 * np.pi * i / 4
            atoms.append({"element": "F", "cartesian": [center[0] + 1.3*np.cos(angle), center[1] + 1.3*np.sin(angle), center[2]]})
    elif anion_type == "NTf2":
        atoms.append({"element": "N", "cartesian": center})
        atoms.append({"element": "S", "cartesian": [center[0] - 1.5, center[1], center[2]]})
        atoms.append({"element": "S", "cartesian": [center[0] + 1.5, center[1], center[2]]})
    elif anion_type == "OTf":
        atoms.append({"element": "S", "cartesian": center})
        for i in range(3):
            angle = 2 * np.pi * i / 3
            atoms.append({"element": "O", "cartesian": [center[0] + 1.5*np.cos(angle), center[1] + 1.5*np.sin(angle), center[2]]})
    
    return atoms


def generate_il_monolayer(
    il_type: str = "[BMIM][PF6]",
    substrate: str = "graphite",
    coverage: float = 1.0,
    size: List[int] = [4, 4],
    orientation: str = "parallel"
) -> Dict[str, Any]:
    """
    Generate ionic liquid monolayer on surface.
    
    Args:
        il_type: Ionic liquid from database
        substrate: Surface type (graphite, graphene, Au)
        coverage: Coverage (0-1)
        size: Surface size
        orientation: Cation orientation (parallel, perpendicular)
    
    Returns:
        IL monolayer structure
    """
    if il_type not in IL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_IL", "message": f"Unknown IL '{il_type}'",
                      "available": list(IL_DATABASE.keys())}
        }
    
    info = IL_DATABASE[il_type]
    cation_size = info["cation_size"]
    anion_size = info["anion_size"]
    anion_type = info["anion"]
    
    nx, ny = size
    atoms = []
    n_pairs = 0
    
    np.random.seed(42)
    
    # Generate substrate
    spacing = 4.26  # Graphite lattice
    for i in range(nx + 2):
        for j in range(ny + 2):
            x = i * spacing
            y = j * spacing * np.sqrt(3) / 2
            if substrate == "graphite":
                atoms.append({"element": "C", "cartesian": [x, y, 0], "layer": "substrate"})
                atoms.append({"element": "C", "cartesian": [x, y, -3.35], "layer": "substrate"})
            elif substrate == "graphene":
                atoms.append({"element": "C", "cartesian": [x, y, 0], "layer": "substrate"})
            else:  # Au
                atoms.append({"element": "Au", "cartesian": [x, y, 0], "layer": "substrate"})
    
    # Generate IL pairs
    il_spacing = cation_size + anion_size
    z_cation = 3.5 if orientation == "parallel" else 5.0
    z_anion = 3.5
    
    for i in range(nx):
        for j in range(ny):
            if np.random.random() > coverage:
                continue
            
            x = i * il_spacing + il_spacing / 2
            y = j * il_spacing + (i % 2) * il_spacing / 2
            
            # Add cation
            rotation = np.random.uniform(0, 2*np.pi)
            cation_atoms = generate_imidazolium_cation([x, y, z_cation], rotation)
            for a in cation_atoms:
                a["ion"] = "cation"
            atoms.extend(cation_atoms)
            
            # Add anion (offset from cation)
            anion_x = x + cation_size * 0.7
            anion_y = y
            anion_atoms = generate_anion(anion_type, [anion_x, anion_y, z_anion])
            for a in anion_atoms:
                a["ion"] = "anion"
            atoms.extend(anion_atoms)
            
            n_pairs += 1
    
    return {
        "success": True,
        "ionic_liquid": il_type,
        "description": info["description"],
        "substrate": substrate,
        "coverage": coverage,
        "orientation": orientation,
        "n_ion_pairs": n_pairs,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
