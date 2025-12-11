"""
molecule/porphyrins.py - Porphyrins and Phthalocyanines

Generates cyclic tetrapyrrole structures:
- Porphyrins (with various metals)
- Phthalocyanines
- Corroles
- Sub-phthalocyanines
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Metal-porphyrin database
PORPHYRIN_DATABASE = {
    "H2-porphyrin": {"metal": None, "beta_substituents": None, "meso_substituents": None},
    "Fe-porphyrin": {"metal": "Fe", "spin": 2, "oxidation": 2},
    "Co-porphyrin": {"metal": "Co", "spin": 0.5, "oxidation": 2},
    "Ni-porphyrin": {"metal": "Ni", "spin": 0, "oxidation": 2},
    "Cu-porphyrin": {"metal": "Cu", "spin": 0.5, "oxidation": 2},
    "Zn-porphyrin": {"metal": "Zn", "spin": 0, "oxidation": 2},
    "Mg-porphyrin": {"metal": "Mg", "spin": 0, "oxidation": 2},
    "Mn-porphyrin": {"metal": "Mn", "spin": 2.5, "oxidation": 2},
    "Pd-porphyrin": {"metal": "Pd", "spin": 0, "oxidation": 2},
    "Pt-porphyrin": {"metal": "Pt", "spin": 0, "oxidation": 2},
    "Ru-porphyrin": {"metal": "Ru", "spin": 0, "oxidation": 2},
    "chlorophyll-a": {"metal": "Mg", "special": True, "description": "Chlorophyll a"},
    "heme": {"metal": "Fe", "special": True, "description": "Heme b"},
    "vitamin-B12": {"metal": "Co", "special": True, "description": "Cobalamin core"},
}


# Phthalocyanine database
PHTHALOCYANINE_DATABASE = {
    "H2-Pc": {"metal": None, "substituents": None},
    "Cu-Pc": {"metal": "Cu", "color": "blue"},
    "Zn-Pc": {"metal": "Zn", "color": "blue-green"},
    "Fe-Pc": {"metal": "Fe", "spin": 1},
    "Co-Pc": {"metal": "Co", "spin": 0.5},
    "Ni-Pc": {"metal": "Ni"},
    "Mg-Pc": {"metal": "Mg"},
    "Al-Pc-Cl": {"metal": "Al", "axial": "Cl"},
    "Si-Pc-Cl2": {"metal": "Si", "axial": ["Cl", "Cl"]},
    "Ti-Pc-O": {"metal": "Ti", "axial": "O"},
    "V-Pc-O": {"metal": "V", "axial": "O"},
}


# Corrole database
CORROLE_DATABASE = {
    "H3-corrole": {"metal": None},
    "Fe-corrole": {"metal": "Fe", "oxidation": 3},
    "Co-corrole": {"metal": "Co", "oxidation": 3},
    "Mn-corrole": {"metal": "Mn", "oxidation": 3},
    "Cu-corrole": {"metal": "Cu", "oxidation": 3},
    "Ga-corrole": {"metal": "Ga", "oxidation": 3},
    "Au-corrole": {"metal": "Au", "oxidation": 3},
}


def generate_porphyrin(
    metalloporphyrin: str = "Fe-porphyrin",
    meso_substituents: Optional[List[str]] = None,
    axial_ligands: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Generate porphyrin or metalloporphyrin structure.
    
    Args:
        metalloporphyrin: Porphyrin type from database
        meso_substituents: Substituents at meso positions (4 values)
        axial_ligands: Axial ligands on metal
    
    Returns:
        Porphyrin structure
    """
    if metalloporphyrin not in PORPHYRIN_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_NAME", "message": f"Unknown porphyrin '{metalloporphyrin}'",
                      "available": list(PORPHYRIN_DATABASE.keys())}
        }
    
    info = PORPHYRIN_DATABASE[metalloporphyrin]
    metal = info.get("metal")
    
    atoms = []
    
    # Core porphyrin ring parameters
    r_metal_n = 2.0  # Metal-N distance
    r_n_c_alpha = 1.38  # N-C(alpha) bond
    r_meso = 4.2  # Meso carbon distance from center
    
    # Central metal (if any)
    if metal:
        atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    else:
        # Two inner H for free-base
        atoms.append({"element": "H", "cartesian": [0.8, 0, 0]})
        atoms.append({"element": "H", "cartesian": [-0.8, 0, 0]})
    
    # Four pyrrole rings
    for i in range(4):
        angle = np.pi / 4 + i * np.pi / 2
        
        # Nitrogen
        nx = r_metal_n * np.cos(angle)
        ny = r_metal_n * np.sin(angle)
        atoms.append({"element": "N", "cartesian": [nx, ny, 0]})
        
        # Alpha carbons (Cα)
        for da in [-0.4, 0.4]:
            angle_c = angle + da
            cx = (r_metal_n + 1.3) * np.cos(angle_c)
            cy = (r_metal_n + 1.3) * np.sin(angle_c)
            atoms.append({"element": "C", "cartesian": [cx, cy, 0]})
        
        # Beta carbons (Cβ)
        for da in [-0.25, 0.25]:
            angle_c = angle + da
            cx = (r_metal_n + 2.4) * np.cos(angle_c)
            cy = (r_metal_n + 2.4) * np.sin(angle_c)
            atoms.append({"element": "C", "cartesian": [cx, cy, 0]})
    
    # Meso carbons (bridging)
    for i in range(4):
        angle = i * np.pi / 2
        mx = r_meso * np.cos(angle)
        my = r_meso * np.sin(angle)
        atoms.append({"element": "C", "cartesian": [mx, my, 0]})
        
        # Meso H or substituent
        if meso_substituents and i < len(meso_substituents):
            sub = meso_substituents[i]
            if sub == "Ph":
                # Simplified phenyl as single C
                atoms.append({"element": "C", "cartesian": [mx * 1.3, my * 1.3, 0]})
            elif sub == "H":
                atoms.append({"element": "H", "cartesian": [mx * 1.15, my * 1.15, 0]})
        else:
            atoms.append({"element": "H", "cartesian": [mx * 1.15, my * 1.15, 0]})
    
    # Beta hydrogens
    for i in range(4):
        angle = np.pi / 4 + i * np.pi / 2
        for da in [-0.25, 0.25]:
            angle_h = angle + da
            hx = (r_metal_n + 3.3) * np.cos(angle_h)
            hy = (r_metal_n + 3.3) * np.sin(angle_h)
            atoms.append({"element": "H", "cartesian": [hx, hy, 0]})
    
    # Axial ligands
    if axial_ligands and metal:
        z_axial = 2.0
        for i, lig in enumerate(axial_ligands):
            z = z_axial if i == 0 else -z_axial
            if lig == "O2":
                atoms.append({"element": "O", "cartesian": [0, 0, z]})
                atoms.append({"element": "O", "cartesian": [0, 0, z + 1.2]})
            elif lig == "CO":
                atoms.append({"element": "C", "cartesian": [0, 0, z]})
                atoms.append({"element": "O", "cartesian": [0, 0, z + 1.1]})
            elif lig in ["Cl", "F", "Br", "I"]:
                atoms.append({"element": lig, "cartesian": [0, 0, z]})
            elif lig == "Im":
                atoms.append({"element": "N", "cartesian": [0, 0, z]})
            elif lig == "pyridine":
                atoms.append({"element": "N", "cartesian": [0, 0, z]})
    
    return {
        "success": True,
        "name": metalloporphyrin,
        "metal": metal,
        "spin": info.get("spin", 0),
        "meso_substituents": meso_substituents or ["H", "H", "H", "H"],
        "axial_ligands": axial_ligands,
        "is_special": info.get("special", False),
        "description": info.get("description", f"{metalloporphyrin} complex"),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_phthalocyanine(
    name: str = "Cu-Pc",
    peripheral_substituents: Optional[str] = None,
    axial_ligand: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate phthalocyanine structure.
    
    Args:
        name: Phthalocyanine from database
        peripheral_substituents: Substituent type (t-Bu, OC8, F16, etc.)
        axial_ligand: Axial ligand
    
    Returns:
        Phthalocyanine structure
    """
    if name not in PHTHALOCYANINE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_NAME", "message": f"Unknown phthalocyanine '{name}'",
                      "available": list(PHTHALOCYANINE_DATABASE.keys())}
        }
    
    info = PHTHALOCYANINE_DATABASE[name]
    metal = info.get("metal")
    
    atoms = []
    
    # Similar to porphyrin but with benzene-fused pyrroles
    r_metal_n = 1.95
    
    if metal:
        atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    else:
        atoms.append({"element": "H", "cartesian": [0.8, 0, 0]})
        atoms.append({"element": "H", "cartesian": [-0.8, 0, 0]})
    
    # Four isoindole units
    for i in range(4):
        angle = np.pi / 4 + i * np.pi / 2
        
        # Pyrrole N
        nx = r_metal_n * np.cos(angle)
        ny = r_metal_n * np.sin(angle)
        atoms.append({"element": "N", "cartesian": [nx, ny, 0]})
        
        # Aza-N (bridging)
        angle_aza = i * np.pi / 2
        ax = 2.8 * np.cos(angle_aza)
        ay = 2.8 * np.sin(angle_aza)
        atoms.append({"element": "N", "cartesian": [ax, ay, 0]})
        
        # Benzene ring (simplified - 6C)
        bz_center_x = 4.5 * np.cos(angle)
        bz_center_y = 4.5 * np.sin(angle)
        for j in range(6):
            bz_angle = j * np.pi / 3
            bx = bz_center_x + 1.4 * np.cos(bz_angle)
            by = bz_center_y + 1.4 * np.sin(bz_angle)
            atoms.append({"element": "C", "cartesian": [bx, by, 0]})
        
        # Benzene H
        for j in range(4):  # Only outer H (4 per benzene)
            bz_angle = np.pi / 6 + j * np.pi / 3
            if j != 1 and j != 2:  # Skip fused positions
                hx = bz_center_x + 2.3 * np.cos(bz_angle)
                hy = bz_center_y + 2.3 * np.sin(bz_angle)
                atoms.append({"element": "H", "cartesian": [hx, hy, 0]})
    
    # Axial ligand
    if axial_ligand and metal:
        atoms.append({"element": axial_ligand, "cartesian": [0, 0, 2.0]})
    
    return {
        "success": True,
        "name": name,
        "metal": metal,
        "peripheral_substituents": peripheral_substituents,
        "axial_ligand": axial_ligand,
        "color": info.get("color", "blue"),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_corrole(
    name: str = "Fe-corrole",
    meso_substituents: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Generate corrole (contracted porphyrin with 3 meso positions).
    
    Args:
        name: Corrole from database
        meso_substituents: 3 meso substituents
    
    Returns:
        Corrole structure
    """
    if name not in CORROLE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_NAME", "message": f"Unknown corrole '{name}'",
                      "available": list(CORROLE_DATABASE.keys())}
        }
    
    info = CORROLE_DATABASE[name]
    metal = info.get("metal")
    
    atoms = []
    
    r_metal_n = 1.9  # Slightly smaller cavity than porphyrin
    
    if metal:
        atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    else:
        # Three inner H for free-base corrole
        for i in range(3):
            angle = 2 * np.pi * i / 3
            atoms.append({"element": "H", "cartesian": [0.6 * np.cos(angle), 0.6 * np.sin(angle), 0]})
    
    # Four pyrrole N (corrole has direct pyrrole-pyrrole link)
    n_angles = [0.4, 1.2, 2.0, 2.8]  # Asymmetric distribution
    for i, angle in enumerate(n_angles):
        nx = r_metal_n * np.cos(angle)
        ny = r_metal_n * np.sin(angle)
        atoms.append({"element": "N", "cartesian": [nx, ny, 0]})
    
    # Three meso positions (not four like porphyrin)
    for i in range(3):
        angle = np.pi / 3 + i * 2 * np.pi / 3
        mx = 3.8 * np.cos(angle)
        my = 3.8 * np.sin(angle)
        atoms.append({"element": "C", "cartesian": [mx, my, 0]})
    
    return {
        "success": True,
        "name": name,
        "metal": metal,
        "oxidation_state": info.get("oxidation", 3),
        "meso_substituents": meso_substituents or ["H", "H", "H"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_subphthalocyanine(
    metal: str = "B",
    axial_ligand: str = "Cl"
) -> Dict[str, Any]:
    """
    Generate subphthalocyanine (bowl-shaped with 3 isoindole units).
    
    Args:
        metal: Central atom (usually B)
        axial_ligand: Axial ligand
    
    Returns:
        Subphthalocyanine structure
    """
    atoms = []
    
    # Boron center (slightly out of plane due to bowl shape)
    atoms.append({"element": metal, "cartesian": [0, 0, 0.5]})
    
    # Three isoindole units arranged in bowl
    for i in range(3):
        angle = i * 2 * np.pi / 3
        
        # Pyrrole N
        r_n = 1.6
        z_n = 0.2  # Bowl curvature
        atoms.append({"element": "N", "cartesian": [r_n * np.cos(angle), r_n * np.sin(angle), z_n]})
        
        # Aza-N (bridging)
        angle_aza = angle + np.pi / 3
        atoms.append({"element": "N", "cartesian": [2.3 * np.cos(angle_aza), 2.3 * np.sin(angle_aza), 0]})
        
        # Benzene (simplified)
        bz_r = 4.0
        atoms.append({"element": "C", "cartesian": [bz_r * np.cos(angle), bz_r * np.sin(angle), -0.3]})
    
    # Axial ligand
    atoms.append({"element": axial_ligand, "cartesian": [0, 0, 2.2]})
    
    return {
        "success": True,
        "name": f"SubPc-{metal}-{axial_ligand}",
        "metal": metal,
        "axial_ligand": axial_ligand,
        "n_isoindole_units": 3,
        "bowl_shaped": True,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
