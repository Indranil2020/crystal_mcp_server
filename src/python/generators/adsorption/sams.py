"""
adsorption/sams.py - Self-Assembled Monolayers

Generates SAM structures:
- Alkanethiols on Au(111)
- Silanes on SiO2
- Phosphonates on TiO2
- Various chain lengths and functional groups
"""

from typing import Dict, Any, List, Optional
import numpy as np


# SAM database
SAM_DATABASE = {
    # Thiol SAMs on Au
    "C6-thiol/Au": {
        "anchor": "S", "substrate": "Au", "chain_length": 6,
        "tilt_angle": 30, "spacing": 4.99, "description": "Hexanethiol on Au(111)"
    },
    "C12-thiol/Au": {
        "anchor": "S", "substrate": "Au", "chain_length": 12,
        "tilt_angle": 30, "spacing": 4.99, "description": "Dodecanethiol on Au(111)"
    },
    "C18-thiol/Au": {
        "anchor": "S", "substrate": "Au", "chain_length": 18,
        "tilt_angle": 30, "spacing": 4.99, "description": "Octadecanethiol on Au(111)"
    },
    "MUA/Au": {
        "anchor": "S", "substrate": "Au", "chain_length": 11,
        "terminal": "COOH", "tilt_angle": 30, "spacing": 4.99,
        "description": "Mercaptoundecanoic acid"
    },
    "MUO/Au": {
        "anchor": "S", "substrate": "Au", "chain_length": 11,
        "terminal": "OH", "tilt_angle": 30, "spacing": 4.99,
        "description": "Mercaptoundecanol"
    },
    
    # Silane SAMs on SiO2
    "OTS/SiO2": {
        "anchor": "Si-O", "substrate": "SiO2", "chain_length": 18,
        "tilt_angle": 10, "spacing": 4.5, "description": "Octadecyltrichlorosilane"
    },
    "APTES/SiO2": {
        "anchor": "Si-O", "substrate": "SiO2", "chain_length": 3,
        "terminal": "NH2", "tilt_angle": 15, "spacing": 4.5,
        "description": "Aminopropyltriethoxysilane"
    },
    "MPTMS/SiO2": {
        "anchor": "Si-O", "substrate": "SiO2", "chain_length": 3,
        "terminal": "SH", "spacing": 4.5, "description": "Mercaptopropyltrimethoxysilane"
    },
    
    # Phosphonate SAMs on TiO2
    "ODPA/TiO2": {
        "anchor": "P-O", "substrate": "TiO2", "chain_length": 18,
        "tilt_angle": 25, "spacing": 5.0, "description": "Octadecylphosphonic acid"
    },
    "PPA/TiO2": {
        "anchor": "P-O", "substrate": "TiO2", "chain_length": 3,
        "terminal": "Ph", "spacing": 5.0, "description": "Phenylphosphonic acid"
    },
}


def generate_sam(
    sam_type: str = "C12-thiol/Au",
    size: List[int] = [4, 4],
    defect_density: float = 0.0
) -> Dict[str, Any]:
    """
    Generate self-assembled monolayer.
    
    Args:
        sam_type: SAM type from database
        size: Surface size [nx, ny]
        defect_density: Fraction of missing molecules
    
    Returns:
        SAM structure
    """
    if sam_type not in SAM_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_SAM", "message": f"Unknown SAM '{sam_type}'",
                      "available": list(SAM_DATABASE.keys())}
        }
    
    info = SAM_DATABASE[sam_type]
    chain_length = info["chain_length"]
    tilt_angle = np.radians(info.get("tilt_angle", 30))
    spacing = info["spacing"]
    anchor = info["anchor"]
    terminal = info.get("terminal", "CH3")
    
    nx, ny = size
    atoms = []
    n_molecules = 0
    
    np.random.seed(42)
    
    # C-C bond length
    c_c = 1.54
    c_s = 1.82 if anchor == "S" else 1.90
    
    # Generate substrate (simplified)
    substrate = info["substrate"]
    for i in range(nx + 2):
        for j in range(ny + 2):
            x = i * spacing * 0.87
            y = j * spacing
            
            if substrate == "Au":
                atoms.append({"element": "Au", "cartesian": [x, y, 0], "layer": "substrate"})
            elif substrate == "SiO2":
                atoms.append({"element": "Si", "cartesian": [x, y, 0], "layer": "substrate"})
                atoms.append({"element": "O", "cartesian": [x, y, 1.6], "layer": "substrate"})
            else:  # TiO2
                atoms.append({"element": "Ti", "cartesian": [x, y, 0], "layer": "substrate"})
                atoms.append({"element": "O", "cartesian": [x, y, 2.0], "layer": "substrate"})
    
    # Generate SAM molecules
    z_anchor = {"Au": 2.0, "SiO2": 3.0, "TiO2": 3.5}[substrate]
    
    for i in range(nx):
        for j in range(ny):
            # Skip defect positions
            if np.random.random() < defect_density:
                continue
            
            x_base = i * spacing * 0.87 + spacing * 0.87 / 2
            y_base = j * spacing + (i % 2) * spacing / 2
            
            # Anchor atom (S for thiols)
            if anchor == "S":
                atoms.append({"element": "S", "cartesian": [x_base, y_base, z_anchor], "molecule": n_molecules})
            elif anchor == "Si-O":
                atoms.append({"element": "Si", "cartesian": [x_base, y_base, z_anchor], "molecule": n_molecules})
            else:  # P-O
                atoms.append({"element": "P", "cartesian": [x_base, y_base, z_anchor], "molecule": n_molecules})
            
            # Alkyl chain
            dx = c_c * np.sin(tilt_angle)
            dz = c_c * np.cos(tilt_angle)
            
            for c in range(chain_length):
                z = z_anchor + c_s + c * dz
                x = x_base + c * dx * (1 if c % 2 == 0 else -1) * 0.5
                atoms.append({"element": "C", "cartesian": [x, y_base, z], "molecule": n_molecules})
            
            # Terminal group
            z_terminal = z_anchor + c_s + chain_length * dz
            x_terminal = x_base + chain_length * dx * 0.3
            
            if terminal == "CH3":
                atoms.append({"element": "C", "cartesian": [x_terminal, y_base, z_terminal], "molecule": n_molecules})
            elif terminal == "COOH":
                atoms.append({"element": "C", "cartesian": [x_terminal, y_base, z_terminal], "molecule": n_molecules})
                atoms.append({"element": "O", "cartesian": [x_terminal + 0.8, y_base, z_terminal + 0.6], "molecule": n_molecules})
                atoms.append({"element": "O", "cartesian": [x_terminal - 0.5, y_base, z_terminal + 1.0], "molecule": n_molecules})
                atoms.append({"element": "H", "cartesian": [x_terminal - 1.0, y_base, z_terminal + 1.7], "molecule": n_molecules})
            elif terminal == "OH":
                atoms.append({"element": "O", "cartesian": [x_terminal, y_base, z_terminal], "molecule": n_molecules})
                atoms.append({"element": "H", "cartesian": [x_terminal + 0.5, y_base, z_terminal + 0.5], "molecule": n_molecules})
            elif terminal == "NH2":
                atoms.append({"element": "N", "cartesian": [x_terminal, y_base, z_terminal], "molecule": n_molecules})
                atoms.append({"element": "H", "cartesian": [x_terminal + 0.5, y_base + 0.3, z_terminal + 0.5], "molecule": n_molecules})
                atoms.append({"element": "H", "cartesian": [x_terminal - 0.5, y_base - 0.3, z_terminal + 0.5], "molecule": n_molecules})
            elif terminal == "SH":
                atoms.append({"element": "S", "cartesian": [x_terminal, y_base, z_terminal], "molecule": n_molecules})
                atoms.append({"element": "H", "cartesian": [x_terminal + 0.5, y_base, z_terminal + 0.8], "molecule": n_molecules})
            
            n_molecules += 1
    
    # Calculate properties
    monolayer_thickness = c_s + chain_length * c_c * np.cos(tilt_angle) + 2
    
    return {
        "success": True,
        "sam_type": sam_type,
        "description": info["description"],
        "substrate": substrate,
        "anchor_group": anchor,
        "chain_length": chain_length,
        "terminal_group": terminal,
        "tilt_angle_deg": info.get("tilt_angle", 30),
        "n_molecules": n_molecules,
        "defect_density": defect_density,
        "monolayer_thickness_angstrom": round(monolayer_thickness, 1),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_mixed_sam(
    sam1: str = "C12-thiol/Au",
    sam2: str = "MUA/Au",
    ratio: float = 0.5,
    mixing_pattern: str = "random",
    size: List[int] = [6, 6]
) -> Dict[str, Any]:
    """
    Generate mixed SAM with two components.
    
    Args:
        sam1: First SAM type
        sam2: Second SAM type
        ratio: Fraction of sam1 (0-1)
        mixing_pattern: 'random', 'stripe', 'island'
        size: Surface size
    
    Returns:
        Mixed SAM structure
    """
    if sam1 not in SAM_DATABASE or sam2 not in SAM_DATABASE:
        return {"success": False, "error": {"code": "INVALID_SAM", "message": "Unknown SAM type"}}
    
    # Check compatibility (same substrate)
    if SAM_DATABASE[sam1]["substrate"] != SAM_DATABASE[sam2]["substrate"]:
        return {"success": False, "error": {"code": "INCOMPATIBLE", "message": "SAMs must have same substrate"}}
    
    info1 = SAM_DATABASE[sam1]
    info2 = SAM_DATABASE[sam2]
    substrate = info1["substrate"]
    spacing = info1["spacing"]
    
    nx, ny = size
    atoms = []
    n_sam1 = 0
    n_sam2 = 0
    
    np.random.seed(42)
    
    # Substrate
    z_anchor = {"Au": 2.0, "SiO2": 3.0, "TiO2": 3.5}[substrate]
    for i in range(nx + 2):
        for j in range(ny + 2):
            x = i * spacing * 0.87
            y = j * spacing
            atoms.append({"element": substrate.replace("O2", ""), "cartesian": [x, y, 0], "layer": "substrate"})
    
    # Mixed molecules
    for i in range(nx):
        for j in range(ny):
            x_base = i * spacing * 0.87 + spacing * 0.87 / 2
            y_base = j * spacing + (i % 2) * spacing / 2
            
            # Determine which SAM
            if mixing_pattern == "random":
                use_sam1 = np.random.random() < ratio
            elif mixing_pattern == "stripe":
                use_sam1 = i % 2 == 0
            else:  # island
                use_sam1 = (i < nx * ratio) and (j < ny * ratio)
            
            info = info1 if use_sam1 else info2
            chain_length = info["chain_length"]
            tilt = np.radians(info.get("tilt_angle", 30))
            
            atoms.append({"element": "S", "cartesian": [x_base, y_base, z_anchor]})
            
            for c in range(chain_length):
                z = z_anchor + 1.82 + c * 1.54 * np.cos(tilt)
                atoms.append({"element": "C", "cartesian": [x_base, y_base, z]})
            
            if use_sam1:
                n_sam1 += 1
            else:
                n_sam2 += 1
    
    return {
        "success": True,
        "sam1": sam1,
        "sam2": sam2,
        "ratio": ratio,
        "mixing_pattern": mixing_pattern,
        "n_sam1": n_sam1,
        "n_sam2": n_sam2,
        "actual_ratio": n_sam1 / (n_sam1 + n_sam2) if (n_sam1 + n_sam2) > 0 else 0,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
