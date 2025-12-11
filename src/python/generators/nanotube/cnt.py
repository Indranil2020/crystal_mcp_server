"""
nanotube/cnt.py - Carbon Nanotubes

Comprehensive carbon nanotube generation:
- Single-walled CNTs (all chiralities)
- Multi-walled CNTs
- Functionalized CNTs
- Capped and open-ended
- CNT bundles
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


# CNT property database
CNT_DATABASE = {
    # Common CNTs with properties
    "(5,5)": {"diameter": 6.78, "type": "armchair", "metallic": True},
    "(6,6)": {"diameter": 8.14, "type": "armchair", "metallic": True},
    "(10,10)": {"diameter": 13.56, "type": "armchair", "metallic": True},
    "(10,0)": {"diameter": 7.83, "type": "zigzag", "metallic": True},
    "(11,0)": {"diameter": 8.61, "type": "zigzag", "metallic": False},
    "(17,0)": {"diameter": 13.31, "type": "zigzag", "metallic": False},
    "(6,5)": {"diameter": 7.47, "type": "chiral", "metallic": False},
    "(7,5)": {"diameter": 8.17, "type": "chiral", "metallic": False},
    "(8,6)": {"diameter": 9.52, "type": "chiral", "metallic": False},
    "(9,7)": {"diameter": 10.88, "type": "chiral", "metallic": False},
    "(10,5)": {"diameter": 10.36, "type": "chiral", "metallic": True},
}


# Functional groups for CNT modification
FUNCTIONAL_GROUPS = {
    "COOH": {"elements": ["C", "O", "O", "H"], "attachment": "sp3", "effect": "hydrophilic"},
    "OH": {"elements": ["O", "H"], "attachment": "sp3", "effect": "hydrophilic"},
    "NH2": {"elements": ["N", "H", "H"], "attachment": "sp3", "effect": "basic"},
    "F": {"elements": ["F"], "attachment": "sp3", "effect": "fluorinated"},
    "phenyl": {"elements": ["C"]*6 + ["H"]*5, "attachment": "sp3", "effect": "aromatic"},
    "alkyl_C6": {"elements": ["C"]*6 + ["H"]*13, "attachment": "sp3", "effect": "hydrophobic"},
}


def calculate_cnt_properties(n: int, m: int) -> Dict[str, Any]:
    """Calculate CNT properties from chirality indices."""
    a_cc = 1.42  # C-C bond length
    a = a_cc * np.sqrt(3)  # Graphene lattice constant
    
    # Chiral vector magnitude
    Ch = a * np.sqrt(n**2 + m**2 + n*m)
    
    # Diameter
    diameter = Ch / np.pi
    
    # Chiral angle
    theta = np.arctan(np.sqrt(3) * m / (2*n + m))
    
    # Electronic type
    if (n - m) % 3 == 0:
        electronic = "metallic"
    else:
        electronic = "semiconducting"
        # Band gap estimate (empirical)
        bandgap = 0.7 / diameter  # eV
    
    # Chirality type
    if n == m:
        chirality_type = "armchair"
    elif m == 0:
        chirality_type = "zigzag"
    else:
        chirality_type = "chiral"
    
    # Translation vector length
    d_R = np.gcd(2*n + m, 2*m + n)
    T = np.sqrt(3) * Ch / d_R
    
    # Atoms per unit cell
    n_atoms = int(4 * (n**2 + m**2 + n*m) / d_R)
    
    return {
        "diameter_angstrom": round(diameter, 3),
        "chiral_angle_deg": round(np.degrees(theta), 2),
        "electronic_type": electronic,
        "chirality_type": chirality_type,
        "bandgap_eV": round(0.7/diameter, 3) if electronic == "semiconducting" else 0,
        "translation_length_angstrom": round(T, 3),
        "atoms_per_unit_cell": n_atoms,
    }


def generate_cnt(
    n: int = 10,
    m: int = 10,
    length_nm: float = 5.0,
    capped: bool = False,
    open_ended: bool = True
) -> Dict[str, Any]:
    """
    Generate single-walled carbon nanotube.
    
    Args:
        n: First chiral index
        m: Second chiral index
        length_nm: Tube length in nanometers
        capped: Add hemispherical caps
        open_ended: Leave ends open
    
    Returns:
        CNT structure with properties
    """
    if n < 0 or m < 0 or (n == 0 and m == 0):
        return {"success": False, "error": {"code": "INVALID_CHIRALITY", "message": "Invalid (n,m) indices"}}
    
    # Calculate properties
    props = calculate_cnt_properties(n, m)
    
    a_cc = 1.42
    a = a_cc * np.sqrt(3)
    
    diameter = props["diameter_angstrom"]
    radius = diameter / 2
    length = length_nm * 10  # Convert to Angstrom
    
    # Translation vector
    d_R = np.gcd(2*n + m, 2*m + n)
    T_len = props["translation_length_angstrom"]
    
    atoms = []
    
    # Graphene basis vectors
    a1 = np.array([a, 0])
    a2 = np.array([a/2, a * np.sqrt(3)/2])
    
    # Chiral vector
    Ch = n * a1 + m * a2
    Ch_len = np.linalg.norm(Ch)
    
    # Translation vector direction
    T_vec = np.array([(2*m + n), -(2*n + m)]) / d_R
    T_vec = T_vec / np.linalg.norm(T_vec) * T_len
    
    # Number of unit cells along length
    n_units = max(1, int(length / T_len))
    
    # Generate graphene-like positions and roll into tube
    graphene_basis = [
        np.array([0, 0]),
        np.array([a/2, a_cc/2])
    ]
    
    # Sample circumference
    n_circum = max(6, int(Ch_len / a_cc))
    
    for t in range(n_units):
        z_base = t * T_len
        
        for c in range(n_circum):
            phi = 2 * np.pi * c / n_circum
            
            for basis in graphene_basis:
                # Position on tube
                x = radius * np.cos(phi + basis[0] / radius)
                y = radius * np.sin(phi + basis[0] / radius)
                z = z_base + basis[1]
                
                atoms.append({"element": "C", "cartesian": [x, y, z]})
    
    # Add caps if requested
    if capped:
        # C60 hemisphere approximation
        for end_z in [0, length]:
            sign = 1 if end_z == 0 else -1
            for i in range(30):  # Half of C60
                phi = 2 * np.pi * i / 30
                theta = np.pi * (i % 5) / 10
                
                x = radius * 0.8 * np.sin(theta) * np.cos(phi)
                y = radius * 0.8 * np.sin(theta) * np.sin(phi)
                z = end_z + sign * radius * np.cos(theta)
                
                atoms.append({"element": "C", "cartesian": [x, y, z], "cap": True})
    
    return {
        "success": True,
        "chirality": [n, m],
        **props,
        "length_nm": length_nm,
        "capped": capped,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_mwcnt(
    inner_chirality: Tuple[int, int] = (5, 5),
    n_walls: int = 3,
    length_nm: float = 5.0,
    wall_spacing: float = 3.4
) -> Dict[str, Any]:
    """
    Generate multi-walled carbon nanotube.
    
    Args:
        inner_chirality: (n, m) for innermost tube
        n_walls: Number of concentric walls
        length_nm: Tube length
        wall_spacing: Interlayer spacing (3.4 Å typical)
    
    Returns:
        MWCNT structure
    """
    n, m = inner_chirality
    
    atoms = []
    wall_info = []
    
    for wall in range(n_walls):
        # Scale chirality for each wall
        scale = 1 + wall * wall_spacing / calculate_cnt_properties(n, m)["diameter_angstrom"]
        wall_n = max(1, int(n * scale))
        wall_m = max(1, int(m * scale))
        
        # Generate this wall
        wall_result = generate_cnt(wall_n, wall_m, length_nm)
        
        if wall_result["success"]:
            for atom in wall_result["structure"]["atoms"]:
                atom["wall"] = wall
                atoms.append(atom)
            
            wall_info.append({
                "wall": wall,
                "chirality": [wall_n, wall_m],
                "diameter_angstrom": wall_result["diameter_angstrom"],
                "electronic_type": wall_result["electronic_type"]
            })
    
    inner_diameter = wall_info[0]["diameter_angstrom"] if wall_info else 0
    outer_diameter = wall_info[-1]["diameter_angstrom"] if wall_info else 0
    
    return {
        "success": True,
        "type": "MWCNT",
        "n_walls": n_walls,
        "inner_chirality": list(inner_chirality),
        "inner_diameter_angstrom": round(inner_diameter, 2),
        "outer_diameter_angstrom": round(outer_diameter, 2),
        "wall_spacing_angstrom": wall_spacing,
        "walls": wall_info,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_functionalized_cnt(
    n: int = 10,
    m: int = 10,
    length_nm: float = 3.0,
    functional_group: str = "COOH",
    coverage: float = 0.05
) -> Dict[str, Any]:
    """
    Generate functionalized CNT.
    
    Args:
        n, m: Chirality indices
        length_nm: Tube length
        functional_group: Functional group type
        coverage: Fraction of C atoms functionalized
    
    Returns:
        Functionalized CNT
    """
    if functional_group not in FUNCTIONAL_GROUPS:
        return {
            "success": False,
            "error": {"code": "INVALID_GROUP", "message": f"Unknown group '{functional_group}'",
                      "available": list(FUNCTIONAL_GROUPS.keys())}
        }
    
    # Generate base CNT
    cnt = generate_cnt(n, m, length_nm)
    
    if not cnt["success"]:
        return cnt
    
    atoms = list(cnt["structure"]["atoms"])
    n_carbon = len(atoms)
    n_functionalized = max(1, int(n_carbon * coverage))
    
    np.random.seed(42)
    func_indices = np.random.choice(n_carbon, n_functionalized, replace=False)
    
    group_info = FUNCTIONAL_GROUPS[functional_group]
    
    for idx in func_indices:
        base_atom = atoms[idx]
        base_pos = np.array(base_atom["cartesian"])
        
        # Radial direction for attachment
        radial = base_pos[:2]
        radial = radial / (np.linalg.norm(radial) + 1e-10)
        
        # Add functional group atoms
        for i, elem in enumerate(group_info["elements"]):
            offset = 1.5 + i * 0.4
            new_pos = base_pos + np.array([radial[0] * offset, radial[1] * offset, i * 0.3])
            atoms.append({"element": elem, "cartesian": list(new_pos), "functional_group": functional_group})
    
    return {
        "success": True,
        "chirality": [n, m],
        "functional_group": functional_group,
        "effect": group_info["effect"],
        "coverage": coverage,
        "n_functionalized": n_functionalized,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_cnt_bundle(
    chirality: Tuple[int, int] = (10, 10),
    n_tubes: int = 7,
    length_nm: float = 5.0,
    bundle_arrangement: str = "hexagonal"
) -> Dict[str, Any]:
    """
    Generate CNT bundle.
    
    Args:
        chirality: (n, m) for all tubes
        n_tubes: Number of tubes
        length_nm: Tube length
        bundle_arrangement: 'hexagonal' or 'square'
    
    Returns:
        CNT bundle structure
    """
    n, m = chirality
    props = calculate_cnt_properties(n, m)
    diameter = props["diameter_angstrom"]
    
    # Inter-tube spacing (van der Waals gap)
    tube_spacing = diameter + 3.4
    
    atoms = []
    tube_centers = []
    
    if bundle_arrangement == "hexagonal" and n_tubes >= 7:
        # Hexagonal close-packed arrangement
        # Center tube
        tube_centers.append([0, 0])
        # First ring (6 tubes)
        for i in range(6):
            angle = i * np.pi / 3
            x = tube_spacing * np.cos(angle)
            y = tube_spacing * np.sin(angle)
            tube_centers.append([x, y])
    else:
        # Square arrangement
        n_side = int(np.ceil(np.sqrt(n_tubes)))
        for i in range(n_tubes):
            x = (i % n_side) * tube_spacing
            y = (i // n_side) * tube_spacing
            tube_centers.append([x, y])
    
    for i, center in enumerate(tube_centers[:n_tubes]):
        tube = generate_cnt(n, m, length_nm)
        if tube["success"]:
            for atom in tube["structure"]["atoms"]:
                new_pos = [
                    atom["cartesian"][0] + center[0],
                    atom["cartesian"][1] + center[1],
                    atom["cartesian"][2]
                ]
                atoms.append({"element": "C", "cartesian": new_pos, "tube_id": i})
    
    return {
        "success": True,
        "chirality": list(chirality),
        "n_tubes": n_tubes,
        "arrangement": bundle_arrangement,
        "tube_spacing_angstrom": round(tube_spacing, 2),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_doped_cnt(
    n: int = 10,
    m: int = 10,
    length_nm: float = 5.0,
    dopant: str = "N",
    doping_level: float = 0.02
) -> Dict[str, Any]:
    """
    Generate doped CNT.
    
    Args:
        n, m: Chirality
        length_nm: Length
        dopant: Dopant element (N, B, P)
        doping_level: Fraction of substituted atoms
    
    Returns:
        Doped CNT structure
    """
    dopants = {
        "N": {"effect": "n-type", "bond_length_change": -0.02},
        "B": {"effect": "p-type", "bond_length_change": 0.02},
        "P": {"effect": "n-type", "bond_length_change": 0.05},
        "S": {"effect": "n-type", "bond_length_change": 0.10},
    }
    
    if dopant not in dopants:
        return {"success": False, "error": {"code": "INVALID_DOPANT", "message": f"Unknown dopant '{dopant}'"}}
    
    cnt = generate_cnt(n, m, length_nm)
    if not cnt["success"]:
        return cnt
    
    atoms = list(cnt["structure"]["atoms"])
    n_carbon = len(atoms)
    n_doped = max(1, int(n_carbon * doping_level))
    
    np.random.seed(42)
    dope_indices = np.random.choice(n_carbon, n_doped, replace=False)
    
    for idx in dope_indices:
        atoms[idx]["element"] = dopant
        atoms[idx]["doped"] = True
    
    return {
        "success": True,
        "chirality": [n, m],
        "dopant": dopant,
        "doping_level": doping_level,
        "n_dopant_atoms": n_doped,
        "electronic_effect": dopants[dopant]["effect"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
