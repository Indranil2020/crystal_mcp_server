"""
quantum/dots.py - Quantum Dots

Comprehensive quantum dot generation per structure_catalogue.md Category 8:
(i) Quantum dots (CdSe, InP, CsPbBr3, Si, Ge, C, MoS2) with ligand passivation
(ii) Core-shell dots (CdSe/CdS, InP/ZnSe, CsPbBr3/CsPbCl3) with gradient alloy shells
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Quantum dot database
QUANTUM_DOT_DATABASE = {
    # II-VI semiconductors
    "CdSe": {
        "structure": "wurtzite", "a": 4.30, "c": 7.01, "bandgap_bulk_eV": 1.74,
        "effective_mass_e": 0.13, "effective_mass_h": 0.45, "dielectric": 9.7
    },
    "CdS": {
        "structure": "wurtzite", "a": 4.14, "c": 6.72, "bandgap_bulk_eV": 2.42,
        "effective_mass_e": 0.21, "effective_mass_h": 0.80, "dielectric": 8.9
    },
    "CdTe": {
        "structure": "zincblende", "a": 6.48, "bandgap_bulk_eV": 1.50,
        "effective_mass_e": 0.096, "effective_mass_h": 0.35, "dielectric": 10.2
    },
    "ZnSe": {
        "structure": "zincblende", "a": 5.67, "bandgap_bulk_eV": 2.70,
        "effective_mass_e": 0.17, "effective_mass_h": 0.75, "dielectric": 8.1
    },
    "ZnS": {
        "structure": "zincblende", "a": 5.41, "bandgap_bulk_eV": 3.54,
        "effective_mass_e": 0.34, "effective_mass_h": 0.58, "dielectric": 8.3
    },
    
    # III-V semiconductors
    "InP": {
        "structure": "zincblende", "a": 5.87, "bandgap_bulk_eV": 1.35,
        "effective_mass_e": 0.077, "effective_mass_h": 0.60, "dielectric": 12.5
    },
    "InAs": {
        "structure": "zincblende", "a": 6.06, "bandgap_bulk_eV": 0.35,
        "effective_mass_e": 0.026, "effective_mass_h": 0.40, "dielectric": 15.1
    },
    "GaAs": {
        "structure": "zincblende", "a": 5.65, "bandgap_bulk_eV": 1.42,
        "effective_mass_e": 0.067, "effective_mass_h": 0.47, "dielectric": 13.1
    },
    
    # Perovskites
    "CsPbBr3": {
        "structure": "perovskite", "a": 5.87, "bandgap_bulk_eV": 2.30,
        "effective_mass_e": 0.22, "effective_mass_h": 0.24, "dielectric": 25.0
    },
    "CsPbI3": {
        "structure": "perovskite", "a": 6.29, "bandgap_bulk_eV": 1.73,
        "effective_mass_e": 0.15, "effective_mass_h": 0.17, "dielectric": 30.0
    },
    "CsPbCl3": {
        "structure": "perovskite", "a": 5.60, "bandgap_bulk_eV": 2.94,
        "effective_mass_e": 0.20, "effective_mass_h": 0.25, "dielectric": 20.0
    },
    
    # Group IV
    "Si": {
        "structure": "diamond", "a": 5.43, "bandgap_bulk_eV": 1.12,
        "effective_mass_e": 1.08, "effective_mass_h": 0.56, "dielectric": 11.7
    },
    "Ge": {
        "structure": "diamond", "a": 5.66, "bandgap_bulk_eV": 0.66,
        "effective_mass_e": 0.55, "effective_mass_h": 0.37, "dielectric": 16.0
    },
    
    # Carbon
    "C_diamond": {
        "structure": "diamond", "a": 3.57, "bandgap_bulk_eV": 5.47,
        "effective_mass_e": 0.28, "effective_mass_h": 0.48, "dielectric": 5.7
    },
    
    # 2D TMD QDs
    "MoS2_QD": {
        "structure": "2D", "a": 3.16, "c": 6.15, "bandgap_bulk_eV": 1.90,
        "effective_mass_e": 0.48, "effective_mass_h": 0.60, "dielectric": 4.0
    },
}


# Core-shell combinations
CORESHELL_DATABASE = {
    "CdSe/CdS": {"core": "CdSe", "shell": "CdS", "type": "Type_I", "offset_eV": 0.25},
    "CdSe/ZnS": {"core": "CdSe", "shell": "ZnS", "type": "Type_I", "offset_eV": 1.0},
    "CdSe/ZnSe": {"core": "CdSe", "shell": "ZnSe", "type": "Type_I", "offset_eV": 0.4},
    "InP/ZnSe": {"core": "InP", "shell": "ZnSe", "type": "Type_I", "offset_eV": 0.6},
    "InP/ZnS": {"core": "InP", "shell": "ZnS", "type": "Type_I", "offset_eV": 1.0},
    "CsPbBr3/CsPbCl3": {"core": "CsPbBr3", "shell": "CsPbCl3", "type": "Type_I", "offset_eV": 0.3},
    "CdTe/CdSe": {"core": "CdTe", "shell": "CdSe", "type": "Type_II", "offset_eV": 0.4},
    "ZnSe/ZnS": {"core": "ZnSe", "shell": "ZnS", "type": "Type_I", "offset_eV": 0.6},
}


# Ligand database
LIGAND_DATABASE = {
    "oleic_acid": {"formula": "C18H34O2", "length_nm": 2.0, "binding": "carboxylate"},
    "oleylamine": {"formula": "C18H37N", "length_nm": 2.0, "binding": "amine"},
    "TOP": {"formula": "C24H51P", "length_nm": 1.2, "binding": "phosphine"},
    "TOPO": {"formula": "C24H51PO", "length_nm": 1.2, "binding": "phosphine_oxide"},
    "thiol_C12": {"formula": "C12H26S", "length_nm": 1.5, "binding": "thiol"},
    "MPA": {"formula": "C3H6O2S", "length_nm": 0.5, "binding": "thiol"},
    "citrate": {"formula": "C6H5O7", "length_nm": 0.6, "binding": "carboxylate"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def calculate_bandgap_qd(material: str, diameter_nm: float) -> float:
    """Calculate size-dependent bandgap using Brus equation."""
    info = QUANTUM_DOT_DATABASE.get(material)
    if not info:
        return 0
    
    Eg_bulk = info["bandgap_bulk_eV"]
    me = info["effective_mass_e"]
    mh = info["effective_mass_h"]
    eps = info["dielectric"]
    
    # Brus equation
    h = 4.136e-15  # eV*s
    m0 = 9.109e-31  # kg
    e = 1.602e-19  # C
    eps0 = 8.854e-12  # F/m
    
    r = diameter_nm / 2 * 1e-9  # meters
    
    # Kinetic term
    kinetic = (h**2 / (8 * r**2)) * (1/(me * m0) + 1/(mh * m0)) * 1e20 * 6.242e18  # Convert to eV
    
    # Coulomb term
    coulomb = -1.786 * e / (4 * np.pi * eps * eps0 * r) / e  # eV
    
    Eg_qd = Eg_bulk + kinetic + coulomb
    
    return max(Eg_bulk, Eg_qd)  # Can't be less than bulk


def generate_quantum_dot(
    material: str = "CdSe",
    diameter_nm: float = 4.0,
    shape: str = "sphere",
    ligand: Optional[str] = "oleic_acid",
    faceted: bool = False
) -> Dict[str, Any]:
    """
    Generate quantum dot structure.
    
    Args:
        material: Material from database
        diameter_nm: Dot diameter in nm
        shape: 'sphere', 'cube', 'rod', 'tetrapod', 'platelet'
        ligand: Ligand for surface passivation
        faceted: Use faceted structure
    
    Returns:
        Quantum dot structure
    """
    if material not in QUANTUM_DOT_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material '{material}'",
                      "available": list(QUANTUM_DOT_DATABASE.keys())}
        }
    
    info = QUANTUM_DOT_DATABASE[material]
    a = info["a"]
    
    diameter = diameter_nm * 10  # Angstrom
    radius = diameter / 2
    
    # Box size
    box_size = diameter + 20  # Extra space for ligands
    lattice = Lattice.cubic(box_size)
    
    species = []
    coords = []
    
    # Extract elements
    if "Cd" in material:
        elem1, elem2 = "Cd", "Se" if "Se" in material else ("S" if "S" in material else "Te")
    elif "Zn" in material:
        elem1, elem2 = "Zn", "Se" if "Se" in material else "S"
    elif "InP" in material:
        elem1, elem2 = "In", "P"
    elif "InAs" in material:
        elem1, elem2 = "In", "As"
    elif "GaAs" in material:
        elem1, elem2 = "Ga", "As"
    elif "CsPb" in material:
        elem1, elem2, elem3 = "Cs", "Pb", "Br" if "Br" in material else ("I" if "I" in material else "Cl")
    elif material in ["Si", "Ge"]:
        elem1, elem2 = material, material
    else:
        elem1, elem2 = "X", "X"
    
    # Generate atom positions based on shape
    if shape == "sphere":
        n_cells = int(diameter / a) + 1
        for i in range(-n_cells, n_cells + 1):
            for j in range(-n_cells, n_cells + 1):
                for k in range(-n_cells, n_cells + 1):
                    # Zincblende basis
                    positions = [
                        np.array([i, j, k]) * a,
                        np.array([i + 0.5, j + 0.5, k]) * a,
                        np.array([i + 0.5, j, k + 0.5]) * a,
                        np.array([i, j + 0.5, k + 0.5]) * a,
                        np.array([i + 0.25, j + 0.25, k + 0.25]) * a,
                        np.array([i + 0.75, j + 0.75, k + 0.25]) * a,
                        np.array([i + 0.75, j + 0.25, k + 0.75]) * a,
                        np.array([i + 0.25, j + 0.75, k + 0.75]) * a,
                    ]
                    elements = [elem1, elem1, elem1, elem1, elem2, elem2, elem2, elem2]
                    
                    for pos, elem in zip(positions, elements):
                        r = np.linalg.norm(pos)
                        if r <= radius:
                            frac_x = (pos[0] + box_size/2) / box_size
                            frac_y = (pos[1] + box_size/2) / box_size
                            frac_z = (pos[2] + box_size/2) / box_size
                            species.append(elem)
                            coords.append([frac_x, frac_y, frac_z])
    
    elif shape == "cube":
        n_cells = int(diameter / a) + 1
        for i in range(n_cells):
            for j in range(n_cells):
                for k in range(n_cells):
                    x = (i + 0.25) / n_cells * 0.5 + 0.25
                    y = (j + 0.25) / n_cells * 0.5 + 0.25
                    z = (k + 0.25) / n_cells * 0.5 + 0.25
                    species.extend([elem1, elem2])
                    coords.extend([[x, y, z], [x + 0.02, y + 0.02, z + 0.02]])
    
    elif shape == "rod":
        aspect_ratio = 3
        length = diameter * aspect_ratio
        n_z = int(length / a)
        n_xy = int(diameter / a)
        
        for i in range(-n_xy, n_xy + 1):
            for j in range(-n_xy, n_xy + 1):
                for k in range(n_z):
                    x = i * a
                    y = j * a
                    z = k * a - length/2
                    
                    r_xy = np.sqrt(x**2 + y**2)
                    if r_xy <= radius:
                        fx = (x + box_size/2) / box_size
                        fy = (y + box_size/2) / box_size
                        fz = (z + box_size/2) / box_size
                        species.extend([elem1, elem2])
                        coords.extend([[fx, fy, fz], [fx + 0.01, fy + 0.01, fz]])
    
    else:  # Generic
        n = int((diameter / a) ** 3 * 4)
        np.random.seed(42)
        for _ in range(n):
            r = np.random.uniform(0, radius)
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(0, 2 * np.pi)
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            fx = (x + box_size/2) / box_size
            fy = (y + box_size/2) / box_size
            fz = (z + box_size/2) / box_size
            species.append(elem1 if _ % 2 == 0 else elem2)
            coords.append([fx, fy, fz])
    
    # Add ligands
    n_ligands = 0
    if ligand and ligand in LIGAND_DATABASE:
        lig_info = LIGAND_DATABASE[ligand]
        lig_length = lig_info["length_nm"] * 10  # Angstrom
        
        # Find surface atoms and add ligands
        surface_atoms = []
        for i, (sp, c) in enumerate(zip(species, coords)):
            pos = np.array(c) * box_size - box_size/2
            r = np.linalg.norm(pos)
            if r > radius - a:  # Surface atom
                surface_atoms.append((i, pos))
        
        # Add ligand atoms at surface
        for i, pos in surface_atoms[::4]:  # Every 4th surface atom
            direction = pos / (np.linalg.norm(pos) + 0.001)
            lig_pos = pos + direction * (lig_length / 2)
            fx = (lig_pos[0] + box_size/2) / box_size
            fy = (lig_pos[1] + box_size/2) / box_size
            fz = (lig_pos[2] + box_size/2) / box_size
            
            if 0 < fx < 1 and 0 < fy < 1 and 0 < fz < 1:
                species.append("C")  # Simplified ligand
                coords.append([fx, fy, fz])
                n_ligands += 1
    
    structure = Structure(lattice, species, coords)
    
    # Calculate bandgap
    bandgap = calculate_bandgap_qd(material, diameter_nm)
    
    return {
        "success": True,
        "material": material,
        "diameter_nm": diameter_nm,
        "shape": shape,
        "bandgap_eV": round(bandgap, 3),
        "bandgap_bulk_eV": info["bandgap_bulk_eV"],
        "blueshift_eV": round(bandgap - info["bandgap_bulk_eV"], 3),
        "n_core_atoms": len(structure) - n_ligands,
        "n_ligands": n_ligands,
        "ligand": ligand,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_coreshell_qd(
    system: str = "CdSe/CdS",
    core_diameter_nm: float = 3.0,
    shell_thickness_nm: float = 1.0,
    gradient_shell: bool = False
) -> Dict[str, Any]:
    """
    Generate core-shell quantum dot.
    
    Args:
        system: Core/shell material combination
        core_diameter_nm: Core diameter
        shell_thickness_nm: Shell thickness
        gradient_shell: Use gradient alloy shell
    
    Returns:
        Core-shell QD structure
    """
    if system not in CORESHELL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", "message": f"Unknown system '{system}'",
                      "available": list(CORESHELL_DATABASE.keys())}
        }
    
    info = CORESHELL_DATABASE[system]
    core_mat = info["core"]
    shell_mat = info["shell"]
    
    total_diameter = core_diameter_nm + 2 * shell_thickness_nm
    core_radius = core_diameter_nm / 2 * 10  # Angstrom
    total_radius = total_diameter / 2 * 10
    
    box_size = total_diameter * 10 + 20
    lattice = Lattice.cubic(box_size)
    
    species = []
    coords = []
    
    # Get elements
    core_info = QUANTUM_DOT_DATABASE[core_mat]
    shell_info = QUANTUM_DOT_DATABASE[shell_mat]
    a_core = core_info["a"]
    a_shell = shell_info["a"]
    
    # Core atoms
    n_cells = int(core_diameter_nm * 10 / a_core) + 1
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            for k in range(-n_cells, n_cells + 1):
                x = i * a_core
                y = j * a_core
                z = k * a_core
                r = np.sqrt(x**2 + y**2 + z**2)
                
                if r <= core_radius:
                    fx = (x + box_size/2) / box_size
                    fy = (y + box_size/2) / box_size
                    fz = (z + box_size/2) / box_size
                    
                    elem1 = core_mat[:2] if len(core_mat) > 2 else core_mat
                    elem2 = core_mat[-2:] if len(core_mat) > 3 else core_mat[-1]
                    species.extend([elem1, elem2])
                    coords.extend([[fx, fy, fz], [fx + 0.01, fy + 0.01, fz]])
    
    n_core = len(species)
    
    # Shell atoms
    n_cells = int(total_diameter * 10 / a_shell) + 1
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            for k in range(-n_cells, n_cells + 1):
                x = i * a_shell
                y = j * a_shell
                z = k * a_shell
                r = np.sqrt(x**2 + y**2 + z**2)
                
                if core_radius < r <= total_radius:
                    fx = (x + box_size/2) / box_size
                    fy = (y + box_size/2) / box_size
                    fz = (z + box_size/2) / box_size
                    
                    elem1 = shell_mat[:2] if len(shell_mat) > 2 else shell_mat
                    elem2 = shell_mat[-2:] if len(shell_mat) > 3 else shell_mat[-1]
                    species.extend([elem1, elem2])
                    coords.extend([[fx, fy, fz], [fx + 0.01, fy + 0.01, fz]])
    
    n_shell = len(species) - n_core
    
    structure = Structure(lattice, species, coords)
    
    # Calculate bandgaps
    core_bandgap = calculate_bandgap_qd(core_mat, core_diameter_nm)
    
    return {
        "success": True,
        "system": system,
        "core_material": core_mat,
        "shell_material": shell_mat,
        "core_diameter_nm": core_diameter_nm,
        "shell_thickness_nm": shell_thickness_nm,
        "total_diameter_nm": total_diameter,
        "band_alignment": info["type"],
        "core_bandgap_eV": round(core_bandgap, 3),
        "n_core_atoms": n_core,
        "n_shell_atoms": n_shell,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
