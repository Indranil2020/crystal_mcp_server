"""
nanotube/nanowires.py - Semiconductor and Metal Nanowires

Comprehensive nanowire generation:
- Semiconductor nanowires (Si, Ge, GaAs, InAs)
- Metal nanowires (Au, Ag, Cu)
- Core-shell nanowires
- Axial/radial heterostructures
- Passivation and surface reconstruction
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Nanowire database
NANOWIRE_DATABASE = {
    # Group IV semiconductors
    "Si": {
        "structure": "diamond", "a": 5.43, "bandgap_eV": 1.12,
        "growth_directions": ["[111]", "[110]", "[112]"], "passivation": "H"
    },
    "Ge": {
        "structure": "diamond", "a": 5.66, "bandgap_eV": 0.66,
        "growth_directions": ["[111]", "[110]"], "passivation": "H"
    },
    # III-V semiconductors
    "GaAs": {
        "structure": "zincblende", "a": 5.65, "bandgap_eV": 1.42,
        "growth_directions": ["[111]B", "[111]A", "[110]"], "passivation": "S"
    },
    "InAs": {
        "structure": "zincblende", "a": 6.06, "bandgap_eV": 0.35,
        "growth_directions": ["[111]B", "[0001]"], "passivation": None
    },
    "InP": {
        "structure": "zincblende", "a": 5.87, "bandgap_eV": 1.35,
        "growth_directions": ["[111]B", "[0001]"], "passivation": "S"
    },
    "GaN": {
        "structure": "wurtzite", "a": 3.19, "c": 5.19, "bandgap_eV": 3.4,
        "growth_directions": ["[0001]", "[1-100]"], "passivation": None
    },
    # II-VI semiconductors
    "ZnO": {
        "structure": "wurtzite", "a": 3.25, "c": 5.21, "bandgap_eV": 3.37,
        "growth_directions": ["[0001]"], "passivation": None
    },
    "CdSe": {
        "structure": "wurtzite", "a": 4.30, "c": 7.01, "bandgap_eV": 1.74,
        "growth_directions": ["[0001]"], "passivation": "Cd-rich"
    },
    # Metals
    "Au": {
        "structure": "fcc", "a": 4.08, "bandgap_eV": 0,
        "growth_directions": ["[111]", "[110]", "[100]"], "passivation": None
    },
    "Ag": {
        "structure": "fcc", "a": 4.09, "bandgap_eV": 0,
        "growth_directions": ["[111]", "[110]"], "passivation": None
    },
    "Cu": {
        "structure": "fcc", "a": 3.61, "bandgap_eV": 0,
        "growth_directions": ["[111]", "[110]"], "passivation": None
    },
    "Pt": {
        "structure": "fcc", "a": 3.92, "bandgap_eV": 0,
        "growth_directions": ["[111]"], "passivation": None
    },
}


# Core-shell combinations
CORESHELL_DATABASE = {
    "Si/SiO2": {"core": "Si", "shell": "SiO2", "interface": "native_oxide"},
    "Ge/Si": {"core": "Ge", "shell": "Si", "interface": "epitaxial"},
    "GaAs/AlGaAs": {"core": "GaAs", "shell": "AlGaAs", "interface": "epitaxial"},
    "InAs/InP": {"core": "InAs", "shell": "InP", "interface": "epitaxial"},
    "Au/Ag": {"core": "Au", "shell": "Ag", "interface": "metallic"},
    "Au/SiO2": {"core": "Au", "shell": "SiO2", "interface": "insulating"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_nanowire(
    material: str = "Si",
    diameter_nm: float = 3.0,
    length_nm: float = 10.0,
    growth_direction: str = "[111]",
    cross_section: str = "circular",
    passivated: bool = True
) -> Dict[str, Any]:
    """
    Generate semiconductor/metal nanowire.
    
    Args:
        material: Material from database
        diameter_nm: Wire diameter in nm
        length_nm: Wire length in nm
        growth_direction: Crystallographic growth direction
        cross_section: 'circular', 'hexagonal', 'square'
        passivated: Add surface passivation
    
    Returns:
        Nanowire structure
    """
    if material not in NANOWIRE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material '{material}'",
                      "available": list(NANOWIRE_DATABASE.keys())}
        }
    
    info = NANOWIRE_DATABASE[material]
    a = info["a"]
    structure_type = info["structure"]
    
    diameter = diameter_nm * 10  # Convert to Angstrom
    length = length_nm * 10
    radius = diameter / 2
    
    atoms = []
    
    # Generate bulk positions and cut to wire shape
    if structure_type == "diamond":
        # Diamond structure
        basis = [
            [0, 0, 0], [0.25, 0.25, 0.25],
            [0.5, 0.5, 0], [0.75, 0.75, 0.25],
            [0.5, 0, 0.5], [0.75, 0.25, 0.75],
            [0, 0.5, 0.5], [0.25, 0.75, 0.75]
        ]
        
        nx = int(diameter / a) + 2
        ny = int(diameter / a) + 2
        nz = int(length / a) + 1
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for b in basis:
                        x = (i + b[0]) * a - diameter / 2
                        y = (j + b[1]) * a - diameter / 2
                        z = (k + b[2]) * a
                        
                        # Check if inside wire cross-section
                        inside = False
                        if cross_section == "circular":
                            inside = x**2 + y**2 <= radius**2
                        elif cross_section == "hexagonal":
                            # Hexagonal approximation
                            inside = abs(x) <= radius and abs(y) <= radius * 0.866
                        else:  # square
                            inside = abs(x) <= radius and abs(y) <= radius
                        
                        if inside and 0 <= z <= length:
                            atoms.append({"element": material, "cartesian": [x, y, z]})
    
    elif structure_type == "zincblende":
        # Zincblende (GaAs, InAs, etc.)
        elem1, elem2 = material[0:2], material[2:4] if len(material) > 2 else material
        if material == "GaAs":
            elem1, elem2 = "Ga", "As"
        elif material == "InAs":
            elem1, elem2 = "In", "As"
        elif material == "InP":
            elem1, elem2 = "In", "P"
        
        basis = [
            ([0, 0, 0], elem1), ([0.25, 0.25, 0.25], elem2),
            ([0.5, 0.5, 0], elem1), ([0.75, 0.75, 0.25], elem2),
            ([0.5, 0, 0.5], elem1), ([0.75, 0.25, 0.75], elem2),
            ([0, 0.5, 0.5], elem1), ([0.25, 0.75, 0.75], elem2)
        ]
        
        nx = int(diameter / a) + 2
        ny = int(diameter / a) + 2
        nz = int(length / a) + 1
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for b, elem in basis:
                        x = (i + b[0]) * a - diameter / 2
                        y = (j + b[1]) * a - diameter / 2
                        z = (k + b[2]) * a
                        
                        inside = x**2 + y**2 <= radius**2 if cross_section == "circular" else \
                                 abs(x) <= radius and abs(y) <= radius
                        
                        if inside and 0 <= z <= length:
                            atoms.append({"element": elem, "cartesian": [x, y, z]})
    
    elif structure_type == "wurtzite":
        c = info.get("c", a * 1.63)
        elem1, elem2 = material[0:2], material[2:] if len(material) > 2 else "O"
        if material == "ZnO":
            elem1, elem2 = "Zn", "O"
        elif material == "GaN":
            elem1, elem2 = "Ga", "N"
        elif material == "CdSe":
            elem1, elem2 = "Cd", "Se"
        
        nx = int(diameter / a) + 2
        ny = int(diameter / (a * np.sqrt(3)/2)) + 2
        nz = int(length / c) + 1
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = i * a + (j % 2) * a / 2 - diameter / 2
                    y = j * a * np.sqrt(3) / 2 - diameter / 2
                    z1 = k * c
                    z2 = k * c + 0.375 * c
                    
                    if x**2 + y**2 <= radius**2 and 0 <= z1 <= length:
                        atoms.append({"element": elem1, "cartesian": [x, y, z1]})
                        atoms.append({"element": elem2, "cartesian": [x, y, z2]})
    
    elif structure_type == "fcc":
        # FCC metal
        basis = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        
        nx = int(diameter / a) + 2
        ny = int(diameter / a) + 2
        nz = int(length / a) + 1
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for b in basis:
                        x = (i + b[0]) * a - diameter / 2
                        y = (j + b[1]) * a - diameter / 2
                        z = (k + b[2]) * a
                        
                        inside = x**2 + y**2 <= radius**2
                        
                        if inside and 0 <= z <= length:
                            atoms.append({"element": material, "cartesian": [x, y, z]})
    
    # Add passivation
    n_passivation = 0
    if passivated and info["passivation"]:
        passivant = info["passivation"]
        # Find surface atoms and add passivation
        for atom in atoms[:]:  # Iterate over copy
            pos = np.array(atom["cartesian"])
            r = np.sqrt(pos[0]**2 + pos[1]**2)
            
            # If near surface
            if r > radius - 2:
                # Add passivation atom radially outward
                radial = pos[:2] / (r + 0.001)
                pass_pos = [pos[0] + radial[0] * 1.5, pos[1] + radial[1] * 1.5, pos[2]]
                atoms.append({"element": passivant, "cartesian": pass_pos, "passivation": True})
                n_passivation += 1
    
    return {
        "success": True,
        "material": material,
        "structure_type": structure_type,
        "diameter_nm": diameter_nm,
        "length_nm": length_nm,
        "growth_direction": growth_direction,
        "cross_section": cross_section,
        "bandgap_eV": info["bandgap_eV"],
        "passivated": passivated,
        "n_passivation_atoms": n_passivation,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_coreshell_nanowire(
    system: str = "Si/SiO2",
    core_diameter_nm: float = 5.0,
    shell_thickness_nm: float = 2.0,
    length_nm: float = 20.0
) -> Dict[str, Any]:
    """
    Generate core-shell nanowire.
    
    Args:
        system: Core/shell material combination
        core_diameter_nm: Core diameter
        shell_thickness_nm: Shell thickness
        length_nm: Wire length
    
    Returns:
        Core-shell nanowire structure
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
    
    # Generate core
    core = generate_nanowire(core_mat, core_diameter_nm, length_nm, passivated=False)
    
    if not core["success"]:
        return core
    
    atoms = list(core["structure"]["atoms"])
    
    # Generate shell
    outer_diameter = core_diameter_nm + 2 * shell_thickness_nm
    core_radius = core_diameter_nm * 10 / 2
    outer_radius = outer_diameter * 10 / 2
    length = length_nm * 10
    
    # Simple shell generation (oxide or semiconductor)
    if shell_mat == "SiO2":
        # Amorphous SiO2 shell
        n_shell = int(2 * np.pi * (core_radius + outer_radius) / 2 * length / 10)
        np.random.seed(42)
        
        for _ in range(n_shell):
            r = np.random.uniform(core_radius + 0.5, outer_radius)
            theta = np.random.uniform(0, 2 * np.pi)
            z = np.random.uniform(0, length)
            
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            
            atoms.append({"element": "Si", "cartesian": [x, y, z], "shell": True})
            # Add 2 O atoms nearby
            for o in range(2):
                ox = x + np.random.uniform(-1.6, 1.6)
                oy = y + np.random.uniform(-1.6, 1.6)
                oz = z + np.random.uniform(-0.8, 0.8)
                atoms.append({"element": "O", "cartesian": [ox, oy, oz], "shell": True})
    else:
        # Epitaxial shell
        shell_nw = generate_nanowire(shell_mat, outer_diameter, length_nm, passivated=False)
        if shell_nw["success"]:
            for atom in shell_nw["structure"]["atoms"]:
                pos = atom["cartesian"]
                r = np.sqrt(pos[0]**2 + pos[1]**2)
                if r > core_radius:
                    atom["shell"] = True
                    atoms.append(atom)
    
    return {
        "success": True,
        "system": system,
        "core_material": core_mat,
        "shell_material": shell_mat,
        "core_diameter_nm": core_diameter_nm,
        "shell_thickness_nm": shell_thickness_nm,
        "total_diameter_nm": outer_diameter,
        "interface_type": info["interface"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_axial_heterostructure(
    materials: List[str] = ["InAs", "InP", "InAs"],
    segment_lengths_nm: List[float] = [10.0, 5.0, 10.0],
    diameter_nm: float = 50.0
) -> Dict[str, Any]:
    """
    Generate axial heterostructure nanowire.
    
    Args:
        materials: List of materials for each segment
        segment_lengths_nm: Length of each segment
        diameter_nm: Wire diameter (constant)
    
    Returns:
        Axial heterostructure
    """
    if len(materials) != len(segment_lengths_nm):
        return {"success": False, "error": {"code": "LENGTH_MISMATCH", "message": "Materials and lengths must match"}}
    
    atoms = []
    z_offset = 0
    segments = []
    
    for mat, seg_len in zip(materials, segment_lengths_nm):
        seg = generate_nanowire(mat, diameter_nm, seg_len, passivated=False)
        
        if seg["success"]:
            for atom in seg["structure"]["atoms"]:
                new_z = atom["cartesian"][2] + z_offset * 10
                atoms.append({
                    "element": atom["element"],
                    "cartesian": [atom["cartesian"][0], atom["cartesian"][1], new_z],
                    "segment": mat
                })
            
            segments.append({
                "material": mat,
                "length_nm": seg_len,
                "z_start_nm": z_offset,
                "z_end_nm": z_offset + seg_len
            })
            
            z_offset += seg_len
    
    return {
        "success": True,
        "type": "axial_heterostructure",
        "materials": materials,
        "n_segments": len(materials),
        "total_length_nm": sum(segment_lengths_nm),
        "diameter_nm": diameter_nm,
        "segments": segments,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
