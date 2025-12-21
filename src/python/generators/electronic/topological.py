"""
electronic/topological.py - Topological Insulators and Semimetals

Comprehensive topological materials per structure_catalogue.md Category 8:
(viii) Topological insulator nano-structures (Bi2Se3, Sb2Te3)
(ix) Majorana nanowires (InSb, InAs) with epitaxial Al shell
Plus: Weyl/Dirac semimetals, nodal-line materials
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Topological insulator database
TOPOLOGICAL_INSULATOR_DATABASE = {
    # 3D TIs
    "Bi2Se3": {
        "class": "3D_TI", "structure": "rhombohedral", "spacegroup": "R-3m",
        "a": 4.14, "c": 28.64, "bandgap_eV": 0.3,
        "surface_state": "single_Dirac_cone", "spin_texture": "helical"
    },
    "Bi2Te3": {
        "class": "3D_TI", "structure": "rhombohedral", "spacegroup": "R-3m",
        "a": 4.38, "c": 30.49, "bandgap_eV": 0.15,
        "surface_state": "single_Dirac_cone", "thermoelectric": True
    },
    "Sb2Te3": {
        "class": "3D_TI", "structure": "rhombohedral", "spacegroup": "R-3m",
        "a": 4.26, "c": 30.35, "bandgap_eV": 0.21,
        "surface_state": "single_Dirac_cone", "applications": ["memory"]
    },
    "Bi2Te2Se": {
        "class": "3D_TI", "structure": "rhombohedral", "spacegroup": "R-3m",
        "a": 4.30, "c": 29.86, "bandgap_eV": 0.25,
        "surface_state": "single", "ordered_ternary": True
    },
    "SmB6": {
        "class": "Kondo_TI", "structure": "cubic", "spacegroup": "Pm-3m",
        "a": 4.13, "bandgap_eV": 0.02,
        "surface_state": "metallic", "strong_correlation": True
    },
    
    # 2D TIs
    "WTe2_1Tp": {
        "class": "2D_TI", "structure": "monoclinic", 
        "a": 6.28, "b": 3.48, "c": 6.28, "bandgap_eV": 0.05,
        "edge_state": "helical", "type": "quantum_spin_Hall"
    },
    "Bi_monolayer": {
        "class": "2D_TI", "structure": "hexagonal",
        "a": 4.54, "bandgap_eV": 0.4,
        "edge_state": "helical", "buckling": 1.6
    },
}


# Weyl and Dirac semimetal database
SEMIMETAL_DATABASE = {
    # Weyl semimetals
    "TaAs": {
        "class": "Weyl_type_I", "structure": "tetragonal", "spacegroup": "I41md",
        "a": 3.44, "c": 11.64, "n_Weyl_points": 24,
        "Fermi_arc": True, "inversion_breaking": True
    },
    "NbAs": {
        "class": "Weyl_type_I", "structure": "tetragonal", "spacegroup": "I41md",
        "a": 3.45, "c": 11.68, "n_Weyl_points": 24
    },
    "WTe2": {
        "class": "Weyl_type_II", "structure": "orthorhombic", "spacegroup": "Pmn21",
        "a": 3.48, "b": 6.27, "c": 14.07, "tilted_cones": True
    },
    "MoTe2": {
        "class": "Weyl_type_II", "structure": "monoclinic",
        "a": 6.33, "b": 3.47, "c": 13.86, "superconducting": True
    },
    
    # Dirac semimetals
    "Cd3As2": {
        "class": "Dirac_3D", "structure": "tetragonal", "spacegroup": "I41/acd",
        "a": 12.67, "c": 25.48, "n_Dirac_points": 2,
        "linear_dispersion": True, "protected_by": "C4_rotation"
    },
    "Na3Bi": {
        "class": "Dirac_3D", "structure": "hexagonal", "spacegroup": "P63/mmc",
        "a": 5.45, "c": 9.66, "n_Dirac_points": 2
    },
    "ZrTe5": {
        "class": "Dirac_quasi_2D", "structure": "orthorhombic",
        "a": 3.99, "b": 14.5, "c": 13.7, "layered": True
    },
    
    # Nodal-line semimetals
    "CaP3": {
        "class": "nodal_line", "structure": "monoclinic",
        "a": 5.57, "b": 5.57, "c": 9.15, "nodal_ring": True
    },
    "Cu3PdN": {
        "class": "nodal_line", "structure": "cubic", "spacegroup": "Pm-3m",
        "a": 3.85, "protected_by": "mirror"
    },
}


# Majorana platform materials
MAJORANA_DATABASE = {
    "InSb_Al": {
        "semiconductor": "InSb", "superconductor": "Al",
        "g_factor": 50, "spin_orbit_meV": 1.0,
        "gap_induced_meV": 0.2, "optimal_field_T": 0.3
    },
    "InAs_Al": {
        "semiconductor": "InAs", "superconductor": "Al",
        "g_factor": 15, "spin_orbit_meV": 0.2,
        "gap_induced_meV": 0.25, "optimal_field_T": 0.5
    },
    "InAs_epitaxial_Al": {
        "semiconductor": "InAs", "superconductor": "Al",
        "g_factor": 8, "spin_orbit_meV": 0.3,
        "gap_induced_meV": 0.2, "epitaxial": True, "hard_gap": True
    },
    "HgTe_QW": {
        "semiconductor": "HgTe", "superconductor": "Nb",
        "quantum_well": True, "inverted_band": True
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_topological_insulator(
    material: str = None,
    thickness_QL: int = 6,
    supercell: List[int] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate topological insulator structure.

    Args:
        material: TI material from database
        thickness_QL: Thickness in quintuple layers (for thin films)
        supercell: Supercell dimensions
        **kwargs: Accepts spacegroup, elements, composition, a, c for custom structures

    Returns:
        TI structure
    """
    if supercell is None:
        supercell = kwargs.get('supercell', [2, 2, 1])

    # Handle custom structure specification via space group
    if kwargs.get('elements') or kwargs.get('spacegroup'):
        elements = kwargs.get('elements', ['Bi', 'Se'])
        composition = kwargs.get('composition', [2, 3])
        a_param = kwargs.get('a', 4.14)
        c_param = kwargs.get('c', 28.64)
        sg = kwargs.get('spacegroup', 166)

        # R-3m (166) rhombohedral structure for Bi2Se3-type TIs
        if sg == 166:
            # Create hexagonal setting of R-3m
            lattice = Lattice.hexagonal(a_param, c_param)

            # Bi2Se3 has 5 atoms per formula unit in quintuple layer
            # Use proper Wyckoff positions for R-3m
            species = []
            coords = []

            # Simplified model: create quintuple layer structure
            # Se-Bi-Se-Bi-Se pattern along c-axis
            for i, (elem, count) in enumerate(zip(elements, composition)):
                for j in range(count):
                    z_frac = (i * 3 + j) / (sum(composition) * 3)
                    species.append(elem)
                    coords.append([0.0, 0.0, z_frac])

            structure = Structure(lattice, species, coords)

            return {
                "success": True,
                "material": f"custom_{''.join(elements)}",
                "spacegroup": sg,
                "n_atoms": len(structure),
                "topological_class": "3D_TI",
                "bandgap_eV": 0.3,
                "surface_state": "single_Dirac_cone",
                "structure": structure_to_dict(structure)
            }

        return {
            "success": False,
            "error": {"code": "UNSUPPORTED_SG", "message": f"Space group {sg} not yet supported",
                      "supported": [166]}
        }

    # Default to Bi2Se3
    if material is None:
        material = "Bi2Se3"

    if material not in TOPOLOGICAL_INSULATOR_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown TI '{material}'",
                      "available": list(TOPOLOGICAL_INSULATOR_DATABASE.keys())}
        }
    
    info = TOPOLOGICAL_INSULATOR_DATABASE[material]
    
    if info["structure"] == "rhombohedral":
        a = info["a"]
        c = info["c"]
        lattice = Lattice.hexagonal(a, c / 3 * thickness_QL)
        
        species = []
        coords = []
        
        # Quintuple layer structure: Se-Bi-Se-Bi-Se
        if "Bi2Se3" in material:
            elements = ["Se", "Bi", "Se", "Bi", "Se"]
            z_positions = [0, 0.2, 0.4, 0.6, 0.8]
        elif "Bi2Te3" in material:
            elements = ["Te", "Bi", "Te", "Bi", "Te"]
            z_positions = [0, 0.2, 0.4, 0.6, 0.8]
        elif "Sb2Te3" in material:
            elements = ["Te", "Sb", "Te", "Sb", "Te"]
            z_positions = [0, 0.2, 0.4, 0.6, 0.8]
        else:
            elements = ["X"] * 5
            z_positions = [0.1 * i for i in range(5)]
        
        nx, ny, nz = supercell
        
        for ql in range(thickness_QL):
            for i in range(nx):
                for j in range(ny):
                    for k, (elem, z_pos) in enumerate(zip(elements, z_positions)):
                        x = (i + (k % 3) * 0.33) / nx
                        y = (j + (k % 3) * 0.33) / ny
                        z = (ql + z_pos) / thickness_QL
                        
                        species.append(elem)
                        coords.append([x % 1, y % 1, z])
        
    elif info["structure"] == "cubic":
        a = info["a"]
        lattice = Lattice.cubic(a * supercell[0])
        
        # SmB6 structure
        species = ["Sm"]
        coords = [[0, 0, 0]]
        
        # B6 octahedra
        b_positions = [
            [0.2, 0.5, 0.5], [0.8, 0.5, 0.5],
            [0.5, 0.2, 0.5], [0.5, 0.8, 0.5],
            [0.5, 0.5, 0.2], [0.5, 0.5, 0.8]
        ]
        for bp in b_positions:
            species.append("B")
            coords.append(bp)
    
    else:
        a = info.get("a", 5.0)
        lattice = Lattice.cubic(a)
        species = [material[:2]]
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "class": info["class"],
        "bandgap_eV": info.get("bandgap_eV", 0),
        "surface_state": info.get("surface_state", ""),
        "spin_texture": info.get("spin_texture", ""),
        "thickness_QL": thickness_QL,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_weyl_semimetal(
    material: str = "TaAs",
    supercell: List[int] = [2, 2, 1]
) -> Dict[str, Any]:
    """
    Generate Weyl or Dirac semimetal structure.
    
    Args:
        material: Semimetal from database
        supercell: Supercell dimensions
    
    Returns:
        Semimetal structure
    """
    if material not in SEMIMETAL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown semimetal",
                      "available": list(SEMIMETAL_DATABASE.keys())}
        }
    
    info = SEMIMETAL_DATABASE[material]
    
    if info["structure"] == "tetragonal":
        a = info["a"]
        c = info.get("c", a * 3)
        lattice = Lattice.tetragonal(a, c)
        
        # TaAs structure (I41md)
        species = []
        coords = []
        
        if "TaAs" in material or "NbAs" in material:
            metal = "Ta" if "Ta" in material else "Nb"
            species = [metal, metal, "As", "As"]
            coords = [
                [0, 0, 0], [0.5, 0.5, 0.5],
                [0, 0.5, 0.25], [0.5, 0, 0.75]
            ]
    else:
        a = info.get("a", 5.0)
        b = info.get("b", a)
        c = info.get("c", a)
        lattice = Lattice.orthorhombic(a, b, c)
        species = [material[:2]]
        coords = [[0, 0, 0]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "material": material,
        "class": info["class"],
        "n_Weyl_points": info.get("n_Weyl_points", 0),
        "n_Dirac_points": info.get("n_Dirac_points", 0),
        "has_Fermi_arc": info.get("Fermi_arc", False),
        "tilted_cones": info.get("tilted_cones", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_majorana_nanowire(
    platform: str = "InSb_Al",
    wire_length_nm: float = 1000,
    wire_diameter_nm: float = 100,
    al_shell_thickness_nm: float = 10,
    tunnel_barrier: bool = True
) -> Dict[str, Any]:
    """
    Generate Majorana nanowire structure.
    
    Args:
        platform: Semiconductor/superconductor platform
        wire_length_nm: Wire length
        wire_diameter_nm: Wire diameter
        al_shell_thickness_nm: Al shell thickness
        tunnel_barrier: Include tunnel barrier
    
    Returns:
        Majorana nanowire structure
    """
    if platform not in MAJORANA_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_PLATFORM", "message": f"Unknown platform",
                      "available": list(MAJORANA_DATABASE.keys())}
        }
    
    info = MAJORANA_DATABASE[platform]
    
    # Create wire geometry (simplified)
    length = wire_length_nm * 10  # Angstrom
    diameter = wire_diameter_nm * 10
    shell_thickness = al_shell_thickness_nm * 10
    
    lattice = Lattice.orthorhombic(diameter + 2 * shell_thickness, 
                                    diameter + 2 * shell_thickness, 
                                    min(length, 500))  # Cap length for structure
    
    species = []
    coords = []
    
    semiconductor = info["semiconductor"]
    superconductor = info["superconductor"]
    
    # Semiconductor core
    n_core = 50
    radius = 0.3
    for i in range(n_core):
        z = (i + 0.5) / n_core
        for j in range(6):
            angle = 2 * np.pi * j / 6
            r = radius * 0.8
            x = 0.5 + r * np.cos(angle)
            y = 0.5 + r * np.sin(angle)
            
            if "InSb" in semiconductor:
                species.extend(["In", "Sb"])
                coords.extend([[x, y, z], [x + 0.02, y, z]])
            else:
                species.extend(["In", "As"])
                coords.extend([[x, y, z], [x + 0.02, y, z]])
    
    # Superconductor shell (partial coverage for gate access)
    n_shell = 30
    shell_radius = 0.4
    for i in range(n_shell):
        z = (i + 0.5) / n_shell
        for j in range(4):  # Partial shell (not full coverage)
            angle = np.pi * j / 4 - np.pi / 4
            x = 0.5 + shell_radius * np.cos(angle)
            y = 0.5 + shell_radius * np.sin(angle)
            
            species.append(superconductor)
            coords.append([x, y, z])
    
    # Tunnel barrier (if requested)
    if tunnel_barrier:
        # Add thin oxide layer at ends
        for end_z in [0.02, 0.98]:
            for j in range(4):
                angle = 2 * np.pi * j / 4
                x = 0.5 + radius * 0.5 * np.cos(angle)
                y = 0.5 + radius * 0.5 * np.sin(angle)
                species.extend(["Al", "O"])
                coords.extend([[x, y, end_z], [x + 0.01, y, end_z]])
    
    structure = Structure(lattice, species, coords)
    
    # Estimate critical magnetic field
    B_c = info.get("optimal_field_T", 0.3)
    
    return {
        "success": True,
        "platform": platform,
        "semiconductor": semiconductor,
        "superconductor": superconductor,
        "wire_length_nm": wire_length_nm,
        "wire_diameter_nm": wire_diameter_nm,
        "shell_thickness_nm": al_shell_thickness_nm,
        "has_tunnel_barrier": tunnel_barrier,
        "g_factor": info.get("g_factor", 0),
        "spin_orbit_meV": info.get("spin_orbit_meV", 0),
        "induced_gap_meV": info.get("gap_induced_meV", 0),
        "optimal_field_T": B_c,
        "is_epitaxial": info.get("epitaxial", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
