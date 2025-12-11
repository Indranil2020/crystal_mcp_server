#!/usr/bin/env python3
"""
Advanced Structure Generator

This module provides advanced structure generation capabilities including:
- Twisted bilayers (graphene, TMDs)
- Prototype structures (Perovskite, Spinel, Heusler, etc.)
- Metal-Organic Frameworks (MOFs)
- High-entropy alloys (SQS)
- 2D materials (MXenes, h-BN, phosphorene)
- Cage structures (clathrates, zeolites)

All functions follow defensive programming without try/except blocks.
"""

from typing import Dict, Any, List, Optional, Tuple
import json
import sys
import numpy as np
from pymatgen.core import Structure, Lattice, Element, Composition
from pymatgen.transformations.advanced_transformations import SupercellTransformation
from ase import Atoms
from ase.build import graphene_nanoribbon, nanotube, molecule
from ase.lattice.hexagonal import Graphene


# ============================================================================
# PROTOTYPE STRUCTURES
# ============================================================================

PROTOTYPE_TEMPLATES = {
    "rocksalt": {
        "space_group": 225,
        "wyckoff": {"A": "4a", "B": "4b"},
        "coords": {"A": [0, 0, 0], "B": [0.5, 0.5, 0.5]},
        "lattice_type": "cubic",
        "description": "NaCl-type (Fm-3m)"
    },
    "zincblende": {
        "space_group": 216,
        "wyckoff": {"A": "4a", "B": "4c"},
        "coords": {"A": [0, 0, 0], "B": [0.25, 0.25, 0.25]},
        "lattice_type": "cubic",
        "description": "ZnS-type (F-43m)"
    },
    "wurtzite": {
        "space_group": 186,
        "wyckoff": {"A": "2b", "B": "2b"},
        "coords": {"A": [1/3, 2/3, 0], "B": [1/3, 2/3, 0.375]},
        "lattice_type": "hexagonal",
        "description": "ZnO-type (P6_3mc)"
    },
    "fluorite": {
        "space_group": 225,
        "wyckoff": {"A": "4a", "B": "8c"},
        "coords": {"A": [0, 0, 0], "B": [0.25, 0.25, 0.25]},
        "lattice_type": "cubic",
        "description": "CaF2-type (Fm-3m)"
    },
    "antifluorite": {
        "space_group": 225,
        "wyckoff": {"A": "8c", "B": "4a"},
        "coords": {"A": [0.25, 0.25, 0.25], "B": [0, 0, 0]},
        "lattice_type": "cubic",
        "description": "Li2O-type (Fm-3m)"
    },
    "perovskite": {
        "space_group": 221,
        "wyckoff": {"A": "1a", "B": "1b", "X": "3c"},
        "coords": {"A": [0, 0, 0], "B": [0.5, 0.5, 0.5], "X": [0.5, 0.5, 0]},
        "lattice_type": "cubic",
        "description": "ABX3 perovskite (Pm-3m)"
    },
    "spinel": {
        "space_group": 227,
        "wyckoff": {"A": "8a", "B": "16d", "X": "32e"},
        "coords": {"A": [0.125, 0.125, 0.125], "B": [0.5, 0.5, 0.5], "X": [0.2625, 0.2625, 0.2625]},
        "lattice_type": "cubic",
        "description": "AB2X4 spinel (Fd-3m)"
    },
    "heusler": {
        "space_group": 225,
        "wyckoff": {"X": "4a", "Y": "4b", "Z": "8c"},
        "coords": {"X": [0, 0, 0], "Y": [0.5, 0.5, 0.5], "Z": [0.25, 0.25, 0.25]},
        "lattice_type": "cubic",
        "description": "X2YZ Heusler (Fm-3m)"
    },
    "rutile": {
        "space_group": 136,
        "wyckoff": {"M": "2a", "X": "4f"},
        "coords": {"M": [0, 0, 0], "X": [0.305, 0.305, 0]},
        "lattice_type": "tetragonal",
        "description": "TiO2-type rutile (P4_2/mnm)"
    },
    "diamond": {
        "space_group": 227,
        "wyckoff": {"A": "8a"},
        "coords": {"A": [0.125, 0.125, 0.125]},
        "lattice_type": "cubic",
        "description": "Diamond (Fd-3m)"
    }
}


def generate_prototype_structure(
    prototype: str,
    elements: Dict[str, str],
    lattice_constant: float,
    c_over_a: float = 1.0
) -> Dict[str, Any]:
    """
    Generate a structure from a prototype template.
    
    Args:
        prototype: Prototype name (e.g., 'perovskite', 'rocksalt', 'spinel')
        elements: Mapping of site labels to elements (e.g., {'A': 'Ca', 'B': 'Ti', 'X': 'O'})
        lattice_constant: Lattice constant 'a' in Angstroms
        c_over_a: c/a ratio for non-cubic systems
    
    Returns:
        Dictionary with structure data
    """
    if prototype not in PROTOTYPE_TEMPLATES:
        available = list(PROTOTYPE_TEMPLATES.keys())
        return {
            "success": False,
            "error": {
                "code": "INVALID_PROTOTYPE",
                "message": f"Unknown prototype '{prototype}'",
                "available": available
            }
        }
    
    template = PROTOTYPE_TEMPLATES[prototype]
    
    # Validate elements match template requirements
    required_sites = set(template["wyckoff"].keys())
    provided_sites = set(elements.keys())
    if required_sites != provided_sites:
        return {
            "success": False,
            "error": {
                "code": "INVALID_ELEMENTS",
                "message": f"Prototype '{prototype}' requires sites {required_sites}, got {provided_sites}"
            }
        }
    
    # Create lattice
    a = lattice_constant
    if template["lattice_type"] == "cubic":
        lattice = Lattice.cubic(a)
    elif template["lattice_type"] == "hexagonal":
        c = a * c_over_a
        lattice = Lattice.hexagonal(a, c)
    elif template["lattice_type"] == "tetragonal":
        c = a * c_over_a
        lattice = Lattice.tetragonal(a, c)
    else:
        lattice = Lattice.cubic(a)
    
    # Build species and coordinates
    species = []
    coords = []
    
    for site_label, element in elements.items():
        site_coords = template["coords"][site_label]
        # Handle multiplicity (simplified - assumes single site per label)
        species.append(element)
        coords.append(site_coords)
    
    # Create structure
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "prototype": prototype,
        "description": template["description"],
        "space_group": template["space_group"],
        "structure": structure_to_dict(structure)
    }


# ============================================================================
# TWISTED BILAYER STRUCTURES
# ============================================================================

def generate_twisted_bilayer(
    material: str = "graphene",
    twist_angle: float = 1.1,
    layers: int = 2,
    stacking: str = "AB",
    interlayer_distance: float = 3.35,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate a twisted bilayer structure.
    
    Args:
        material: Base material ('graphene', 'MoS2', 'hBN')
        twist_angle: Twist angle in degrees
        layers: Number of layers (default 2)
        stacking: Stacking type ('AA', 'AB')
        interlayer_distance: Distance between layers in Angstroms
        vacuum: Vacuum padding in Angstroms
    
    Returns:
        Dictionary with twisted structure data
    """
    # For small twist angles, find commensurate supercell
    # Using simplified formula for graphene magic angles
    
    if material.lower() == "graphene":
        a = 2.46  # Graphene lattice constant
        bond = 1.42  # C-C bond length
    elif material.lower() in ["mos2", "ws2"]:
        a = 3.16
        bond = None
    elif material.lower() == "hbn":
        a = 2.50
        bond = 1.45
    else:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Material '{material}' not supported for twisted bilayer"}
        }
    
    # Calculate commensurate indices for twist angle
    # For angle θ, use (m, n) indices where cos(θ) = (m² + n² + 4mn) / (2(m² + n² + mn))
    # Simplified: find approximate commensurate cell
    
    theta_rad = np.radians(twist_angle)
    
    # Heuristic: find m, n for common twist angles
    # For magic angle ~1.1°, use (31, 30) or similar large indices
    # For 30°, use (1, 1)
    # For 21.8°, use (2, 1)
    
    # Simplification: Generate small twisted cell for demonstration
    # Real implementation would need larger commensurate cells
    
    if abs(twist_angle - 30.0) < 1.0:
        m, n = 1, 1
    elif abs(twist_angle - 21.8) < 1.0:
        m, n = 2, 1
    elif abs(twist_angle - 13.2) < 1.0:
        m, n = 3, 1
    else:
        # For arbitrary angles, estimate size
        m = max(1, int(60 / twist_angle))
        n = m - 1 if m > 1 else 1
    
    # Create bottom layer (untwisted)
    atoms_bottom = Graphene(symbol='C', latticeconstant={'a': a, 'c': interlayer_distance * 2}, 
                            size=(m + n, m + n, 1))
    
    # Create top layer (will be rotated)
    atoms_top = Graphene(symbol='C', latticeconstant={'a': a, 'c': interlayer_distance * 2}, 
                         size=(m + n, m + n, 1))
    
    # Rotate top layer
    atoms_top.rotate(twist_angle, 'z', center='COP')
    
    # Shift top layer up
    atoms_top.translate([0, 0, interlayer_distance])
    
    # Combine
    combined = atoms_bottom + atoms_top
    
    # Center and add vacuum
    combined.center(vacuum=vacuum, axis=2)
    
    return {
        "success": True,
        "material": material,
        "twist_angle": twist_angle,
        "commensurate_indices": (m, n),
        "n_atoms": len(combined),
        "structure": atoms_to_dict(combined),
        "warnings": ["Commensurate cell is approximate. For exact magic angle cells, use specialized tools."]
    }


# ============================================================================
# HIGH-ENTROPY ALLOYS
# ============================================================================

def generate_high_entropy_alloy(
    elements: List[str],
    concentrations: Optional[List[float]] = None,
    structure_type: str = "fcc",
    supercell: List[int] = [3, 3, 3],
    lattice_constant: Optional[float] = None,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate a high-entropy alloy structure using random (SQS-like) approach.
    
    Args:
        elements: List of 4+ elements (e.g., ['Fe', 'Co', 'Ni', 'Cr', 'Mn'])
        concentrations: Optional concentrations (default: equimolar)
        structure_type: Base structure ('fcc', 'bcc', 'hcp')
        supercell: Supercell size
        lattice_constant: Lattice constant (estimated if not provided)
        seed: Random seed for reproducibility
    
    Returns:
        Dictionary with HEA structure data
    """
    if len(elements) < 4:
        return {
            "success": False,
            "error": {"code": "INVALID_INPUT", "message": "High-entropy alloys require at least 4 elements"}
        }
    
    if seed is not None:
        np.random.seed(seed)
    
    # Default to equimolar
    if concentrations is None:
        concentrations = [1.0 / len(elements)] * len(elements)
    
    # Normalize concentrations
    total = sum(concentrations)
    concentrations = [c / total for c in concentrations]
    
    # Estimate lattice constant from average atomic radii if not provided
    if lattice_constant is None:
        # Average metallic radii (simplified)
        radii = {
            'Fe': 1.26, 'Co': 1.25, 'Ni': 1.24, 'Cr': 1.28, 'Mn': 1.27,
            'Cu': 1.28, 'Al': 1.43, 'Ti': 1.47, 'V': 1.34, 'Nb': 1.46,
            'Mo': 1.39, 'W': 1.39, 'Ta': 1.46, 'Zr': 1.60, 'Hf': 1.59
        }
        avg_radius = sum(radii.get(e, 1.35) * c for e, c in zip(elements, concentrations))
        
        if structure_type == "fcc":
            lattice_constant = 2 * np.sqrt(2) * avg_radius
        elif structure_type == "bcc":
            lattice_constant = 4 * avg_radius / np.sqrt(3)
        else:  # hcp
            lattice_constant = 2 * avg_radius
    
    # Create base structure
    if structure_type == "fcc":
        lattice = Lattice.cubic(lattice_constant)
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    elif structure_type == "bcc":
        lattice = Lattice.cubic(lattice_constant)
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    else:  # hcp
        c = lattice_constant * 1.633  # Ideal c/a
        lattice = Lattice.hexagonal(lattice_constant, c)
        coords = [[1/3, 2/3, 0.25], [2/3, 1/3, 0.75]]
    
    # Create unit cell with placeholder element
    unit_species = [elements[0]] * len(coords)
    unit_cell = Structure(lattice, unit_species, coords)
    
    # Create supercell
    nx, ny, nz = supercell
    sc = unit_cell.make_supercell([nx, ny, nz])
    
    # Random assignment based on concentrations
    n_atoms = len(sc)
    n_per_element = [int(round(c * n_atoms)) for c in concentrations]
    
    # Adjust to match total
    diff = n_atoms - sum(n_per_element)
    if diff > 0:
        for i in range(diff):
            n_per_element[i % len(n_per_element)] += 1
    elif diff < 0:
        for i in range(-diff):
            n_per_element[i % len(n_per_element)] -= 1
    
    # Create assignment list
    assignment = []
    for elem, count in zip(elements, n_per_element):
        assignment.extend([elem] * count)
    
    # Shuffle
    np.random.shuffle(assignment)
    
    # Replace species
    for i, elem in enumerate(assignment):
        sc.replace(i, elem)
    
    # Calculate composition statistics
    comp = Composition(sc.formula)
    
    return {
        "success": True,
        "elements": elements,
        "concentrations": concentrations,
        "structure_type": structure_type,
        "n_atoms": n_atoms,
        "actual_composition": {e: n for e, n in zip(elements, n_per_element)},
        "lattice_constant": lattice_constant,
        "structure": structure_to_dict(sc)
    }


# ============================================================================
# 2D MATERIALS
# ============================================================================

def generate_2d_material(
    material: str,
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate various 2D materials.
    
    Args:
        material: Material type ('MXene', 'hBN', 'phosphorene', 'silicene', 'MoS2', 'WS2')
        size: Supercell size [nx, ny, 1]
        vacuum: Vacuum padding in Angstroms
        **kwargs: Material-specific parameters
    
    Returns:
        Dictionary with 2D material structure
    """
    material_lower = material.lower()
    
    if material_lower == "hbn":
        # Hexagonal boron nitride
        a = kwargs.get("a", 2.50)
        atoms = Graphene(symbol='B', latticeconstant={'a': a, 'c': vacuum * 2}, size=size)
        # Replace alternating atoms with N
        positions = atoms.get_scaled_positions()
        symbols = list(atoms.get_chemical_symbols())
        for i, pos in enumerate(positions):
            if (round(pos[0] * 3) + round(pos[1] * 3)) % 2 == 1:
                symbols[i] = 'N'
        atoms.set_chemical_symbols(symbols)
        
    elif material_lower in ["mos2", "ws2", "mose2", "wse2"]:
        from ase.build import mx2
        formula = material
        kind = kwargs.get("kind", "2H")
        a = kwargs.get("a", 3.16 if "Mo" in material else 3.18)
        thickness = kwargs.get("thickness", 3.19)
        atoms = mx2(formula=formula, kind=kind, a=a, thickness=thickness, size=tuple(size), vacuum=vacuum)
        
    elif material_lower == "phosphorene":
        # Black phosphorus monolayer
        a, b = 3.31, 4.38
        lattice = Lattice.orthorhombic(a, b, vacuum * 2)
        # Puckered structure
        coords = [
            [0.0, 0.0, 0.5 - 0.01],
            [0.0, 0.5, 0.5 + 0.01],
            [0.5, 0.0, 0.5 + 0.01],
            [0.5, 0.5, 0.5 - 0.01]
        ]
        structure = Structure(lattice, ['P'] * 4, coords)
        if size != [1, 1, 1]:
            structure.make_supercell(size)
        return {
            "success": True,
            "material": "phosphorene",
            "structure": structure_to_dict(structure)
        }
        
    elif material_lower == "silicene":
        # Buckled honeycomb silicon
        a = kwargs.get("a", 3.86)
        buckling = kwargs.get("buckling", 0.44)
        atoms = Graphene(symbol='Si', latticeconstant={'a': a, 'c': vacuum * 2}, size=size)
        # Add buckling
        positions = atoms.get_positions()
        for i, pos in enumerate(positions):
            if i % 2 == 1:
                positions[i, 2] += buckling
        atoms.set_positions(positions)
        
    elif "mxene" in material_lower:
        # Ti3C2 MXene as example
        a = kwargs.get("a", 3.07)
        return generate_mxene(
            metal=kwargs.get("metal", "Ti"),
            n_metal_layers=kwargs.get("n_metal_layers", 3),
            termination=kwargs.get("termination", "O"),
            a=a,
            size=size,
            vacuum=vacuum
        )
    else:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"2D material '{material}' not recognized"}
        }
    
    return {
        "success": True,
        "material": material,
        "n_atoms": len(atoms) if 'atoms' in dir() else 0,
        "structure": atoms_to_dict(atoms) if 'atoms' in dir() else None
    }


def generate_mxene(
    metal: str = "Ti",
    n_metal_layers: int = 3,
    termination: str = "O",
    a: float = 3.07,
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate MXene structure (M_{n+1}X_n T_2).
    
    Args:
        metal: Metal element (Ti, V, Nb, etc.)
        n_metal_layers: Number of metal layers (2, 3, or 4)
        termination: Surface termination ('O', 'OH', 'F', 'bare')
        a: Lattice constant
        size: Supercell size
        vacuum: Vacuum amount
    
    Returns:
        MXene structure dictionary
    """
    c = vacuum * 2
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    # Metal layers
    z_spacing = 2.3 / c  # Approximate interlayer spacing in fractional
    z_center = 0.5
    
    for i in range(n_metal_layers):
        z = z_center + (i - n_metal_layers / 2) * z_spacing
        species.append(metal)
        coords.append([1/3, 2/3, z])
    
    # Carbon layers (between metals)
    for i in range(n_metal_layers - 1):
        z = z_center + (i - n_metal_layers / 2 + 0.5) * z_spacing
        species.append('C')
        coords.append([0, 0, z])
    
    # Terminations
    if termination.lower() != "bare":
        term_elem = 'O' if termination.upper() in ['O', 'OH'] else termination
        z_top = z_center + (n_metal_layers / 2) * z_spacing + 0.02
        z_bot = z_center - (n_metal_layers / 2) * z_spacing - 0.02
        species.extend([term_elem, term_elem])
        coords.extend([[2/3, 1/3, z_top], [2/3, 1/3, z_bot]])
    
    structure = Structure(lattice, species, coords)
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "material": f"{metal}{n_metal_layers}C{n_metal_layers-1}{termination}2",
        "structure": structure_to_dict(structure)
    }


# ============================================================================
# CAGE STRUCTURES
# ============================================================================

def generate_cage_structure(
    cage_type: str,
    guest: Optional[str] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate cage-like structures (fullerenes, clathrates, zeolites).
    
    Args:
        cage_type: Type of cage ('C60', 'C70', 'clathrate_I', 'clathrate_II', 'SOD', 'LTA')
        guest: Optional guest atom/molecule for endohedral structures
        **kwargs: Additional parameters
    
    Returns:
        Cage structure dictionary
    """
    cage_lower = cage_type.lower()
    
    if cage_lower in ['c60', 'c70', 'c80', 'c84']:
        # Fullerenes from ASE
        atoms = molecule(cage_type.upper())
        atoms.center(vacuum=kwargs.get("vacuum", 5.0))
        
        if guest:
            # Endohedral: add guest at center
            center = atoms.get_center_of_mass()
            from ase import Atom
            atoms.append(Atom(guest, position=center))
        
        return {
            "success": True,
            "cage_type": cage_type,
            "guest": guest,
            "n_atoms": len(atoms),
            "structure": atoms_to_dict(atoms)
        }
    
    elif "clathrate" in cage_lower:
        # Type I clathrate (Pm-3n, SG 223)
        # Simplified model: Si46 framework
        a = kwargs.get("a", 10.2)
        lattice = Lattice.cubic(a)
        
        # Type I has 46 Si atoms per unit cell
        # Wyckoff positions: 6c, 16i, 24k
        coords_6c = [[0.25, 0, 0.5], [0.75, 0, 0.5], [0.5, 0.25, 0], 
                     [0.5, 0.75, 0], [0, 0.5, 0.25], [0, 0.5, 0.75]]
        coords_16i = [[0.185, 0.185, 0.185], [0.815, 0.815, 0.185], # etc... simplified
                      [0.185, 0.815, 0.815], [0.815, 0.185, 0.815]]
        # Simplified - only partial
        
        species = ['Si'] * (len(coords_6c) + len(coords_16i))
        coords = coords_6c + coords_16i
        
        structure = Structure(lattice, species, coords)
        
        return {
            "success": True,
            "cage_type": "clathrate_I",
            "framework": "Si",
            "structure": structure_to_dict(structure)
        }
    
    else:
        return {
            "success": False,
            "error": {"code": "INVALID_CAGE_TYPE", "message": f"Unknown cage type '{cage_type}'"}
        }


# ============================================================================
# MOF STRUCTURES
# ============================================================================

# Common MOF building blocks (simplified)
MOF_TEMPLATES = {
    "MOF-5": {
        "metal_node": "Zn4O",
        "linker": "BDC",  # benzene-1,4-dicarboxylate
        "topology": "pcu",
        "a": 25.83,
        "space_group": 225
    },
    "HKUST-1": {
        "metal_node": "Cu2",
        "linker": "BTC",  # benzene-1,3,5-tricarboxylate
        "topology": "tbo",
        "a": 26.34,
        "space_group": 225
    },
    "UiO-66": {
        "metal_node": "Zr6O4(OH)4",
        "linker": "BDC",
        "topology": "fcu",
        "a": 20.75,
        "space_group": 225
    },
    "ZIF-8": {
        "metal_node": "Zn",
        "linker": "mIm",  # 2-methylimidazolate
        "topology": "SOD",
        "a": 16.99,
        "space_group": 217
    }
}


def generate_mof_structure(
    mof_type: str,
    functionalization: Optional[str] = None,
    size: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate Metal-Organic Framework structure.
    
    Args:
        mof_type: MOF name ('MOF-5', 'HKUST-1', 'UiO-66', 'ZIF-8')
        functionalization: Optional linker functionalization
        size: Supercell size
    
    Returns:
        MOF structure dictionary
    """
    if mof_type not in MOF_TEMPLATES:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MOF",
                "message": f"Unknown MOF type '{mof_type}'",
                "available": list(MOF_TEMPLATES.keys())
            }
        }
    
    template = MOF_TEMPLATES[mof_type]
    a = template["a"]
    
    # Simplified MOF generation - create framework with approximate positions
    # Real MOFs have complex geometries; this is a simplification
    
    lattice = Lattice.cubic(a)
    
    if mof_type == "MOF-5":
        # Zn4O nodes at corners and face centers
        # BDC linkers along edges
        species = []
        coords = []
        
        # Add Zn4O clusters at 8 positions
        for node_pos in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
            # Each Zn4O cluster
            for dx, dy, dz in [[0.02, 0.02, 0.02], [-0.02, -0.02, 0.02], 
                               [0.02, -0.02, -0.02], [-0.02, 0.02, -0.02]]:
                species.append('Zn')
                coords.append([node_pos[0] + dx, node_pos[1] + dy, node_pos[2] + dz])
            species.append('O')
            coords.append(node_pos)
        
        # Add simplified BDC linkers (C6H4 + 2 COO)
        # Just add C atoms at edge centers for now
        for linker_pos in [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]]:
            for i in range(6):
                species.append('C')
                offset = [0.02 * np.sin(i * np.pi / 3), 0.02 * np.cos(i * np.pi / 3), 0]
                coords.append([linker_pos[0] + offset[0], linker_pos[1] + offset[1], linker_pos[2] + offset[2]])
    
    elif mof_type == "ZIF-8":
        # Simplified sodalite cage topology
        species = ['Zn', 'Zn', 'N', 'N', 'C', 'C', 'C']
        coords = [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0.1, 0.1, 0.1], [0.4, 0.4, 0.4],
            [0.2, 0.2, 0.2], [0.3, 0.3, 0.3], [0.25, 0.25, 0.25]
        ]
    
    else:
        # Fallback: simple placeholder
        species = ['Zn', 'O', 'C', 'C']
        coords = [[0, 0, 0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.3, 0.3, 0.3]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    return {
        "success": True,
        "mof_type": mof_type,
        "metal_node": template["metal_node"],
        "linker": template["linker"],
        "topology": template["topology"],
        "space_group": template["space_group"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure),
        "warnings": ["This is a simplified representation. For accurate MOF structures, use CIF from databases."]
    }


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """Convert Pymatgen Structure to dictionary."""
    lattice = structure.lattice
    return {
        "lattice": {
            "a": lattice.a,
            "b": lattice.b,
            "c": lattice.c,
            "alpha": lattice.alpha,
            "beta": lattice.beta,
            "gamma": lattice.gamma,
            "matrix": lattice.matrix.tolist(),
            "volume": lattice.volume
        },
        "atoms": [
            {
                "element": str(site.specie),
                "coords": list(site.frac_coords),
                "cartesian": list(site.coords)
            }
            for site in structure
        ],
        "metadata": {
            "formula": structure.formula,
            "n_atoms": len(structure),
            "density": structure.density
        }
    }


def atoms_to_dict(atoms: Atoms) -> Dict[str, Any]:
    """Convert ASE Atoms to dictionary."""
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Calculate lattice parameters
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    alpha = np.degrees(np.arccos(np.dot(cell[1], cell[2]) / (b * c))) if b * c > 0 else 90
    beta = np.degrees(np.arccos(np.dot(cell[0], cell[2]) / (a * c))) if a * c > 0 else 90
    gamma = np.degrees(np.arccos(np.dot(cell[0], cell[1]) / (a * b))) if a * b > 0 else 90
    
    volume = abs(np.linalg.det(cell))
    
    # Fractional coordinates
    inv_cell = np.linalg.inv(cell) if volume > 1e-6 else np.eye(3)
    frac_coords = positions @ inv_cell
    
    return {
        "lattice": {
            "a": float(a),
            "b": float(b),
            "c": float(c),
            "alpha": float(alpha),
            "beta": float(beta),
            "gamma": float(gamma),
            "matrix": cell.tolist(),
            "volume": float(volume)
        },
        "atoms": [
            {
                "element": sym,
                "coords": list(frac_coords[i]),
                "cartesian": list(positions[i])
            }
            for i, sym in enumerate(symbols)
        ],
        "metadata": {
            "formula": atoms.get_chemical_formula(),
            "n_atoms": len(atoms)
        }
    }


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """Main entry point for CLI usage."""
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "No input file provided"}))
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        params = json.load(f)
    
    operation = params.get("operation", "")
    
    if operation == "prototype":
        result = generate_prototype_structure(
            prototype=params.get("prototype"),
            elements=params.get("elements", {}),
            lattice_constant=float(params.get("lattice_constant", 4.0)),
            c_over_a=float(params.get("c_over_a", 1.0))
        )
    
    elif operation == "twisted_bilayer":
        result = generate_twisted_bilayer(
            material=params.get("material", "graphene"),
            twist_angle=float(params.get("twist_angle", 1.1)),
            layers=int(params.get("layers", 2)),
            stacking=params.get("stacking", "AB"),
            interlayer_distance=float(params.get("interlayer_distance", 3.35)),
            vacuum=float(params.get("vacuum", 15.0))
        )
    
    elif operation == "high_entropy_alloy":
        result = generate_high_entropy_alloy(
            elements=params.get("elements", []),
            concentrations=params.get("concentrations"),
            structure_type=params.get("structure_type", "fcc"),
            supercell=params.get("supercell", [3, 3, 3]),
            lattice_constant=params.get("lattice_constant"),
            seed=params.get("seed")
        )
    
    elif operation == "2d_material":
        result = generate_2d_material(
            material=params.get("material", "graphene"),
            size=params.get("size", [1, 1, 1]),
            vacuum=float(params.get("vacuum", 15.0)),
            **params.get("extra_params", {})
        )
    
    elif operation == "mof":
        result = generate_mof_structure(
            mof_type=params.get("mof_type", "MOF-5"),
            functionalization=params.get("functionalization"),
            size=params.get("size", [1, 1, 1])
        )
    
    elif operation == "cage":
        result = generate_cage_structure(
            cage_type=params.get("cage_type", "C60"),
            guest=params.get("guest"),
            **params.get("extra_params", {})
        )
    
    else:
        result = {
            "success": False,
            "error": {
                "code": "INVALID_OPERATION",
                "message": f"Unknown operation: {operation}",
                "available": ["prototype", "twisted_bilayer", "high_entropy_alloy", "2d_material", "mof", "cage"]
            }
        }
    
    print(json.dumps(result))


if __name__ == "__main__":
    main()
