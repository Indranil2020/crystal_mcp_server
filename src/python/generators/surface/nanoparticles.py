"""
surface/nanoparticles.py - Nanoparticle Generation

Generates nanoparticle structures:
- Wulff construction (equilibrium shape)
- Core-shell nanoparticles
- Various shapes: sphere, cube, octahedron, icosahedron
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

try:
    from ase.cluster import FaceCenteredCubic, Icosahedron, Octahedron, Decahedron
    HAS_ASE_CLUSTER = True
except ImportError:
    HAS_ASE_CLUSTER = False


def atoms_to_dict(atoms) -> Dict[str, Any]:
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell()
    
    return {
        "lattice": {"matrix": cell.tolist(), "a": np.linalg.norm(cell[0])},
        "atoms": [{"element": s, "coords": list(p / 10), "cartesian": list(p)} for s, p in zip(symbols, positions)],
        "metadata": {"n_atoms": len(atoms), "formula": atoms.get_chemical_formula()}
    }


def generate_wulff_nanoparticle(
    element: str = "Au",
    size: int = 3,
    shape: str = "fcc_truncated_octahedron"
) -> Dict[str, Any]:
    """
    Generate Wulff-shaped nanoparticle.
    
    Args:
        element: Element symbol
        size: Size parameter (number of atomic layers)
        shape: Shape type
    
    Returns:
        Nanoparticle structure
    """
    lattice_constants = {"Au": 4.08, "Ag": 4.09, "Cu": 3.61, "Pt": 3.92, "Pd": 3.89, "Ni": 3.52}
    a = lattice_constants.get(element, 4.0)
    
    if HAS_ASE_CLUSTER:
        if shape == "icosahedron":
            atoms = Icosahedron(element, size)
        elif shape == "octahedron":
            atoms = Octahedron(element, size)
        elif shape == "decahedron":
            atoms = Decahedron(element, p=size, q=size, r=0)
        else:  # FCC truncated octahedron
            surfaces = [(1, 0, 0), (1, 1, 1)]
            layers = [size, size]
            atoms = FaceCenteredCubic(element, surfaces, layers, latticeconstant=a)
        
        return {
            "success": True,
            "element": element,
            "shape": shape,
            "size_parameter": size,
            "n_atoms": len(atoms),
            "structure": atoms_to_dict(atoms)
        }
    else:
        # Simple cubic cluster without ASE
        n_per_side = size * 2 + 1
        positions = []
        center = np.array([a * size, a * size, a * size])
        
        for i in range(n_per_side):
            for j in range(n_per_side):
                for k in range(n_per_side):
                    pos = np.array([i, j, k]) * a
                    if np.linalg.norm(pos - center) < a * size:
                        positions.append(pos)
        
        return {
            "success": True,
            "element": element,
            "shape": "sphere",
            "n_atoms": len(positions),
            "structure": {
                "atoms": [{"element": element, "cartesian": list(p)} for p in positions]
            }
        }


def generate_core_shell(
    core_element: str = "Au",
    shell_element: str = "Pt",
    core_size: int = 2,
    shell_thickness: int = 1,
    shape: str = "sphere"
) -> Dict[str, Any]:
    """
    Generate core-shell nanoparticle.
    
    Args:
        core_element: Core element
        shell_element: Shell element
        core_size: Core size parameter
        shell_thickness: Shell thickness in layers
        shape: Particle shape
    
    Returns:
        Core-shell nanoparticle
    """
    a_core = {"Au": 4.08, "Ag": 4.09, "Cu": 3.61}.get(core_element, 4.0)
    a_shell = {"Au": 4.08, "Ag": 4.09, "Pt": 3.92, "Pd": 3.89}.get(shell_element, 4.0)
    
    total_size = core_size + shell_thickness
    core_radius = core_size * a_core
    
    positions = []
    elements = []
    
    n_per_side = total_size * 2 + 1
    center = np.array([total_size, total_size, total_size]) * a_shell
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                pos = np.array([i, j, k]) * a_shell
                dist = np.linalg.norm(pos - center)
                
                if dist < total_size * a_shell:
                    if dist < core_radius:
                        elements.append(core_element)
                    else:
                        elements.append(shell_element)
                    positions.append(pos - center)
    
    return {
        "success": True,
        "core_element": core_element,
        "shell_element": shell_element,
        "core_size": core_size,
        "shell_thickness": shell_thickness,
        "n_atoms": len(positions),
        "n_core_atoms": sum(1 for e in elements if e == core_element),
        "n_shell_atoms": sum(1 for e in elements if e == shell_element),
        "structure": {
            "atoms": [{"element": e, "cartesian": list(p)} for e, p in zip(elements, positions)]
        }
    }
