"""
twist/multilayers.py - Twisted Multilayer Structures

Generates twisted multilayer stacks:
- Trilayers (mirror twisted, helical twisted)
- Quad-layers
- Pentalayers
- Arbitrary N-layer stacks

Scientific basis:
- Mirror vs helical twist arrangements
- Layer-dependent electronic properties
- Flat band engineering in multilayers
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import rotation_matrix_2d, structure_to_dict


def generate_twisted_trilayer(
    material: str = "graphene",
    theta12: float = 1.5,
    theta23: float = -1.5,
    arrangement: str = "mirror",
    supercell_size: int = 5,
    vacuum: float = 25.0
) -> Dict[str, Any]:
    """
    Generate twisted trilayer structure.
    
    Args:
        material: Material type
        theta12: Twist angle between layers 1 and 2
        theta23: Twist angle between layers 2 and 3
        arrangement: 'mirror' (θ, 0, -θ) or 'helical' (0, θ, 2θ)
        supercell_size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Twisted trilayer structure
    """
    if material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep = 3.35
    else:
        a = 3.16
        element = "Mo"
        z_sep = 6.5
    
    if arrangement == "mirror":
        angles = [theta12, 0, -theta12]
    elif arrangement == "helical":
        angles = [0, theta23, 2 * theta23]
    else:
        angles = [0, theta12, theta12 + theta23]
    
    supercell_size = min(supercell_size, 20)
    c = 2 * z_sep + vacuum
    lattice = Lattice.hexagonal(a * supercell_size, c)
    
    species = []
    coords = []
    
    for layer_idx in range(3):
        z_layer = (layer_idx * z_sep + vacuum/2) / c
        rot = rotation_matrix_2d(angles[layer_idx])
        
        for i in range(supercell_size):
            for j in range(supercell_size):
                for (dx, dy) in [(1/3, 2/3), (2/3, 1/3)]:
                    x = (dx + i) / supercell_size
                    y = (dy + j) / supercell_size
                    
                    pos = np.array([x - 0.5, y - 0.5])
                    pos_rot = rot @ pos + 0.5
                    
                    species.append(element)
                    coords.append([pos_rot[0] % 1, pos_rot[1] % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "n_layers": 3,
        "arrangement": arrangement,
        "twist_angles_deg": angles,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_twisted_multilayer(
    material: str = "graphene",
    n_layers: int = 4,
    twist_angles: Optional[List[float]] = None,
    twist_pattern: str = "alternating",
    base_angle: float = 1.5,
    supercell_size: int = 5,
    vacuum: float = 30.0
) -> Dict[str, Any]:
    """
    Generate twisted N-layer structure.
    
    Args:
        material: Material type
        n_layers: Number of layers (2-10)
        twist_angles: Explicit list of angles for each layer
        twist_pattern: 'alternating', 'helical', 'random'
        base_angle: Base twist angle for patterns
        supercell_size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Twisted multilayer structure
    """
    if n_layers < 2 or n_layers > 10:
        return {"success": False, "error": {"code": "INVALID_N", "message": "n_layers must be 2-10"}}
    
    if material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep = 3.35
    else:
        a = 3.16
        element = "Mo"
        z_sep = 6.5
    
    # Generate twist angles
    if twist_angles:
        angles = twist_angles[:n_layers]
        while len(angles) < n_layers:
            angles.append(0)
    elif twist_pattern == "alternating":
        angles = [base_angle if i % 2 == 0 else -base_angle for i in range(n_layers)]
    elif twist_pattern == "helical":
        angles = [base_angle * i for i in range(n_layers)]
    elif twist_pattern == "random":
        np.random.seed(42)
        angles = list(np.random.uniform(-base_angle, base_angle, n_layers))
    else:
        angles = [0] * n_layers
    
    supercell_size = min(supercell_size, 15)
    c = (n_layers - 1) * z_sep + vacuum
    lattice = Lattice.hexagonal(a * supercell_size, c)
    
    species = []
    coords = []
    
    for layer_idx in range(n_layers):
        z_layer = (layer_idx * z_sep + vacuum/2) / c
        rot = rotation_matrix_2d(angles[layer_idx])
        
        for i in range(supercell_size):
            for j in range(supercell_size):
                for (dx, dy) in [(1/3, 2/3), (2/3, 1/3)]:
                    x = (dx + i) / supercell_size
                    y = (dy + j) / supercell_size
                    
                    pos = np.array([x - 0.5, y - 0.5])
                    pos_rot = rot @ pos + 0.5
                    
                    species.append(element)
                    coords.append([pos_rot[0] % 1, pos_rot[1] % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "n_layers": n_layers,
        "twist_pattern": twist_pattern if not twist_angles else "custom",
        "twist_angles_deg": [round(a, 3) for a in angles],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure, vacuum)
    }
