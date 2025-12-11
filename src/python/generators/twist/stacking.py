"""
twist/stacking.py - 3D Bulk Stacking Structures

Generates bulk stacked layered materials:
- ABA (Bernal) graphite
- ABC (rhombohedral) graphite
- AAA stacking
- Turbostratic stacking
- Custom stacking sequences

Scientific basis:
- Graphite polytypes and stability
- Stacking-dependent band structure
- h-BN polymorphs
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Stacking sequences for different materials
STACKING_SEQUENCES = {
    "ABA": {
        "sequence": [0, 1, 0],  # 0=A, 1=B
        "description": "Bernal (2H) stacking - most stable graphite",
        "periodicity": 2
    },
    "ABC": {
        "sequence": [0, 1, 2],  # A, B, C
        "description": "Rhombohedral (3R) stacking",
        "periodicity": 3
    },
    "AAA": {
        "sequence": [0, 0, 0],
        "description": "AA stacking - metastable",
        "periodicity": 1
    },
    "AAB": {
        "sequence": [0, 0, 1],
        "description": "Intermediate stacking",
        "periodicity": 3
    },
    "ABAB": {
        "sequence": [0, 1, 0, 1],
        "description": "4-layer Bernal",
        "periodicity": 2
    },
    "ABCABC": {
        "sequence": [0, 1, 2, 0, 1, 2],
        "description": "6-layer rhombohedral",
        "periodicity": 3
    },
    "ABBA": {
        "sequence": [0, 1, 1, 0],
        "description": "Mirror symmetric 4-layer",
        "periodicity": 4
    }
}

# Stacking vectors
STACKING_VECTORS = {
    0: [0, 0],      # A position
    1: [1/3, 0],    # B position (shifted)
    2: [2/3, 0],    # C position (shifted)
}


def generate_stacked_layers(
    material: str = "graphene",
    stacking: str = "ABA",
    n_layers: int = 4,
    vacuum: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate bulk stacked layered structure.
    
    Args:
        material: Material type (graphene, hBN, MoS2)
        stacking: Stacking sequence (ABA, ABC, AAA, or custom)
        n_layers: Number of layers
        vacuum: Vacuum padding (None for bulk)
    
    Returns:
        Stacked structure
    """
    if material.lower() == "graphene":
        a = 2.46
        element = "C"
        z_sep = 3.35
    elif material.lower() == "hbn":
        a = 2.50
        element = ["B", "N"]
        z_sep = 3.33
    else:
        a = 3.16
        element = "Mo"
        z_sep = 6.5
    
    # Get stacking sequence
    if stacking in STACKING_SEQUENCES:
        seq_info = STACKING_SEQUENCES[stacking]
        base_seq = seq_info["sequence"]
    else:
        # Custom stacking string (e.g., "ABAC")
        char_map = {"A": 0, "B": 1, "C": 2}
        base_seq = [char_map.get(c.upper(), 0) for c in stacking]
    
    # Extend sequence to n_layers
    full_seq = []
    for i in range(n_layers):
        full_seq.append(base_seq[i % len(base_seq)])
    
    # Create structure
    if vacuum:
        c = (n_layers - 1) * z_sep + vacuum
        is_bulk = False
    else:
        c = n_layers * z_sep
        is_bulk = True
    
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    for layer_idx, stack_pos in enumerate(full_seq):
        shift = STACKING_VECTORS[stack_pos]
        z_layer = (layer_idx * z_sep + (vacuum/2 if vacuum else 0)) / c
        
        # Two atoms per layer for graphene/hBN
        x1 = (1/3 + shift[0]) % 1
        y1 = (2/3 + shift[1]) % 1
        x2 = (2/3 + shift[0]) % 1
        y2 = (1/3 + shift[1]) % 1
        
        if isinstance(element, list):
            species.extend(element)
        else:
            species.extend([element, element])
        coords.extend([[x1, y1, z_layer], [x2, y2, z_layer]])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "stacking": stacking,
        "stacking_sequence": "".join(["ABC"[s] for s in full_seq]),
        "n_layers": n_layers,
        "is_bulk": is_bulk,
        "interlayer_distance_angstrom": z_sep,
        "structure": structure_to_dict(structure, vacuum or 0)
    }


def generate_turbostratic_graphite(
    n_layers: int = 6,
    disorder_type: str = "rotational",
    disorder_amount: float = 5.0,
    vacuum: Optional[float] = None,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate turbostratic (disordered stacking) graphite.
    
    Args:
        n_layers: Number of layers
        disorder_type: 'rotational' or 'translational'
        disorder_amount: Rotation in degrees or translation in Angstrom
        vacuum: Vacuum padding
        seed: Random seed
    
    Returns:
        Turbostratic graphite structure
    """
    np.random.seed(seed)
    
    a = 2.46
    z_sep = 3.40  # Slightly larger for turbostratic
    
    if vacuum:
        c = (n_layers - 1) * z_sep + vacuum
    else:
        c = n_layers * z_sep
    
    lattice = Lattice.hexagonal(a, c)
    
    species = []
    coords = []
    
    for layer_idx in range(n_layers):
        z_layer = (layer_idx * z_sep + (vacuum/2 if vacuum else 0)) / c
        
        if disorder_type == "rotational":
            theta = np.random.uniform(-disorder_amount, disorder_amount)
            theta_rad = np.radians(theta)
            cos_t, sin_t = np.cos(theta_rad), np.sin(theta_rad)
        else:
            theta = 0
            cos_t, sin_t = 1, 0
        
        # Base positions
        for (x0, y0) in [(1/3, 2/3), (2/3, 1/3)]:
            # Apply rotation
            x_centered = x0 - 0.5
            y_centered = y0 - 0.5
            x_rot = x_centered * cos_t - y_centered * sin_t + 0.5
            y_rot = x_centered * sin_t + y_centered * cos_t + 0.5
            
            if disorder_type == "translational":
                x_rot += np.random.uniform(-disorder_amount/a, disorder_amount/a)
                y_rot += np.random.uniform(-disorder_amount/a, disorder_amount/a)
            
            species.append("C")
            coords.append([x_rot % 1, y_rot % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": "graphite",
        "stacking": "turbostratic",
        "n_layers": n_layers,
        "disorder_type": disorder_type,
        "disorder_amount": disorder_amount,
        "structure": structure_to_dict(structure, vacuum or 0)
    }
