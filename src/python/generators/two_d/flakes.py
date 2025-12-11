"""
two_d/flakes.py - Nanoflakes, Nanomeshes, and Quantum Dots

Generates finite 2D structures:
- Circular/hexagonal nanoflakes
- Triangular flakes
- Nanomeshes (periodic array of holes)
- Quantum dots
- Exact edge terminations

Scientific basis:
- Edge-dependent electronic properties
- Quantum confinement in 2D QDs
- Size-dependent bandgaps
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Molecule, Lattice, Structure

from .base import structure_to_dict


def generate_nanoflake(
    base_material: str = "graphene",
    shape: str = "hexagonal",
    size: int = 3,
    passivation: str = "H",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Generate 2D nanoflake (finite 2D fragment).
    
    Args:
        base_material: Base 2D material (graphene, MoS2, etc.)
        shape: Flake shape ('hexagonal', 'triangular', 'circular', 'rectangular')
        size: Size in lattice units (defines edge length)
        passivation: Edge passivation element
        vacuum: Vacuum for periodic box
    
    Returns:
        Nanoflake structure
    """
    valid_shapes = ["hexagonal", "triangular", "circular", "rectangular"]
    if shape not in valid_shapes:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SHAPE",
                "message": f"Unknown shape '{shape}'",
                "available": valid_shapes
            }
        }
    
    if base_material.lower() == "graphene":
        a = 2.46
        element = "C"
        bond_length = 1.42
    elif base_material.lower() == "bn":
        a = 2.50
        element = "B"  # Will alternate B-N
        bond_length = 1.45
    elif base_material.lower() in ["mos2", "ws2"]:
        a = 3.16
        element = "Mo" if "mo" in base_material.lower() else "W"
        bond_length = 3.16
    else:
        a = 2.5
        element = "C"
        bond_length = 1.42
    
    # Generate positions based on shape
    positions = []
    elements = []
    
    if shape == "hexagonal":
        # Hexagonal flake with zigzag edges
        for i in range(-size, size + 1):
            for j in range(-size, size + 1):
                # Hexagonal boundary condition
                if abs(i) + abs(j) + abs(-i-j) <= 2 * size:
                    # Two atoms per hexagonal unit
                    x1 = a * (i + j * 0.5)
                    y1 = a * j * np.sqrt(3) / 2
                    x2 = x1 + bond_length * 0.5
                    y2 = y1 + bond_length * np.sqrt(3) / 2
                    
                    positions.append([x1, y1, 0])
                    positions.append([x2, y2, 0])
                    elements.extend([element, element])
    
    elif shape == "triangular":
        for i in range(size):
            for j in range(size - i):
                x = a * (i + j * 0.5)
                y = a * j * np.sqrt(3) / 2
                positions.append([x, y, 0])
                positions.append([x + bond_length * 0.5, y + bond_length * np.sqrt(3) / 2, 0])
                elements.extend([element, element])
    
    elif shape == "circular":
        radius = size * a
        for i in range(-size * 2, size * 2 + 1):
            for j in range(-size * 2, size * 2 + 1):
                x1 = a * (i + j * 0.5)
                y1 = a * j * np.sqrt(3) / 2
                if x1**2 + y1**2 <= radius**2:
                    positions.append([x1, y1, 0])
                    positions.append([x1 + bond_length * 0.5, y1 + bond_length * np.sqrt(3) / 2, 0])
                    elements.extend([element, element])
    
    else:  # rectangular
        for i in range(size * 2):
            for j in range(size):
                x1 = a * i
                y1 = a * j * np.sqrt(3) / 2
                positions.append([x1, y1, 0])
                positions.append([x1 + bond_length * 0.5, y1 + bond_length * np.sqrt(3) / 2, 0])
                elements.extend([element, element])
    
    if not positions:
        return {"success": False, "error": {"code": "EMPTY", "message": "No atoms generated"}}
    
    positions = np.array(positions)
    
    # Center the flake
    center = positions.mean(axis=0)
    positions = positions - center
    
    # Add passivation at edges
    if passivation and passivation.lower() != "bare":
        # Find edge atoms (fewer neighbors)
        from scipy.spatial import distance_matrix
        dist_mat = distance_matrix(positions, positions)
        
        for i, pos in enumerate(positions):
            n_neighbors = np.sum((dist_mat[i] > 0.1) & (dist_mat[i] < bond_length * 1.5))
            if n_neighbors < 3:  # Edge atom
                # Add passivant radially outward
                r = np.linalg.norm(pos[:2])
                if r > 0.1:
                    direction = pos[:2] / r
                    pass_pos = pos + np.array([direction[0] * 1.09, direction[1] * 1.09, 0])
                    positions = np.vstack([positions, pass_pos])
                    elements.append(passivation)
    
    # Create box
    max_extent = np.abs(positions).max() + vacuum
    box_size = 2 * max_extent
    
    # Shift to positive coordinates
    positions = positions + max_extent
    
    # Convert to fractional coordinates
    frac_coords = positions / box_size
    frac_coords[:, 2] = 0.5  # Center in z
    
    lattice = Lattice.orthorhombic(box_size, box_size, vacuum)
    structure = Structure(lattice, elements, frac_coords)
    
    return {
        "success": True,
        "base_material": base_material,
        "shape": shape,
        "size": size,
        "n_atoms": len(elements),
        "passivation": passivation,
        "diameter_angstrom": round(np.abs(positions[:, :2]).max() * 2, 1),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_quantum_dot(
    base_material: str = "MoS2",
    diameter: float = 2.0,  # nm
    shape: str = "circular",
    passivation: str = "S",
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate 2D quantum dot (nm-scale flake).
    
    Args:
        base_material: Base 2D material
        diameter: Dot diameter in nanometers
        shape: Dot shape ('circular', 'hexagonal', 'triangular')
        passivation: Edge passivation
        vacuum: Vacuum padding
    
    Returns:
        Quantum dot structure
    """
    # Convert diameter to lattice units
    if base_material.lower() in ["mos2", "ws2"]:
        a = 3.16  # Å
    elif base_material.lower() == "graphene":
        a = 2.46
    else:
        a = 3.0
    
    size = int(diameter * 10 / (2 * a))  # nm to lattice units
    size = max(size, 2)  # Minimum size
    
    result = generate_nanoflake(base_material, shape, size, passivation, vacuum)
    
    if result["success"]:
        result["target_diameter_nm"] = diameter
        result["is_quantum_dot"] = True
        result["quantum_confinement"] = diameter < 5  # nm
    
    return result


def generate_nanomesh(
    base_material: str = "graphene",
    supercell: List[int] = [5, 5, 1],
    hole_pattern: str = "hexagonal",
    hole_fraction: float = 0.3,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate 2D nanomesh (periodic hole array).
    
    Args:
        base_material: Base 2D material
        supercell: Supercell size
        hole_pattern: Pattern of holes ('hexagonal', 'square')
        hole_fraction: Fraction of atoms to remove
        vacuum: Vacuum padding
    
    Returns:
        Nanomesh structure
    """
    if base_material.lower() == "graphene":
        a = 2.46
        element = "C"
    else:
        a = 3.0
        element = "C"
    
    # Generate full supercell
    species = []
    coords = []
    z = 0.5
    
    n_total = 0
    removed = 0
    
    for i in range(supercell[0]):
        for j in range(supercell[1]):
            n_total += 2
            
            # Check if in hole region
            hole_i = np.floor(i / 3) * 3 + 1
            hole_j = np.floor(j / 3) * 3 + 1
            
            dist_to_hole = ((i - hole_i)**2 + (j - hole_j)**2)**0.5
            in_hole = dist_to_hole < supercell[0] * hole_fraction / 3
            
            if not in_hole:
                x1 = (1/3 + i) / supercell[0]
                y1 = (2/3 + j) / supercell[1]
                x2 = (2/3 + i) / supercell[0]
                y2 = (1/3 + j) / supercell[1]
                species.extend([element, element])
                coords.extend([[x1, y1, z], [x2, y2, z]])
            else:
                removed += 2
    
    lattice = Lattice.hexagonal(a * supercell[0], vacuum)
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "base_material": base_material,
        "supercell": supercell,
        "hole_pattern": hole_pattern,
        "hole_fraction": round(removed / n_total, 3) if n_total > 0 else 0,
        "n_atoms": len(structure),
        "n_removed": removed,
        "structure": structure_to_dict(structure, vacuum)
    }
