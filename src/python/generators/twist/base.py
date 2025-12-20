"""
twist/base.py - Twist Base Utilities

Common utilities and databases for twisted/Moiré structures.
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


class TwistStructure:
    """
    Class to hold twisted structure information.
    """
    def __init__(
        self,
        structure: Structure,
        twist_angle: float,
        layer_index: Optional[List[int]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        self.structure = structure
        self.twist_angle = twist_angle
        self.layer_index = layer_index or []
        self.metadata = metadata or {}
        
    def to_dict(self) -> Dict[str, Any]:
        return {
            "structure": structure_to_dict(self.structure),
            "twist_angle": self.twist_angle,
            "metadata": self.metadata
        }



# Twist angle database for common systems
MAGIC_ANGLES = {
    # Graphene
    "twisted_bilayer_graphene": {
        "magic_1": 1.08, "magic_2": 0.50, "magic_3": 0.35,
        "flat_band": True, "correlated": True, "superconducting": True
    },
    "twisted_trilayer_graphene": {
        "magic_1": 1.57, "alternating": True,
        "displacement_field": True, "FQHE": True
    },
    
    # TMDs
    "twisted_MoS2": {
        "flat_band_angle": 1.5, "moire_excitons": True,
        "interlayer_excitons": True
    },
    "twisted_WSe2": {
        "flat_band_angle": 2.0, "moire_potential_meV": 100,
        "Moiré_Hubbard": True
    },
    "MoSe2_WSe2_heterobilayer": {
        "type_II": True, "moire_trapping": True,
        "spatially_indirect_excitons": True
    },
    
    # Other 2D materials
    "twisted_hBN": {
        "ferroelectric": True, "polar_domains": True
    },
    "twisted_phosphorene": {
        "anisotropic": True, "angle_dependent_bandgap": True
    },
}


# Commensurate angle database
COMMENSURATE_ANGLES = {
    # (n, m) indices for twisted bilayer graphene
    (1, 0): 60.0,
    (2, 1): 21.79,
    (3, 2): 13.17,
    (4, 3): 9.43,
    (5, 4): 7.34,
    (6, 5): 6.01,
    (7, 6): 5.09,
    (8, 7): 4.41,
    (9, 8): 3.89,
    (10, 9): 3.48,
    (11, 10): 3.15,
    (21, 20): 1.70,
    (31, 30): 1.16,
    (32, 31): 1.12,  # Near first magic angle
}


# Stacking registry database
STACKING_REGISTRY = {
    "graphene": {
        "AA": {"energy_meV_atom": 20, "interlayer_A": 3.60},
        "AB": {"energy_meV_atom": 0, "interlayer_A": 3.35, "Bernal": True},
        "SP": {"energy_meV_atom": 10, "saddle_point": True},
    },
    "MoS2": {
        "AA": {"energy_meV_atom": 30, "interlayer_A": 6.5},
        "AB": {"energy_meV_atom": 0, "interlayer_A": 6.15, "2H": True},
        "AA'": {"energy_meV_atom": 15, "3R": True},
    },
    "hBN": {
        "AA": {"energy_meV_atom": 25, "polar": True},
        "AB": {"energy_meV_atom": 0, "centrosymmetric": True},
        "BA": {"energy_meV_atom": 0, "centrosymmetric": True},
    },
}


def calculate_moire_period(
    a: float,
    theta_deg: float
) -> float:
    """
    Calculate Moiré period from twist angle.
    
    Args:
        a: Lattice constant in Angstrom
        theta_deg: Twist angle in degrees
    
    Returns:
        Moiré period in Angstrom
    """
    theta = np.radians(theta_deg)
    if np.sin(theta / 2) < 1e-10:
        return float('inf')
    return a / (2 * np.sin(theta / 2))


def calculate_moire_atoms(
    a: float,
    theta_deg: float,
    atoms_per_unit_cell: int = 2
) -> int:
    """
    Calculate number of atoms in Moiré unit cell.
    
    Args:
        a: Lattice constant
        theta_deg: Twist angle
        atoms_per_unit_cell: Atoms in primitive cell
    
    Returns:
        Number of atoms in Moiré cell (per layer)
    """
    L = calculate_moire_period(a, theta_deg)
    if L == float('inf'):
        return 2
    
    # For graphene-like honeycomb
    area_moire = L**2 * np.sqrt(3) / 2
    area_unit = a**2 * np.sqrt(3) / 2
    
    n_unit_cells = int(area_moire / area_unit)
    return n_unit_cells * atoms_per_unit_cell


def find_commensurate(
    target_angle: float,
    max_n: int = 50
) -> Tuple[int, int, float]:
    """
    Find nearest commensurate angle.
    
    Args:
        target_angle: Target angle in degrees
        max_n: Maximum index to search
    
    Returns:
        (n, m, actual_angle)
    """
    best_n, best_m = 1, 0
    best_angle = 60.0
    min_diff = abs(target_angle - 60.0)
    
    for n in range(1, max_n + 1):
        for m in range(n):
            # Commensurate angle formula
            angle = np.degrees(np.arccos((n**2 + 4*n*m + m**2) / (2*(n**2 + n*m + m**2))))
            
            diff = abs(target_angle - angle)
            if diff < min_diff:
                min_diff = diff
                best_n, best_m = n, m
                best_angle = angle
    
    return (best_n, best_m, best_angle)


def calculate_moire_angle(n: int, m: int) -> float:
    """
    Calculate twist angle for commensurate indices (n, m).
    """
    try:
        cos_theta = (n**2 + 4*n*m + m**2) / (2*(n**2 + n*m + m**2))
        return np.degrees(np.arccos(cos_theta))
    except:
        return 0.0


def calculate_moire_wavelength(theta_deg: float, a: float = 2.46) -> float:
    """Wrapper for calculate_moire_period."""
    return calculate_moire_period(a, theta_deg)


def find_commensurate_cell(target_angle: float, max_n: int = 50) -> Tuple[int, int, float]:
    """Alias for find_commensurate."""
    return find_commensurate(target_angle, max_n)


def rotation_matrix_2d(angle_deg: float) -> np.ndarray:
    """Generate 2D rotation matrix."""
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])



def calculate_twist_strain(
    a1: float,
    a2: float,
    theta_deg: float
) -> Dict[str, float]:
    """
    Calculate strain from lattice mismatch in twisted heterostructure.
    
    Args:
        a1: Lattice constant of layer 1
        a2: Lattice constant of layer 2
        theta_deg: Twist angle
    
    Returns:
        Strain information
    """
    mismatch = abs(a1 - a2) / ((a1 + a2) / 2) * 100
    
    L1 = calculate_moire_period(a1, theta_deg)
    L2 = calculate_moire_period(a2, theta_deg)
    
    return {
        "lattice_mismatch_percent": round(mismatch, 2),
        "moire_period_layer1_A": round(L1, 1),
        "moire_period_layer2_A": round(L2, 1),
        "average_moire_period_A": round((L1 + L2) / 2, 1),
        "commensurate": mismatch < 0.5 and theta_deg in [1.08, 1.1, 21.8]
    }


def get_twist_database() -> Dict[str, Any]:
    """Get all twist system information."""
    return {
        "success": True,
        "magic_angle_systems": list(MAGIC_ANGLES.keys()),
        "commensurate_angles": {str(k): v for k, v in COMMENSURATE_ANGLES.items()},
        "stacking_registries": list(STACKING_REGISTRY.keys()),
        "first_magic_angle_graphene": 1.08,
    }


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_twisted_bilayer_simple(
    material: str = "graphene",
    twist_angle_deg: float = 1.08,
    interlayer_distance_A: float = 3.35
) -> Dict[str, Any]:
    """
    Generate simplified twisted bilayer structure (small angle approximation).
    
    Args:
        material: Base material
        twist_angle_deg: Twist angle
        interlayer_distance_A: Interlayer distance
    
    Returns:
        Twisted bilayer structure
    """
    lattice_params = {"graphene": 2.46, "MoS2": 3.16, "hBN": 2.50, "WSe2": 3.28}
    a = lattice_params.get(material, 2.46)
    
    L_moire = calculate_moire_period(a, twist_angle_deg)
    
    # For small angles, use Moiré period as lattice constant
    if L_moire > 100:
        # Too large, use smaller approximation
        L_moire = 50
    
    lattice = Lattice.hexagonal(L_moire, interlayer_distance_A * 2 + 20)
    
    species = []
    coords = []
    
    n_atoms_layer = calculate_moire_atoms(a, twist_angle_deg)
    n_atoms_layer = min(n_atoms_layer, 500)  # Limit for practical structures
    
    # Generate atoms for both layers
    theta = np.radians(twist_angle_deg / 2)
    
    for layer in range(2):
        sign = 1 if layer == 0 else -1
        z = 0.4 if layer == 0 else 0.6
        
        for i in range(int(np.sqrt(n_atoms_layer))):
            for j in range(int(np.sqrt(n_atoms_layer))):
                # Rotate positions
                x = i / np.sqrt(n_atoms_layer)
                y = j / np.sqrt(n_atoms_layer)
                
                x_rot = x * np.cos(sign * theta) - y * np.sin(sign * theta)
                y_rot = x * np.sin(sign * theta) + y * np.cos(sign * theta)
                
                species.append("C" if material == "graphene" else material[:2])
                coords.append([x_rot % 1, y_rot % 1, z])
    
    structure = Structure(lattice, species, coords)
    
    info = MAGIC_ANGLES.get(f"twisted_bilayer_{material}", {})
    
    return {
        "success": True,
        "material": material,
        "twist_angle_deg": twist_angle_deg,
        "moire_period_A": round(L_moire, 1),
        "interlayer_distance_A": interlayer_distance_A,
        "n_atoms": len(structure),
        "is_magic_angle": abs(twist_angle_deg - 1.08) < 0.1,
        "flat_band_expected": info.get("flat_band", False),
        "structure": structure_to_dict(structure)
    }
