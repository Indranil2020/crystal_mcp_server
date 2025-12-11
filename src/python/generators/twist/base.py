"""
twist/base.py - Common utilities for twisted structures

Provides base classes and helper functions for:
- Moiré angle calculations
- Commensurate cell finding
- Layer alignment
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from dataclasses import dataclass, field

from pymatgen.core import Structure, Lattice


@dataclass
class TwistStructure:
    """
    Represents a twisted multilayer structure.
    """
    lattice: Dict[str, Any]
    layers: List[Dict[str, Any]]
    twist_angles: List[float]
    interlayer_distances: List[float]
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "lattice": self.lattice,
            "layers": self.layers,
            "twist_angles": self.twist_angles,
            "interlayer_distances": self.interlayer_distances,
            "metadata": self.metadata
        }


def calculate_moire_angle(theta: float, a1: float, a2: Optional[float] = None) -> Dict[str, Any]:
    """
    Calculate Moiré superlattice periodicity.
    
    Args:
        theta: Twist angle in degrees
        a1: Lattice constant of layer 1
        a2: Lattice constant of layer 2 (if different from a1)
    
    Returns:
        Moiré parameters
    """
    if a2 is None:
        a2 = a1
    
    theta_rad = np.radians(theta)
    
    if abs(a1 - a2) < 1e-6:
        # Homo-bilayer: L_M = a / (2 * sin(theta/2))
        if abs(theta) < 0.1:
            L_M = float('inf')
        else:
            L_M = a1 / (2 * np.sin(theta_rad / 2))
    else:
        # Hetero-bilayer with lattice mismatch
        delta = (a2 - a1) / a1
        L_M = a1 / np.sqrt(delta**2 + (2 * np.sin(theta_rad / 2))**2)
    
    # Number of atoms in Moiré unit cell (approximate)
    n_atoms_approx = 2 * (L_M / a1) ** 2
    
    return {
        "twist_angle_deg": theta,
        "moire_period_angstrom": round(L_M, 2),
        "n_atoms_approx": int(n_atoms_approx),
        "lattice_mismatch_percent": round(abs(a1 - a2) / a1 * 100, 3) if a2 else 0
    }


def find_commensurate_cell(
    theta: float,
    a: float,
    max_size: int = 100
) -> Optional[Tuple[int, int, int, int]]:
    """
    Find commensurate supercell for twisted bilayer.
    
    Uses the algorithm from Moiré literature to find (m, r) such that
    the twist angle is exactly representable.
    
    Args:
        theta: Target twist angle (degrees)
        a: Lattice constant
        max_size: Maximum supercell size to search
    
    Returns:
        Tuple (m, r, nx, ny) for commensurate cell, or None
    """
    theta_rad = np.radians(theta)
    target_cos = np.cos(theta_rad)
    
    best_match = None
    best_error = float('inf')
    
    for m in range(1, max_size):
        for r in range(0, m + 1):
            # Commensurate angle formula
            cos_approx = (3*m**2 + 3*m*r + r**2/2) / (3*m**2 + 3*m*r + r**2)
            
            if abs(cos_approx) <= 1:
                error = abs(cos_approx - target_cos)
                if error < best_error:
                    best_error = error
                    nx = m
                    ny = m + r
                    best_match = (m, r, nx, ny)
                    
                    if error < 0.001:
                        return best_match
    
    return best_match


def rotation_matrix_2d(theta: float) -> np.ndarray:
    """
    2D rotation matrix.
    
    Args:
        theta: Angle in degrees
    
    Returns:
        2x2 rotation matrix
    """
    theta_rad = np.radians(theta)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s], [s, c]])


def structure_to_dict(structure: Structure, vacuum: float = 15.0) -> Dict[str, Any]:
    """Convert pymatgen Structure to dictionary format."""
    lattice = structure.lattice
    
    return {
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "matrix": lattice.matrix.tolist(), "volume": lattice.volume
        },
        "atoms": [
            {"element": str(site.specie), "coords": list(site.frac_coords), "cartesian": list(site.coords)}
            for site in structure
        ],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }
