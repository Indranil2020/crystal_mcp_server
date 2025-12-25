"""
two_d/base.py - Common utilities for 2D structure generation

Provides base classes and helper functions used across all 2D submodules.
Uses ASE, pymatgen, and specialized 2D libraries.
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from dataclasses import dataclass, field

from pymatgen.core import Structure, Lattice
from ase import Atoms


@dataclass
class TwoD_Structure:
    """
    Represents a 2D material structure with full metadata.
    
    Attributes:
        lattice: Lattice parameters (includes vacuum)
        atoms: List of atomic sites
        vacuum: Vacuum thickness in Å
        thickness: Layer thickness in Å
        metadata: Additional properties
    """
    lattice: Dict[str, Any]
    atoms: List[Dict[str, Any]]
    vacuum: float = 15.0
    thickness: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format."""
        return {
            "lattice": self.lattice,
            "atoms": self.atoms,
            "vacuum": self.vacuum,
            "thickness": self.thickness,
            "metadata": self.metadata
        }
    
    def to_pymatgen(self) -> Structure:
        """Convert to pymatgen Structure."""
        lattice = Lattice(self.lattice.get("matrix", [[1,0,0],[0,1,0],[0,0,1]]))
        species = [a["element"] for a in self.atoms]
        coords = [a["coords"] for a in self.atoms]
        return Structure(lattice, species, coords)


def structure_to_dict(structure: Structure, vacuum: float = 15.0) -> Dict[str, Any]:
    """
    Convert pymatgen Structure to dictionary format for 2D materials.
    
    Args:
        structure: Pymatgen Structure object
        vacuum: Vacuum thickness
    
    Returns:
        Dictionary with lattice, atoms, and metadata
    """
    lattice = structure.lattice
    
    # Calculate layer thickness from z-coordinates
    z_coords = [site.frac_coords[2] * lattice.c for site in structure]
    thickness = max(z_coords) - min(z_coords) if z_coords else 0.0
    
    atoms_list = [
        {
            "element": str(site.specie),
            "coords": list(site.frac_coords),
            "cartesian": list(site.coords)
        }
        for site in structure
    ]
    
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
        "atoms": atoms_list,
        "sites": atoms_list,  # Added for test compatibility
        "vacuum": vacuum,
        "thickness": thickness,
        "metadata": {
            "formula": structure.formula,
            "n_atoms": len(structure),
            "is_2d": True
        }
    }


def atoms_to_dict(atoms: Atoms, vacuum: float = 15.0) -> Dict[str, Any]:
    """
    Convert ASE Atoms to dictionary format for 2D materials.
    """
    cell = np.array(atoms.get_cell())
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    # Calculate angles with safety
    alpha = beta = gamma = 90.0
    if b * c > 1e-6:
        alpha = np.degrees(np.arccos(np.clip(np.dot(cell[1], cell[2]) / (b * c), -1, 1)))
    if a * c > 1e-6:
        beta = np.degrees(np.arccos(np.clip(np.dot(cell[0], cell[2]) / (a * c), -1, 1)))
    if a * b > 1e-6:
        gamma = np.degrees(np.arccos(np.clip(np.dot(cell[0], cell[1]) / (a * b), -1, 1)))
    
    volume = abs(np.linalg.det(cell))
    inv_cell = np.linalg.inv(cell) if volume > 1e-6 else np.eye(3)
    frac_coords = positions @ inv_cell
    
    # Layer thickness
    z_pos = positions[:, 2]
    thickness = max(z_pos) - min(z_pos) if len(z_pos) > 0 else 0.0
    
    atoms_list = [
        {"element": sym, "coords": list(frac_coords[i]), "cartesian": list(positions[i])}
        for i, sym in enumerate(symbols)
    ]
    
    return {
        "lattice": {
            "a": float(a), "b": float(b), "c": float(c),
            "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
            "matrix": cell.tolist(), "volume": float(volume)
        },
        "atoms": atoms_list,
        "sites": atoms_list,  # Added for test compatibility
        "vacuum": vacuum,
        "thickness": thickness,
        "metadata": {"formula": atoms.get_chemical_formula(), "n_atoms": len(atoms), "is_2d": True}
    }


def recompute_cartesian_coords(structure_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recompute Cartesian coordinates from fractional coordinates.

    This function ensures Cartesian coordinates are always consistent with
    fractional coordinates and the current lattice. Should be called whenever
    fractional coordinates or lattice parameters are modified.

    Args:
        structure_dict: Structure with 'atoms' and 'lattice' fields

    Returns:
        Structure with updated Cartesian coordinates
    """
    lattice = structure_dict.get("lattice", {})
    lattice_matrix = np.array(lattice.get("matrix", [[1,0,0],[0,1,0],[0,0,1]]))
    atoms = structure_dict.get("atoms", [])

    updated_atoms = []
    for atom in atoms:
        atom_copy = dict(atom)
        frac_coords = np.array(atom["coords"])
        # Cartesian = fractional @ lattice_matrix
        cart_coords = frac_coords @ lattice_matrix
        atom_copy["cartesian"] = cart_coords.tolist()
        updated_atoms.append(atom_copy)

    result = dict(structure_dict)
    result["atoms"] = updated_atoms
    return result


def add_vacuum(structure_dict: Dict[str, Any], vacuum: float) -> Dict[str, Any]:
    """
    Add or modify vacuum in a 2D structure.

    IMPORTANT SEMANTICS:
    - `thickness` = vertical extent of atomic layer (max_z - min_z in Cartesian coords)
    - `vacuum` = desired vacuum spacing
    - `c` (lattice parameter) = thickness + vacuum (total cell height)

    This function:
    1. Calculates new c = thickness + vacuum
    2. Scales the lattice matrix to match new c
    3. Recomputes Cartesian coordinates to reflect the new lattice

    Args:
        structure_dict: 2D structure
        vacuum: New vacuum thickness in Angstroms

    Returns:
        Structure with updated vacuum and consistent Cartesian coordinates
    """
    lattice = structure_dict.get("lattice", {})
    current_c = lattice.get("c", 20.0)
    thickness = structure_dict.get("thickness", 0.0)

    # SEMANTICS: c = thickness + vacuum (total cell height)
    new_c = thickness + vacuum
    scale = new_c / current_c if current_c > 1e-6 else 1.0

    # Update lattice
    new_lattice = dict(lattice)
    new_lattice["c"] = new_c

    if "matrix" in lattice:
        matrix = np.array(lattice["matrix"])
        # Scale only the c-axis (third lattice vector)
        matrix[2] = matrix[2] * scale
        new_lattice["matrix"] = matrix.tolist()
        new_lattice["volume"] = float(abs(np.linalg.det(matrix)))

    result = dict(structure_dict)
    result["lattice"] = new_lattice
    result["vacuum"] = vacuum

    # CRITICAL: Recompute Cartesian coordinates with new lattice
    result = recompute_cartesian_coords(result)

    return result


def center_in_vacuum(structure_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Center 2D layer in the middle of the vacuum.

    Shifts fractional coordinates so the layer is centered at z=0.5 in fractional
    coordinates, then recomputes Cartesian coordinates to maintain consistency.

    Args:
        structure_dict: 2D structure

    Returns:
        Structure with layer centered at z=0.5 and updated Cartesian coordinates
    """
    atoms = structure_dict.get("atoms", [])

    if not atoms:
        return structure_dict

    # Find z-center in fractional coordinates
    z_coords = [a["coords"][2] for a in atoms]
    z_center = (max(z_coords) + min(z_coords)) / 2
    shift = 0.5 - z_center

    # Shift fractional coordinates
    new_atoms = []
    for atom in atoms:
        atom_copy = dict(atom)
        coords = list(atom["coords"])
        coords[2] = coords[2] + shift
        atom_copy["coords"] = coords
        new_atoms.append(atom_copy)

    result = dict(structure_dict)
    result["atoms"] = new_atoms

    # CRITICAL: Recompute Cartesian coordinates after modifying fractional coords
    result = recompute_cartesian_coords(result)

    return result
