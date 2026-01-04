"""
bulk/base.py - Common utilities for bulk structure generation

Provides base classes and helper functions used across all bulk submodules.
Uses PyXtal, ASE, pymatgen, and spglib for scientific accuracy.
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
from dataclasses import dataclass, field

# Core scientific libraries
from pymatgen.core import Structure, Lattice, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase import Atoms
import spglib


@dataclass
class BulkStructure:
    """
    Represents a crystalline bulk structure with full metadata.
    
    Attributes:
        lattice: Lattice parameters and matrix
        atoms: List of atomic sites
        space_group: Space group information
        metadata: Additional properties (formula, density, etc.)
    """
    lattice: Dict[str, Any]
    atoms: List[Dict[str, Any]]
    space_group: Optional[Dict[str, Any]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format."""
        return {
            "lattice": self.lattice,
            "atoms": self.atoms,
            "space_group": self.space_group,
            "metadata": self.metadata
        }
    
    def to_pymatgen(self) -> Structure:
        """Convert to pymatgen Structure."""
        lattice = Lattice(self.lattice.get("matrix", [[1,0,0],[0,1,0],[0,0,1]]))
        species = [a["element"] for a in self.atoms]
        coords = [a["coords"] for a in self.atoms]
        return Structure(lattice, species, coords)
    
    def to_ase(self) -> Atoms:
        """Convert to ASE Atoms."""
        pmg = self.to_pymatgen()
        from pymatgen.io.ase import AseAtomsAdaptor
        return AseAtomsAdaptor.get_atoms(pmg)


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """
    Convert pymatgen Structure to dictionary format.
    
    Args:
        structure: Pymatgen Structure object
    
    Returns:
        Dictionary with lattice, atoms, and metadata
    """
    lattice = structure.lattice
    
    # Get symmetry info
    sga = SpacegroupAnalyzer(structure, symprec=0.1)
    sg_symbol = sga.get_space_group_symbol()
    sg_number = sga.get_space_group_number()
    hall = sga.get_hall()
    point_group = sga.get_point_group_symbol()
    crystal_system = sga.get_crystal_system()
    
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
        "space_group": {
            "number": sg_number,
            "symbol": sg_symbol,
            "hall_symbol": hall,
            "point_group": point_group,
            "crystal_system": crystal_system
        },
        "metadata": {
            "formula": structure.formula,
            "reduced_formula": structure.composition.reduced_formula,
            "n_atoms": len(structure),
            "density": structure.density
        }
    }


def atoms_to_dict(atoms: Atoms) -> Dict[str, Any]:
    """
    Convert ASE Atoms to dictionary format.
    
    Args:
        atoms: ASE Atoms object
    
    Returns:
        Dictionary with lattice, atoms, and metadata
    """
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Compute lattice parameters
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    # Handle degenerate cases
    alpha = beta = gamma = 90.0
    if b * c > 1e-6:
        alpha = np.degrees(np.arccos(np.clip(np.dot(cell[1], cell[2]) / (b * c), -1, 1)))
    if a * c > 1e-6:
        beta = np.degrees(np.arccos(np.clip(np.dot(cell[0], cell[2]) / (a * c), -1, 1)))
    if a * b > 1e-6:
        gamma = np.degrees(np.arccos(np.clip(np.dot(cell[0], cell[1]) / (a * b), -1, 1)))
    
    volume = abs(np.linalg.det(cell))
    
    # Fractional coordinates
    if volume > 1e-6:
        inv_cell = np.linalg.inv(cell)
        frac_coords = positions @ inv_cell
    else:
        frac_coords = positions
    
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


def validate_composition(composition: List[str]) -> Tuple[bool, Optional[str]]:
    """
    Validate a chemical composition.
    
    Args:
        composition: List of element symbols
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not composition:
        return False, "Composition cannot be empty"
    
    for elem in composition:
        if elem not in [e.symbol for e in Element]:
            return False, f"Invalid element symbol: {elem}"
    
    return True, None


def validate_space_group(sg: Union[int, str]) -> Tuple[bool, Optional[str], Optional[int]]:
    """
    Validate and normalize a space group.
    
    Args:
        sg: Space group number (1-230) or symbol
    
    Returns:
        Tuple of (is_valid, error_message, sg_number)
    """
    if isinstance(sg, int):
        if 1 <= sg <= 230:
            return True, None, sg
        return False, f"Space group number must be 1-230, got {sg}", None
    
    if isinstance(sg, str):
        # Try to lookup by symbol
        from pymatgen.symmetry.groups import SpaceGroup
        from pymatgen.symmetry.groups import sg_symbol_from_int_number
        
        for i in range(1, 231):
            sym = sg_symbol_from_int_number(i)
            if sg.replace(" ", "") == sym.replace(" ", ""):
                return True, None, i
        
        return False, f"Unknown space group symbol: {sg}", None
    
    return False, f"Space group must be int or str, got {type(sg)}", None
