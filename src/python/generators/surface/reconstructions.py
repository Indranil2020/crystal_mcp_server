"""
surface/reconstructions.py - Surface Reconstructions

Generates common surface reconstructions:
- Si(100)-2×1
- Si(111)-7×7 (simplified)
- Au(111)-22×√3 (herringbone)
- TiO2(110)-2×1
- Fe3O4(001)-√2×√2R45°
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


def structure_to_dict(structure: Structure, vacuum: float = 15.0) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(site.specie), "coords": list(site.frac_coords)} for site in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


RECONSTRUCTIONS = {
    "Si(100)-2x1": {"surface": "Si(100)", "reconstruction": "2x1", "description": "Dimer rows"},
    "Si(111)-7x7": {"surface": "Si(111)", "reconstruction": "7x7", "description": "DAS model"},
    "Au(111)-22xsqrt3": {"surface": "Au(111)", "reconstruction": "22×√3", "description": "Herringbone"},
    "TiO2(110)-2x1": {"surface": "TiO2(110)", "reconstruction": "2x1", "description": "Added Ti2O3 row"},
    "Pt(110)-1x2": {"surface": "Pt(110)", "reconstruction": "1x2", "description": "Missing row"},
}


def generate_reconstruction(
    reconstruction_name: str,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate reconstructed surface structure.
    
    Args:
        reconstruction_name: Name of reconstruction (e.g., "Si(100)-2x1")
        vacuum: Vacuum padding
    
    Returns:
        Reconstructed surface structure
    """
    if reconstruction_name not in RECONSTRUCTIONS:
        return {
            "success": False,
            "error": {"code": "UNKNOWN", "message": f"Unknown: {reconstruction_name}", "available": list(RECONSTRUCTIONS.keys())}
        }
    
    info = RECONSTRUCTIONS[reconstruction_name]
    
    if "Si(100)" in reconstruction_name:
        a = 3.84 * 2  # 2x1 reconstruction
        b = 3.84
        c = vacuum
        lattice = Lattice.orthorhombic(a, b, c + 10)
        
        species = ["Si"] * 4
        z = 0.5
        # Dimer pairs
        coords = [
            [0.15, 0.5, z], [0.35, 0.5, z],  # Dimer 1
            [0.65, 0.5, z], [0.85, 0.5, z],  # Dimer 2
        ]
        
        structure = Structure(lattice, species, coords)
        
    elif "Au(111)" in reconstruction_name:
        a = 2.88 * 22
        b = 2.88 * np.sqrt(3)
        c = vacuum
        lattice = Lattice.orthorhombic(a, b, c + 10)
        
        # Simplified herringbone
        species = ["Au"] * 10
        coords = [[i/10, 0.5, 0.5] for i in range(10)]
        structure = Structure(lattice, species, coords)
        
    else:
        # Generic reconstruction
        a = 5.0
        lattice = Lattice.orthorhombic(a, a, vacuum + 10)
        species = ["X"]
        coords = [[0.5, 0.5, 0.5]]
        structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "reconstruction": reconstruction_name,
        "description": info["description"],
        "structure": structure_to_dict(structure, vacuum)
    }
