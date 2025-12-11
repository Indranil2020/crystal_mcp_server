"""
molecule/small_molecules.py - Small Molecule Generation

Generates molecular structures:
- Common catalysis molecules
- Adsorbates for surface science
- Organic building blocks
"""

from typing import Dict, Any, List, Optional
import numpy as np

try:
    from ase.build import molecule as ase_molecule
    HAS_ASE = True
except ImportError:
    HAS_ASE = False


COMMON_MOLECULES = {
    "H2": {"atoms": [("H", [0, 0, 0]), ("H", [0.74, 0, 0])]},
    "O2": {"atoms": [("O", [0, 0, 0]), ("O", [1.21, 0, 0])]},
    "N2": {"atoms": [("N", [0, 0, 0]), ("N", [1.10, 0, 0])]},
    "CO": {"atoms": [("C", [0, 0, 0]), ("O", [1.13, 0, 0])]},
    "CO2": {"atoms": [("O", [-1.16, 0, 0]), ("C", [0, 0, 0]), ("O", [1.16, 0, 0])]},
    "H2O": {"atoms": [("O", [0, 0, 0]), ("H", [0.96, 0, 0]), ("H", [-0.24, 0.93, 0])]},
    "NH3": {"atoms": [("N", [0, 0, 0]), ("H", [0, 0.94, 0.38]), ("H", [0.81, -0.47, 0.38]), ("H", [-0.81, -0.47, 0.38])]},
    "CH4": {"atoms": [("C", [0, 0, 0]), ("H", [0.63, 0.63, 0.63]), ("H", [-0.63, -0.63, 0.63]), ("H", [-0.63, 0.63, -0.63]), ("H", [0.63, -0.63, -0.63])]},
    "C2H4": {"atoms": [("C", [-0.67, 0, 0]), ("C", [0.67, 0, 0]), ("H", [-1.24, 0.92, 0]), ("H", [-1.24, -0.92, 0]), ("H", [1.24, 0.92, 0]), ("H", [1.24, -0.92, 0])]},
    "C2H6": {"atoms": [("C", [-0.77, 0, 0]), ("C", [0.77, 0, 0]), ("H", [-1.16, 0.5, 0.87]), ("H", [-1.16, 0.5, -0.87]), ("H", [-1.16, -1.0, 0]), ("H", [1.16, -0.5, 0.87]), ("H", [1.16, -0.5, -0.87]), ("H", [1.16, 1.0, 0])]},
    "CH3OH": {"atoms": [("C", [0, 0, 0]), ("O", [1.43, 0, 0]), ("H", [1.83, 0.87, 0]), ("H", [-0.4, 0.5, 0.87]), ("H", [-0.4, 0.5, -0.87]), ("H", [-0.4, -1.0, 0])]},
    "HCOOH": {"atoms": [("C", [0, 0, 0]), ("O", [1.2, 0, 0]), ("O", [-0.6, 1.1, 0]), ("H", [-1.1, 0.4, 0]), ("H", [-0.1, 1.4, 0])]},
}


def generate_molecule(
    name: str,
    optimize: bool = False
) -> Dict[str, Any]:
    """
    Generate molecular structure.
    
    Args:
        name: Molecule name or formula
        optimize: Whether to optimize geometry
    
    Returns:
        Molecule structure
    """
    name_upper = name.upper()
    
    # Try ASE first
    if HAS_ASE:
        try:
            mol = ase_molecule(name)
            positions = mol.get_positions()
            symbols = mol.get_chemical_symbols()
            
            return {
                "success": True,
                "name": name,
                "formula": mol.get_chemical_formula(),
                "n_atoms": len(mol),
                "structure": {
                    "atoms": [{"element": s, "cartesian": list(p)} for s, p in zip(symbols, positions)]
                }
            }
        except:
            pass
    
    # Fall back to database
    if name_upper in COMMON_MOLECULES:
        mol_data = COMMON_MOLECULES[name_upper]
        atoms = mol_data["atoms"]
        
        return {
            "success": True,
            "name": name,
            "formula": name_upper,
            "n_atoms": len(atoms),
            "structure": {
                "atoms": [{"element": e, "cartesian": list(p)} for e, p in atoms]
            }
        }
    
    return {
        "success": False,
        "error": {"code": "UNKNOWN", "message": f"Unknown molecule: {name}", "available": list(COMMON_MOLECULES.keys())}
    }
