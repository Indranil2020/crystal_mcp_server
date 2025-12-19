"""
molecule/small_molecules.py - Small Molecule Generation

Comprehensive small molecule generation for adsorption and reaction studies.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Molecule


# Small molecule database with accurate geometries
MOLECULE_DATABASE = {
    # Diatomics
    "H2": {"atoms": ["H", "H"], "coords": [[0, 0, 0], [0.74, 0, 0]], "bond_order": 1, "spin": 0},
    "O2": {"atoms": ["O", "O"], "coords": [[0, 0, 0], [1.21, 0, 0]], "bond_order": 2, "spin": 2},
    "N2": {"atoms": ["N", "N"], "coords": [[0, 0, 0], [1.10, 0, 0]], "bond_order": 3, "spin": 0},
    "CO": {"atoms": ["C", "O"], "coords": [[0, 0, 0], [1.13, 0, 0]], "bond_order": 3, "spin": 0},
    "NO": {"atoms": ["N", "O"], "coords": [[0, 0, 0], [1.15, 0, 0]], "bond_order": 2.5, "spin": 1},
    "HF": {"atoms": ["H", "F"], "coords": [[0, 0, 0], [0.92, 0, 0]], "bond_order": 1, "spin": 0},
    "HCl": {"atoms": ["H", "Cl"], "coords": [[0, 0, 0], [1.27, 0, 0]], "bond_order": 1, "spin": 0},
    "Cl2": {"atoms": ["Cl", "Cl"], "coords": [[0, 0, 0], [1.99, 0, 0]], "bond_order": 1, "spin": 0},
    "F2": {"atoms": ["F", "F"], "coords": [[0, 0, 0], [1.42, 0, 0]], "bond_order": 1, "spin": 0},
    
    # Triatomics - bent
    "H2O": {"atoms": ["O", "H", "H"], 
            "coords": [[0, 0, 0], [0.96*0.866, 0.96*0.5, 0], [-0.96*0.866, 0.96*0.5, 0]],
            "angle_deg": 104.5, "spin": 0},
    "H2S": {"atoms": ["S", "H", "H"], 
            "coords": [[0, 0, 0], [1.34*0.819, 1.34*0.574, 0], [-1.34*0.819, 1.34*0.574, 0]],
            "angle_deg": 92.1, "spin": 0},
    "SO2": {"atoms": ["S", "O", "O"], 
            "coords": [[0, 0, 0], [1.43, 0.52, 0], [-1.43, 0.52, 0]],
            "angle_deg": 119.5, "spin": 0},
    "NO2": {"atoms": ["N", "O", "O"], 
            "coords": [[0, 0, 0], [1.19, 0.46, 0], [-1.19, 0.46, 0]],
            "angle_deg": 134.3, "spin": 1},
    "O3": {"atoms": ["O", "O", "O"], 
           "coords": [[0, 0, 0], [1.28, 0, 0], [-0.64, 1.11, 0]],
           "angle_deg": 116.8, "spin": 0},
    
    # Triatomics - linear
    "CO2": {"atoms": ["C", "O", "O"], 
            "coords": [[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
            "linear": True, "spin": 0},
    "N2O": {"atoms": ["N", "N", "O"], 
            "coords": [[0, 0, 0], [1.13, 0, 0], [2.31, 0, 0]],
            "linear": True, "spin": 0},
    "HCN": {"atoms": ["H", "C", "N"], 
            "coords": [[0, 0, 0], [1.06, 0, 0], [2.22, 0, 0]],
            "linear": True, "spin": 0},
    "C2H2": {"atoms": ["C", "C", "H", "H"], 
             "coords": [[0, 0, 0], [1.20, 0, 0], [-1.06, 0, 0], [2.26, 0, 0]],
             "linear": True, "spin": 0},
    
    # Tetrahedral
    "CH4": {"atoms": ["C", "H", "H", "H", "H"], 
            "coords": [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], 
                      [0.63, -0.63, -0.63], [-0.63, 0.63, -0.63]],
            "symmetry": "Td", "spin": 0},
    "NH3": {"atoms": ["N", "H", "H", "H"], 
            "coords": [[0, 0, 0], [0.94, 0.38, 0], [-0.47, 0.38, 0.81], [-0.47, 0.38, -0.81]],
            "symmetry": "C3v", "spin": 0},
    "PH3": {"atoms": ["P", "H", "H", "H"], 
            "coords": [[0, 0, 0], [1.42, 0.30, 0], [-0.71, 0.30, 1.23], [-0.71, 0.30, -1.23]],
            "symmetry": "C3v", "spin": 0},
    "SiH4": {"atoms": ["Si", "H", "H", "H", "H"], 
             "coords": [[0, 0, 0], [0.85, 0.85, 0.85], [-0.85, -0.85, 0.85], 
                       [0.85, -0.85, -0.85], [-0.85, 0.85, -0.85]],
             "symmetry": "Td", "spin": 0},
    
    # Planar
    "C2H4": {"atoms": ["C", "C", "H", "H", "H", "H"], 
             "coords": [[0, 0, 0], [1.34, 0, 0], 
                       [-0.55, 0.93, 0], [-0.55, -0.93, 0],
                       [1.89, 0.93, 0], [1.89, -0.93, 0]],
             "symmetry": "D2h", "spin": 0},
    "HCHO": {"atoms": ["C", "O", "H", "H"], 
             "coords": [[0, 0, 0], [1.21, 0, 0], [-0.54, 0.93, 0], [-0.54, -0.93, 0]],
             "symmetry": "C2v", "spin": 0},
    "HCOOH": {"atoms": ["C", "O", "O", "H", "H"], 
              "coords": [[0, 0, 0], [1.20, 0.15, 0], [-0.35, 1.20, 0], 
                        [-0.65, -0.85, 0], [1.30, 1.05, 0]],
              "spin": 0},
    
    # Alcohols and ethers
    "CH3OH": {"atoms": ["C", "O", "H", "H", "H", "H"], 
              "coords": [[0, 0, 0], [1.43, 0, 0], 
                        [-0.38, 1.03, 0], [-0.38, -0.51, 0.89], [-0.38, -0.51, -0.89],
                        [1.79, 0.87, 0]],
              "spin": 0},
    "C2H5OH": {"atoms": ["C", "C", "O", "H", "H", "H", "H", "H", "H"], 
               "coords": [[0, 0, 0], [1.52, 0, 0], [2.20, 1.30, 0],
                         [-0.38, 1.03, 0], [-0.38, -0.51, 0.89], [-0.38, -0.51, -0.89],
                         [1.90, -0.51, 0.89], [1.90, -0.51, -0.89], [3.10, 1.20, 0]],
               "spin": 0},
    
    # Radicals and intermediates
    "OH": {"atoms": ["O", "H"], "coords": [[0, 0, 0], [0.97, 0, 0]], "spin": 1, "radical": True},
    "CH3": {"atoms": ["C", "H", "H", "H"], 
            "coords": [[0, 0, 0], [1.08, 0, 0], [-0.54, 0.94, 0], [-0.54, -0.94, 0]],
            "spin": 1, "radical": True},
    "CHO": {"atoms": ["C", "H", "O"], 
            "coords": [[0, 0, 0], [-0.60, 0.90, 0], [1.20, 0, 0]],
            "spin": 1, "radical": True},
    "OOH": {"atoms": ["O", "O", "H"], 
            "coords": [[0, 0, 0], [1.33, 0, 0], [1.68, 0.93, 0]],
            "spin": 1, "radical": True},
    "COOH": {"atoms": ["C", "O", "O", "H"], 
             "coords": [[0, 0, 0], [1.20, 0.15, 0], [-0.35, 1.20, 0], [1.30, 1.05, 0]],
             "spin": 1, "radical": True},
}


def generate_molecule(
    formula: str,
    optimize: bool = False
) -> Dict[str, Any]:
    """
    Generate small molecule structure.
    
    Args:
        formula: Molecule formula from database
        optimize: Whether to optimize geometry (placeholder)
    
    Returns:
        Molecule structure
    """
    if formula not in MOLECULE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MOLECULE", "message": f"Unknown molecule '{formula}'",
                      "available": list(MOLECULE_DATABASE.keys())}
        }
    
    info = MOLECULE_DATABASE[formula]
    
    mol = Molecule(info["atoms"], info["coords"])
    
    return {
        "success": True,
        "formula": formula,
        "n_atoms": len(info["atoms"]),
        "atoms": info["atoms"],
        "coords": info["coords"],
        "spin_multiplicity": 2 * info.get("spin", 0) + 1,
        "is_radical": info.get("radical", False),
        "symmetry": info.get("symmetry", ""),
        "is_linear": info.get("linear", False),
        "bond_angle_deg": info.get("angle_deg", 0),
        "molecule": {
            "species": [str(s) for s in mol.species],
            "coords": [list(c) for c in mol.cart_coords]
        }
    }


def get_molecule_database() -> Dict[str, Any]:
    """Get all available molecules."""
    return {
        "success": True,
        "n_molecules": len(MOLECULE_DATABASE),
        "molecules": list(MOLECULE_DATABASE.keys()),
        "by_type": {
            "diatomic": [k for k in MOLECULE_DATABASE if len(MOLECULE_DATABASE[k]["atoms"]) == 2],
            "triatomic": [k for k in MOLECULE_DATABASE if len(MOLECULE_DATABASE[k]["atoms"]) == 3],
            "tetrahedral": [k for k in MOLECULE_DATABASE if MOLECULE_DATABASE[k].get("symmetry") in ["Td", "C3v"]],
            "radicals": [k for k in MOLECULE_DATABASE if MOLECULE_DATABASE[k].get("radical", False)],
        }
    }
