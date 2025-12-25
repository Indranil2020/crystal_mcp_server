#!/usr/bin/env python3
"""
Molecule Generator

Primary backend for molecule generation via MCP tools.
Generates isolated molecular structures using ASE.

Related modular generators (for advanced use):
  - generators.molecule.small_molecules.generate_molecule
  - generators.molecule.conformers.generate_conformers
"""

import sys
import json
import numpy as np
from ase.build import molecule
from ase.collections import g2
from ase.build.molecule import extra
from typing import Dict, Any, List, Optional

def generate_molecule(
    name: str,
    formula: Optional[str] = None,
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Generate an isolated molecule.
    
    Args:
        name: Name of the molecule (e.g., "H2O", "C60", "Benzene")
        formula: Optional chemical formula if name is generic (not always supported by ASE molecule)
        vacuum: Vacuum padding around the molecule
        
    Returns:
        Structure dictionary
    """
    atoms = None
    
    # Check if name is in known molecule databases
    if name not in g2.names and name not in extra:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MOLECULE",
                "message": f"Unknown molecule '{name}'",
                "details": {"available": sorted(list(g2.names))[:20]}
            }
        }

    atoms = molecule(name)
         
    if atoms is None:
         return {
            "success": False,
            "error": {
                "code": "GENERATION_FAILED",
                "message": f"Failed to generate molecule '{name}'"
            }
        }
        
    # Center and add vacuum
    atoms.center(vacuum=vacuum)
    
    # Convert to dictionary
    return {
        "success": True,
        "structure": atoms_to_dict(atoms)
    }

def atoms_to_dict(atoms) -> Dict[str, Any]:
    """Convert ASE atoms to dictionary."""
    cell = atoms.get_cell()
    sites = [
        {
            "element": atom.symbol,
            "coords": atoms.get_scaled_positions()[i].tolist(),
            "cartesian": atom.position.tolist(),
            "species": [{"element": atom.symbol, "occupation": 1.0}]
        }
        for i, atom in enumerate(atoms)
    ]

    return {
        "lattice": {
            "matrix": cell.tolist(),
            "volume": float(atoms.get_volume()),
            "a": float(np.linalg.norm(cell[0])),
            "b": float(np.linalg.norm(cell[1])),
            "c": float(np.linalg.norm(cell[2])),
            "alpha": 90.0,  # Approximate for box
            "beta": 90.0,
            "gamma": 90.0
        },
        "atoms": sites,
        "sites": sites,
        "metadata": {
            "formula": atoms.get_chemical_formula(),
            "natoms": len(atoms),
            "pbc": atoms.pbc.tolist()
        }
    }

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "Usage: python molecule_generator.py <input.json>"}))
        sys.exit(1)
        
    with open(sys.argv[1], 'r') as f:
        params = json.load(f)
        
    result = generate_molecule(
        name=params.get("name"),
        formula=params.get("formula"),
        vacuum=float(params.get("vacuum", 10.0))
    )
    print(json.dumps(result))

if __name__ == "__main__":
    main()
