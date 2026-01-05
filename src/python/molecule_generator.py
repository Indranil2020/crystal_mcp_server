#!/usr/bin/env python3
"""
Molecule Generator

Primary backend for molecule generation via MCP tools.
Supports ANY molecule through:
  - Local database (small molecules with accurate geometries)
  - Local aliases (PTCDA, aspirin, caffeine, etc.)
  - RDKit (SMILES → 3D with force field optimization)
  - PubChem API (~130M+ molecules by name/CID)
  - OPSIN (IUPAC names → SMILES)

Related modular generators:
  - generators.molecule.universal_molecule.generate_molecule_universal
  - generators.molecule.small_molecules.generate_molecule
  - generators.molecule.conformers.generate_conformers
"""

import sys
import json
import numpy as np
import importlib.util
import os
from typing import Dict, Any, Optional

# Robustness: Ensure src/python is in path
# This must happen BEFORE any imports from 'generators'
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Check availability using importlib.util.find_spec (no try/except)
UNIVERSAL_AVAILABLE = importlib.util.find_spec("generators.molecule.universal_molecule") is not None
ASE_AVAILABLE = importlib.util.find_spec("ase") is not None

# Import modules only if available
generate_molecule_universal = None
get_supported_molecules = None
ase_molecule = None
g2 = None
extra = {}

if UNIVERSAL_AVAILABLE:
    from generators.molecule.universal_molecule import generate_molecule_universal, get_supported_molecules

if ASE_AVAILABLE:
    from ase.build import molecule as ase_molecule
    from ase.collections import g2
    from ase.build.molecule import extra



def generate_molecule(
    name: str,
    input_type: str = "auto",
    optimize: bool = True,
    vacuum: float = 10.0,
    formula: Optional[str] = None,  # Legacy parameter
) -> Dict[str, Any]:
    """
    Generate an isolated molecule from any identifier.
    
    Args:
        name: Molecule identifier (name, SMILES, IUPAC, CID, InChI)
        input_type: Type hint ("auto", "name", "smiles", "iupac", "cid")
        optimize: Whether to optimize 3D geometry
        vacuum: Vacuum padding around the molecule
        formula: Legacy parameter (ignored, use name instead)
        
    Returns:
        Structure dictionary with molecule data
    """
    # Try universal molecule generation first (supports ~130M+ molecules)
    if UNIVERSAL_AVAILABLE:
        result = generate_molecule_universal(
            identifier=name,
            input_type=input_type,
            optimize=optimize,
            allow_external=True
        )
        
        if result.get("success"):
            # Convert to MCP structure format
            return {
                "success": True,
                "structure": molecule_to_structure_dict(result, vacuum),
                "source": result.get("source", "universal"),
                "metadata": {
                    "formula": result.get("formula", name),
                    "smiles": result.get("smiles"),
                    "canonical_smiles": result.get("canonical_smiles"),
                    "iupac_name": result.get("iupac_name"),
                    "pubchem_cid": result.get("pubchem_cid"),
                    "optimized": result.get("optimized", optimize),
                }
            }
    
    # Fallback to ASE for molecules in g2 database
    if ASE_AVAILABLE:
        if g2 and name in g2.names:
            try:
                atoms = ase_molecule(name)
                atoms.center(vacuum=vacuum)
                return {
                    "success": True,
                    "structure": atoms_to_dict(atoms),
                    "source": "ase_g2"
                }
            except Exception as e:
                pass
        
        if name in extra:
            try:
                atoms = ase_molecule(name)
                atoms.center(vacuum=vacuum)
                return {
                    "success": True,
                    "structure": atoms_to_dict(atoms),
                    "source": "ase_extra"
                }
            except Exception:
                pass
    
    # All methods failed
    error_msg = f"Unknown molecule '{name}'"
    suggestions = [
        "Check spelling of the molecule name",
        "Try providing a SMILES string (e.g., 'c1ccccc1' for benzene)",
        "Use IUPAC systematic name",
        f"Search PubChem for '{name}' to find the correct identifier",
    ]
    
    available = []
    if ASE_AVAILABLE and g2:
        available = sorted(list(g2.names))[:20]
    
    if UNIVERSAL_AVAILABLE:
        info = get_supported_molecules()
        available.extend(info.get("local_database", []))
        available.extend(info.get("local_aliases", []))
        available = sorted(set(available))[:30]
        suggestions.insert(0, "Universal molecule generation available - try any common molecule name")
    
    return {
        "success": False,
        "error": {
            "code": "MOLECULE_NOT_FOUND",
            "message": error_msg,
            "suggestions": suggestions,
            "available_examples": available,
            "capabilities": {
                "universal_available": UNIVERSAL_AVAILABLE,
                "ase_available": ASE_AVAILABLE,
            }
        }
    }


def molecule_to_structure_dict(mol_result: Dict, vacuum: float = 10.0) -> Dict[str, Any]:
    """Convert universal molecule result to MCP structure format."""
    atoms = mol_result.get("atoms", [])
    coords = mol_result.get("coords", [])
    n_atoms = len(atoms)
    
    if n_atoms == 0:
        return {}
    
    # Calculate bounding box and add vacuum
    coords_array = np.array(coords)
    min_coords = coords_array.min(axis=0)
    max_coords = coords_array.max(axis=0)
    box_size = max_coords - min_coords + 2 * vacuum
    
    # Center molecule in box
    center_offset = box_size / 2 - (max_coords + min_coords) / 2
    centered_coords = coords_array + center_offset
    
    # Create cell matrix
    cell = np.diag(box_size)
    
    # Create sites
    sites = []
    for i, (atom, coord) in enumerate(zip(atoms, centered_coords)):
        frac_coord = coord / box_size
        sites.append({
            "element": atom,
            "coords": frac_coord.tolist(),
            "cartesian": coord.tolist(),
            "species": [{"element": atom, "occupation": 1.0}]
        })
    
    return {
        "lattice": {
            "matrix": cell.tolist(),
            "volume": float(np.prod(box_size)),
            "a": float(box_size[0]),
            "b": float(box_size[1]),
            "c": float(box_size[2]),
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0
        },
        "atoms": sites,
        "sites": sites,
        "space_group": {
            "number": 1,
            "symbol": "P1",
            "crystal_system": "triclinic"
        },
        "metadata": {
            "formula": mol_result.get("formula", ""),
            "natoms": n_atoms,
            "pbc": [False, False, False],
            "molecular_weight": mol_result.get("molecular_weight"),
            "source": mol_result.get("source"),
        }
    }


def atoms_to_dict(atoms) -> Dict[str, Any]:
    """Convert ASE atoms to dictionary (legacy fallback)."""
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
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0
        },
        "atoms": sites,
        "sites": sites,
        "space_group": {
            "number": 1,
            "symbol": "P1",
            "crystal_system": "triclinic"
        },
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
        
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
            print(json.dumps({"success": False, "error": f"Input file not found: {input_file}"}))
            sys.exit(1)

    with open(input_file, 'r') as f:
        params = json.load(f)
        
    result = generate_molecule(
        name=params.get("name"),
        input_type=params.get("input_type", "auto"),
        optimize=params.get("optimize", True),
        vacuum=float(params.get("vacuum", 10.0)),
        formula=params.get("formula"),  # Legacy
    )
    print(json.dumps(result))


if __name__ == "__main__":
    main()


