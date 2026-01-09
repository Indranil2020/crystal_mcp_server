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

PERFORMANCE: Uses LRU caching for repeated lookups.
  - Cache size: 256 entries
  - Cache hit/miss statistics tracked
  - Debug logging for all cache operations

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
import time
from functools import lru_cache
from typing import Dict, Any, Optional, Tuple

# Debug logging for cache operations
def _debug_log(action: str, details: str = "") -> None:
    """Print debug message with timestamp for cache operations."""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[MOLECULE_CACHE] {timestamp} {action}: {details}", file=sys.stderr)

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
    _debug_log("INIT", "Universal molecule generator available")

if ASE_AVAILABLE:
    from ase.build import molecule as ase_molecule
    from ase.collections import g2
    from ase.build.molecule import extra
    _debug_log("INIT", f"ASE available with {len(g2.names) if g2 else 0} G2 molecules")

# Cache statistics tracking
_cache_stats = {"hits": 0, "misses": 0}


@lru_cache(maxsize=256)
def _cached_molecule_lookup(name: str, input_type: str, optimize: bool) -> Tuple[bool, str]:
    """
    Cached molecule lookup - avoids repeated RDKit/PubChem calls.
    
    Returns tuple of (success, json_result) since lru_cache needs hashable return values.
    The JSON is parsed by the caller.
    """
    _debug_log("CACHE_MISS", f"name='{name}' type='{input_type}' optimize={optimize}")
    _cache_stats["misses"] += 1
    
    result = None
    source = None
    
    # Try universal molecule generation first (supports ~130M+ molecules)
    if UNIVERSAL_AVAILABLE:
        result = generate_molecule_universal(
            identifier=name,
            input_type=input_type,
            optimize=optimize,
            allow_external=True
        )
        
        if result.get("success"):
            source = result.get("source", "universal")
            return (True, json.dumps({"result": result, "source": source}))
    
    # Fallback to ASE for molecules in g2 database
    if ASE_AVAILABLE:
        if g2 and name in g2.names:
            atoms = ase_molecule(name)
            atoms.center(vacuum=10.0)
            result = {
                "atoms": [a.symbol for a in atoms],
                "coords": atoms.get_positions().tolist(),
                "formula": atoms.get_chemical_formula(),
            }
            return (True, json.dumps({"result": result, "source": "ase_g2", "ase_atoms": True}))
        
        if name in extra:
            atoms = ase_molecule(name)
            atoms.center(vacuum=10.0)
            result = {
                "atoms": [a.symbol for a in atoms],
                "coords": atoms.get_positions().tolist(),
                "formula": atoms.get_chemical_formula(),
            }
            return (True, json.dumps({"result": result, "source": "ase_extra", "ase_atoms": True}))
    
    # All methods failed
    return (False, json.dumps({"error": f"Unknown molecule '{name}'"}))


def get_cache_stats() -> Dict[str, Any]:
    """Get current cache statistics."""
    cache_info = _cached_molecule_lookup.cache_info()
    return {
        "hits": _cache_stats["hits"],
        "misses": _cache_stats["misses"],
        "cache_size": cache_info.currsize,
        "max_size": cache_info.maxsize,
        "hit_rate": _cache_stats["hits"] / max(1, _cache_stats["hits"] + _cache_stats["misses"]),
    }


def clear_cache() -> Dict[str, Any]:
    """Clear the molecule lookup cache."""
    _cached_molecule_lookup.cache_clear()
    _cache_stats["hits"] = 0
    _cache_stats["misses"] = 0
    _debug_log("CACHE_CLEAR", "Cache cleared")
    return {"success": True, "message": "Cache cleared"}


def generate_molecule(
    name: str,
    input_type: str = "auto",
    optimize: bool = True,
    vacuum: float = 10.0,
    formula: Optional[str] = None,  # Legacy parameter
) -> Dict[str, Any]:
    """
    Generate an isolated molecule from any identifier.
    
    Uses LRU caching for repeated lookups - cache stats available via get_cache_stats().
    
    Args:
        name: Molecule identifier (name, SMILES, IUPAC, CID, InChI)
        input_type: Type hint ("auto", "name", "smiles", "iupac", "cid")
        optimize: Whether to optimize 3D geometry
        vacuum: Vacuum padding around the molecule
        formula: Legacy parameter (ignored, use name instead)
        
    Returns:
        Structure dictionary with molecule data
    """
    _debug_log("REQUEST", f"name='{name}' type='{input_type}' optimize={optimize} vacuum={vacuum}")
    
    # Check cache first (for cache hit tracking)
    cache_info = _cached_molecule_lookup.cache_info()
    cache_hits_before = cache_info.hits
    
    # Use cached lookup (vacuum is not part of cache key since it's post-processing)
    success, result_json = _cached_molecule_lookup(name, input_type, optimize)
    
    # Track cache hits
    cache_info_after = _cached_molecule_lookup.cache_info()
    if cache_info_after.hits > cache_hits_before:
        _cache_stats["hits"] += 1
        _debug_log("CACHE_HIT", f"name='{name}' (total hits: {_cache_stats['hits']})")
    
    if success:
        cached_data = json.loads(result_json)
        result = cached_data["result"]
        source = cached_data["source"]
        
        # Handle ASE-returned results
        if cached_data.get("ase_atoms"):
            # Need to reconstruct structure dict from minimal cached data
            structure = _ase_result_to_structure(result, vacuum)
            return {
                "success": True,
                "structure": structure,
                "source": source,
                "cache_stats": get_cache_stats(),
            }
        
        # Standard universal molecule result
        return {
            "success": True,
            "structure": molecule_to_structure_dict(result, vacuum),
            "source": source,
            "metadata": {
                "formula": result.get("formula", name),
                "smiles": result.get("smiles"),
                "canonical_smiles": result.get("canonical_smiles"),
                "iupac_name": result.get("iupac_name"),
                "pubchem_cid": result.get("pubchem_cid"),
                "optimized": result.get("optimized", optimize),
            },
            "cache_stats": get_cache_stats(),
        }
    
    # All methods failed
    error_msg = f"Unknown molecule '{name}'"
    _debug_log("ERROR", error_msg)
    
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
        },
        "cache_stats": get_cache_stats(),
    }


def _ase_result_to_structure(result: Dict, vacuum: float) -> Dict[str, Any]:
    """Convert cached ASE molecule result to structure dict."""
    atoms = result.get("atoms", [])
    coords = result.get("coords", [])
    n_atoms = len(atoms)
    
    if n_atoms == 0:
        return {}
    
    coords_array = np.array(coords)
    min_coords = coords_array.min(axis=0)
    max_coords = coords_array.max(axis=0)
    box_size = max_coords - min_coords + 2 * vacuum
    
    center_offset = box_size / 2 - (max_coords + min_coords) / 2
    centered_coords = coords_array + center_offset
    
    cell = np.diag(box_size)
    
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
            "formula": result.get("formula", ""),
            "natoms": n_atoms,
            "pbc": [False, False, False],
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


