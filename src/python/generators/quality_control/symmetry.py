"""
quality_control/symmetry.py - Symmetry Analysis

Analyzes crystal symmetry:
- Space group detection (spglib)
- Point group analysis
- Tolerance sweeps
"""

from typing import Dict, Any, List, Optional
import numpy as np


def analyze_symmetry(
    structure: Dict[str, Any],
    tolerance: float = 0.01
) -> Dict[str, Any]:
    """
    Analyze symmetry of structure.
    
    Args:
        structure: Structure dictionary
        tolerance: Symmetry tolerance in Angstrom
    
    Returns:
        Symmetry analysis results
    """
    try:
        from pymatgen.core import Structure, Lattice
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        
        # Reconstruct structure
        lattice_info = structure.get("structure", structure).get("lattice", {})
        atoms = structure.get("structure", structure).get("atoms", [])
        
        if not atoms:
            return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}
        
        # Build lattice
        if "matrix" in lattice_info:
            lattice = Lattice(lattice_info["matrix"])
        else:
            a = lattice_info.get("a", 5.0)
            b = lattice_info.get("b", a)
            c = lattice_info.get("c", a)
            lattice = Lattice.orthorhombic(a, b, c)
        
        species = [a["element"] for a in atoms]
        coords = [a["coords"] for a in atoms]
        
        pmg_structure = Structure(lattice, species, coords)
        
        # Analyze symmetry
        analyzer = SpacegroupAnalyzer(pmg_structure, symprec=tolerance)
        
        return {
            "success": True,
            "space_group_symbol": analyzer.get_space_group_symbol(),
            "space_group_number": analyzer.get_space_group_number(),
            "point_group": analyzer.get_point_group_symbol(),
            "crystal_system": analyzer.get_crystal_system(),
            "hall_symbol": analyzer.get_hall(),
            "tolerance_used": tolerance,
            "n_symmetry_operations": len(analyzer.get_symmetry_operations()),
            "is_laue": analyzer.is_laue(),
        }
    
    except ImportError:
        return {"success": False, "error": {"code": "IMPORT_ERROR", "message": "spglib not available"}}
    except Exception as e:
        return {"success": False, "error": {"code": "ANALYSIS_ERROR", "message": str(e)}}


def find_space_group(
    structure: Dict[str, Any],
    method: str = "spglib"
) -> Dict[str, Any]:
    """
    Find space group of structure.
    
    Args:
        structure: Structure dictionary
        method: 'spglib' or 'findsym'
    
    Returns:
        Space group information
    """
    return analyze_symmetry(structure, tolerance=0.01)


def check_tolerance(
    structure: Dict[str, Any],
    tolerances: List[float] = [0.001, 0.01, 0.1, 0.5]
) -> Dict[str, Any]:
    """
    Sweep symmetry tolerance to find optimal value.
    
    Args:
        structure: Structure dictionary
        tolerances: List of tolerances to try
    
    Returns:
        Results for each tolerance
    """
    results = []
    
    for tol in tolerances:
        result = analyze_symmetry(structure, tolerance=tol)
        if result["success"]:
            results.append({
                "tolerance": tol,
                "space_group_symbol": result["space_group_symbol"],
                "space_group_number": result["space_group_number"],
                "n_ops": result["n_symmetry_operations"]
            })
    
    return {
        "success": True,
        "tolerance_sweep": results,
        "recommended_tolerance": tolerances[1] if len(tolerances) > 1 else tolerances[0]
    }
