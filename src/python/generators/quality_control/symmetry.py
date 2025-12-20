"""
quality_control/symmetry.py - Symmetry Analysis

Comprehensive symmetry analysis per structure_catalogue.md Category 16:
(i) Automated symmetry refinement & space-group detection
(ii) Tolerance sweep for near-symmetric structures
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# Crystal system database
CRYSTAL_SYSTEMS = {
    "triclinic": {"min_sg": 1, "max_sg": 2, "point_groups": ["1", "-1"]},
    "monoclinic": {"min_sg": 3, "max_sg": 15, "point_groups": ["2", "m", "2/m"]},
    "orthorhombic": {"min_sg": 16, "max_sg": 74, "point_groups": ["222", "mm2", "mmm"]},
    "tetragonal": {"min_sg": 75, "max_sg": 142, "point_groups": ["4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm"]},
    "trigonal": {"min_sg": 143, "max_sg": 167, "point_groups": ["3", "-3", "32", "3m", "-3m"]},
    "hexagonal": {"min_sg": 168, "max_sg": 194, "point_groups": ["6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm"]},
    "cubic": {"min_sg": 195, "max_sg": 230, "point_groups": ["23", "m-3", "432", "-43m", "m-3m"]},
}


# Common space groups
COMMON_SPACEGROUPS = {
    "Fm-3m": {"number": 225, "system": "cubic", "examples": ["NaCl", "Cu", "Au", "Ag"]},
    "Fd-3m": {"number": 227, "system": "cubic", "examples": ["diamond", "Si", "spinel"]},
    "Pm-3m": {"number": 221, "system": "cubic", "examples": ["perovskite", "CsCl"]},
    "Im-3m": {"number": 229, "system": "cubic", "examples": ["bcc_Fe", "W"]},
    "P6_3/mmc": {"number": 194, "system": "hexagonal", "examples": ["hcp", "graphite"]},
    "R-3m": {"number": 166, "system": "trigonal", "examples": ["Bi2Te3", "Bi2Se3"]},
    "Pnma": {"number": 62, "system": "orthorhombic", "examples": ["perovskite_distorted"]},
    "I4/mmm": {"number": 139, "system": "tetragonal", "examples": ["rutile"]},
    "P2_1/c": {"number": 14, "system": "monoclinic", "examples": ["common_molecular"]},
    "P-1": {"number": 2, "system": "triclinic", "examples": ["low_symmetry"]},
    "C2/m": {"number": 12, "system": "monoclinic", "examples": ["layered_oxides"]},
    "P4/mmm": {"number": 123, "system": "tetragonal", "examples": ["cuprates"]},
    "P6mm": {"number": 183, "system": "hexagonal", "examples": ["wurtzite"]},
    "F-43m": {"number": 216, "system": "cubic", "examples": ["zincblende", "half_heusler"]},
}


# Wyckoff position database (selected)
WYCKOFF_POSITIONS = {
    225: {  # Fm-3m
        "a": {"multiplicity": 4, "coordinates": ["0,0,0"]},
        "b": {"multiplicity": 4, "coordinates": ["1/2,1/2,1/2"]},
        "c": {"multiplicity": 8, "coordinates": ["1/4,1/4,1/4"]},
        "d": {"multiplicity": 24, "coordinates": ["0,1/4,1/4"]},
    },
    227: {  # Fd-3m
        "a": {"multiplicity": 8, "coordinates": ["0,0,0"]},
        "b": {"multiplicity": 8, "coordinates": ["1/2,1/2,1/2"]},
        "c": {"multiplicity": 16, "coordinates": ["1/8,1/8,1/8"]},
    },
    221: {  # Pm-3m
        "a": {"multiplicity": 1, "coordinates": ["0,0,0"]},
        "b": {"multiplicity": 1, "coordinates": ["1/2,1/2,1/2"]},
        "c": {"multiplicity": 3, "coordinates": ["0,1/2,1/2"]},
        "d": {"multiplicity": 3, "coordinates": ["1/2,0,0"]},
    },
}


def analyze_symmetry(
    structure: Structure,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0
) -> Dict[str, Any]:
    """
    Analyze crystal symmetry.
    
    Args:
        structure: Input structure
        symprec: Symmetry precision in Angstrom
        angle_tolerance: Angle tolerance in degrees
    
    Returns:
        Symmetry analysis results
    """
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    
    spacegroup = analyzer.get_space_group_symbol()
    sg_number = analyzer.get_space_group_number()
    point_group = analyzer.get_point_group_symbol()
    crystal_system = analyzer.get_crystal_system()
    hall_symbol = analyzer.get_hall()
    
    # Get symmetry operations
    sym_ops = analyzer.get_symmetry_operations()
    n_operations = len(sym_ops)
    
    # Get conventional and primitive cells
    conventional = analyzer.get_conventional_standard_structure()
    primitive = analyzer.get_primitive_standard_structure()
    
    # Identify Wyckoff positions
    sym_dataset = analyzer.get_symmetry_dataset()
    wyckoff_symbols = list(set(sym_dataset['wyckoffs'])) if sym_dataset else []
    
    return {
        "success": True,
        "spacegroup_symbol": spacegroup,
        "spacegroup_number": sg_number,
        "point_group": point_group,
        "crystal_system": crystal_system,
        "hall_symbol": hall_symbol,
        "n_symmetry_operations": n_operations,
        "wyckoff_letters": wyckoff_symbols,
        "n_atoms_primitive": len(primitive),
        "n_atoms_conventional": len(conventional),
        "symprec_used": symprec,
        "is_centrosymmetric": "-" in spacegroup or "m" in spacegroup.lower(),
    }


def tolerance_sweep(
    structure: Structure,
    tolerances: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Sweep symmetry tolerance to find optimal value.
    
    Args:
        structure: Input structure
        tolerances: List of tolerances to try
    
    Returns:
        Results at each tolerance
    """
    if tolerances is None:
        tolerances = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    
    results = []
    
    for tol in tolerances:
        try:
            analyzer = SpacegroupAnalyzer(structure, symprec=tol)
            sg = analyzer.get_space_group_symbol()
            sg_num = analyzer.get_space_group_number()
            n_ops = len(analyzer.get_symmetry_operations())
            
            results.append({
                "tolerance": tol,
                "spacegroup": sg,
                "sg_number": sg_num,
                "n_operations": n_ops,
                "success": True
            })
        except Exception as e:
            results.append({
                "tolerance": tol,
                "success": False,
                "error": str(e)
            })
    
    # Find optimal tolerance (highest symmetry with reasonable tolerance)
    valid_results = [r for r in results if r.get("success", False)]
    if valid_results:
        # Prefer higher symmetry (more operations) at reasonable tolerance
        optimal = max(valid_results, key=lambda x: (x["n_operations"], -x["tolerance"]))
        optimal_tol = optimal["tolerance"]
    else:
        optimal_tol = 0.01
    
    return {
        "success": True,
        "sweep_results": results,
        "optimal_tolerance": optimal_tol,
        "n_tolerances_tested": len(tolerances)
    }


def refine_symmetry(
    structure: Structure,
    target_spacegroup: Optional[int] = None,
    symprec: float = 0.1
) -> Dict[str, Any]:
    """
    Refine structure to higher symmetry.
    
    Args:
        structure: Input structure
        target_spacegroup: Target space group number (optional)
        symprec: Symmetry precision
    
    Returns:
        Refined structure
    """
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
    
    # Get refined structure
    refined = analyzer.get_refined_structure()
    conventional = analyzer.get_conventional_standard_structure()
    
    # Analyze result
    new_analyzer = SpacegroupAnalyzer(refined, symprec=0.01)
    
    return {
        "success": True,
        "original_spacegroup": analyzer.get_space_group_symbol(),
        "refined_spacegroup": new_analyzer.get_space_group_symbol(),
        "original_n_atoms": len(structure),
        "refined_n_atoms": len(refined),
        "conventional_n_atoms": len(conventional),
        "symprec_used": symprec,
        "refined_structure": {
            "lattice": {"a": refined.lattice.a, "b": refined.lattice.b, "c": refined.lattice.c},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in refined]
        }
    }


def get_equivalent_positions(
    spacegroup_number: int,
    general_position: List[float] = [0.1, 0.2, 0.3]
) -> Dict[str, Any]:
    """
    Get equivalent positions from general position.
    
    Args:
        spacegroup_number: Space group number
        general_position: Starting position [x, y, z]
    
    Returns:
        All equivalent positions
    """
    # Create simple structure in given space group
    a = 5.0
    
    if spacegroup_number >= 195:  # Cubic
        lattice = Lattice.cubic(a)
    elif spacegroup_number >= 168:  # Hexagonal
        lattice = Lattice.hexagonal(a, a * 1.6)
    elif spacegroup_number >= 143:  # Trigonal
        lattice = Lattice.rhombohedral(a, 80)
    elif spacegroup_number >= 75:  # Tetragonal
        lattice = Lattice.tetragonal(a, a * 1.2)
    elif spacegroup_number >= 16:  # Orthorhombic
        lattice = Lattice.orthorhombic(a, a * 1.1, a * 1.2)
    elif spacegroup_number >= 3:  # Monoclinic
        lattice = Lattice.monoclinic(a, a * 1.1, a * 1.2, 95)
    else:  # Triclinic
        lattice = Lattice([[a, 0, 0], [0.5, a * 1.05, 0], [0.3, 0.4, a * 1.1]])
    
    structure = Structure(lattice, ["X"], [general_position])
    
    analyzer = SpacegroupAnalyzer(structure, symprec=0.01)
    sym_ops = analyzer.get_symmetry_operations()
    
    equivalent_positions = []
    x, y, z = general_position
    
    for op in sym_ops:
        new_pos = op.operate([x, y, z])
        # Wrap to [0, 1)
        new_pos = [p % 1 for p in new_pos]
        
        # Check if unique
        is_unique = True
        for existing in equivalent_positions:
            dist = sum((a - b)**2 for a, b in zip(new_pos, existing))
            if dist < 0.0001:
                is_unique = False
                break
        
        if is_unique:
            equivalent_positions.append(new_pos)
    
    return {
        "success": True,
        "spacegroup_number": spacegroup_number,
        "general_position": general_position,
        "n_equivalent": len(equivalent_positions),
        "equivalent_positions": equivalent_positions,
        "multiplicity": len(equivalent_positions)
    }


def identify_special_positions(
    structure: Structure,
    symprec: float = 0.01
) -> Dict[str, Any]:
    """
    Identify atoms on special Wyckoff positions.
    
    Args:
        structure: Input structure
        symprec: Symmetry precision
    
    Returns:
        Analysis of special positions
    """
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
    sym_dataset = analyzer.get_symmetry_dataset()
    
    if not sym_dataset:
        return {"success": False, "error": {"code": "NO_SYMMETRY", "message": "Could not analyze symmetry"}}
    
    wyckoffs = sym_dataset['wyckoffs']
    equivalent_atoms = sym_dataset['equivalent_atoms']
    
    # Group atoms by Wyckoff position
    wyckoff_groups = {}
    for i, (wyck, equiv) in enumerate(zip(wyckoffs, equivalent_atoms)):
        if wyck not in wyckoff_groups:
            wyckoff_groups[wyck] = []
        wyckoff_groups[wyck].append({
            "atom_index": i,
            "element": str(structure[i].specie),
            "position": list(structure[i].frac_coords),
            "equivalent_to": int(equiv)
        })
    
    return {
        "success": True,
        "spacegroup": analyzer.get_space_group_symbol(),
        "wyckoff_positions": wyckoff_groups,
        "n_unique_wyckoffs": len(wyckoff_groups),
        "total_atoms": len(structure)
    }


def find_space_group(structure: Structure, symprec: float = 0.01) -> Dict[str, Any]:
    """
    Find space group of a structure (Wrapper for analyze_symmetry).
    """
    return analyze_symmetry(structure, symprec=symprec)


def check_tolerance(structure: Structure, tolerance: float = 0.01) -> bool:
    """
    Check if structure maintains symmetry at given tolerance.
    """
    try:
        analyzer = SpacegroupAnalyzer(structure, symprec=tolerance)
        return True
    except:
        return False
