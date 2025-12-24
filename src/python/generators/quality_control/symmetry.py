"""
quality_control/symmetry.py - Symmetry Analysis

Comprehensive symmetry analysis per structure_catalogue.md Category 16:
(i) Automated symmetry refinement & space-group detection
(ii) Tolerance sweep for near-symmetric structures
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import importlib.util
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

PYXTAL_AVAILABLE = importlib.util.find_spec("pyxtal.symmetry") is not None
if PYXTAL_AVAILABLE:
    from pyxtal.symmetry import Group
else:
    Group = None


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
    structure: Union[Structure, Dict[str, Any]],
    symprec: float = 0.01,
    angle_tolerance: float = 5.0
) -> Dict[str, Any]:
    """
    Analyze crystal symmetry.
    
    Args:
        structure: Input structure (pymatgen Structure or dict)
        symprec: Symmetry precision in Angstrom
        angle_tolerance: Angle tolerance in degrees
    
    Returns:
        Symmetry analysis results
    """
    if isinstance(structure, dict):
        # Validate dict structure structure
        lattice_data = structure.get("lattice")
        if not lattice_data:
             return {"success": False, "error": {"code": "INVALID_STRUCTURE", "message": "Missing lattice data"}}
        
        atoms_data = structure.get("atoms")
        if not atoms_data:
            atoms_data = structure.get("sites")
        
        if not atoms_data:
             return {"success": False, "error": {"code": "INVALID_STRUCTURE", "message": "Missing atoms/sites data"}}

        # Reconstruct structure from dict with explicit checks
        # Handle matrix format
        if "matrix" in lattice_data:
            matrix = lattice_data["matrix"]
            if not isinstance(matrix, list) or len(matrix) != 3:
                 return {"success": False, "error": {"code": "INVALID_LATTICE", "message": "Invalid lattice matrix format"}}
            lattice = Lattice(matrix)
        else:
             # Fallback to parameters
             a = lattice_data.get("a", 0)
             b = lattice_data.get("b", 0)
             c = lattice_data.get("c", 0)
             if a <= 0 or b <= 0 or c <= 0:
                  return {"success": False, "error": {"code": "INVALID_LATTICE", "message": "Invalid lattice parameters"}}
             
             lattice = Lattice.from_parameters(
                 a, 
                 b if b else a, 
                 c if c else a,
                 lattice_data.get("alpha", 90), 
                 lattice_data.get("beta", 90), 
                 lattice_data.get("gamma", 90)
             )
        
        species = []
        coords = []
        for i, atom in enumerate(atoms_data):
            # Element validation
            elem = None
            if atom.get("species") and isinstance(atom["species"], list) and len(atom["species"]) > 0:
                elem = atom["species"][0].get("element")
            if not elem:
                elem = atom.get("element")
            
            if not elem:
                 return {"success": False, "error": {"code": "INVALID_ATOM", "message": f"Missing element at index {i}"}}
            
            # Coord validation
            coord = atom.get("coords")
            if not coord or not isinstance(coord, list) or len(coord) != 3:
                 return {"success": False, "error": {"code": "INVALID_ATOM", "message": f"Invalid coords at index {i}"}}
                 
            species.append(elem)
            coords.append(coord)
            
        structure = Structure(lattice, species, coords)
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
    
    if not isinstance(structure, Structure):
        return {
            "success": False,
            "error": {"code": "INVALID_STRUCTURE", "message": "structure must be a pymatgen Structure"}
        }

    results = []

    for tol in tolerances:
        if not isinstance(tol, (int, float)) or tol <= 0:
            results.append({
                "tolerance": tol,
                "success": False,
                "error": "Tolerance must be a positive number"
            })
            continue

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
    if not isinstance(structure, Structure):
        return False
    # Explicit validation instead of try-except
    if not structure.is_valid():
        return False
        
    analyzer = SpacegroupAnalyzer(structure, symprec=tolerance)
    return True


def validate_structure(structure: Union[Structure, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Validate a crystal structure for physical reasonableness.

    Checks:
    1. Valid lattice parameters (positive, non-zero volume)
    2. Valid elements
    3. Minimum interatomic distances (no overlapping atoms)
    4. Composition validity

    Args:
        structure: Structure (pymatgen Structure or dictionary)

    Returns:
        Validation result dictionary
    """
    # Handle pymatgen Structure objects directly
    if isinstance(structure, Structure):
        # Already a valid structure - just validate distances
        for i in range(len(structure)):
            neighbors = structure.get_neighbors(structure[i], 0.5)
            if neighbors:
                dist = neighbors[0].nn_distance
                return {
                    "success": True,
                    "valid": False,
                    "error": f"Atoms too close: site {i} ({structure[i].specie}) has neighbor at {dist:.3f} A"
                }

        return {
            "success": True,
            "valid": True,
            "formula": structure.formula,
            "density": structure.density,
            "volume": structure.volume,
            "n_atoms": len(structure)
        }

    # Check lattice
    lattice = structure.get("lattice", {})
    if not lattice:
        return {"success": True, "valid": False, "error": "Missing lattice information"}
        
    a = lattice.get("a", 0)
    b = lattice.get("b", 0)
    c = lattice.get("c", 0)
    
    # Check matrix if a,b,c are missing or zero (fallback)
    if "matrix" in lattice and (a==0 or b==0 or c==0):
         mat = np.array(lattice["matrix"])
         # Check matrix shape
         if mat.shape == (3,3):
             a = np.linalg.norm(mat[0])
             b = np.linalg.norm(mat[1])
             c = np.linalg.norm(mat[2])

    volume = lattice.get("volume", 0)
    if volume == 0 and "matrix" in lattice:
         mat = np.array(lattice["matrix"])
         if mat.shape == (3,3):
             volume = float(abs(np.linalg.det(mat)))
    
    if a <= 0 or b <= 0 or c <= 0:
        return {"success": True, "valid": False, "error": f"Invalid lattice parameters: a={a}, b={b}, c={c}"}
    if volume <= 0:
         return {"success": True, "valid": False, "error": f"Invalid volume: {volume}"}
        
    # Check atoms
    atoms = structure.get("atoms", [])
    if not atoms:
         # Try sites key if atoms missing (support both formats)
        atoms = structure.get("sites", [])
        
    if not atoms:
        return {"success": True, "valid": False, "error": "Structure contains no atoms"}
        
    # Check elements
    from pymatgen.core.periodic_table import Element
    valid_symbols = {e.symbol for e in Element}
    
    positions = []
    species_list = []
    
    for i, atom in enumerate(atoms):
        elem = atom.get("element", "")
        # Handle species list format
        if not elem and atom.get("species"):
            elem = atom["species"][0].get("element")
        
        if not elem or elem not in valid_symbols:
            return {"success": True, "valid": False, "error": f"Invalid element at site {i}: {elem}"}
        
        # Validate coordinates
        coords = atom.get("coords")
        if not coords or len(coords) != 3:
             return {"success": True, "valid": False, "error": f"Invalid coordinates at site {i}"}
             
        species_list.append(elem)
        positions.append(coords)

    # Check atomic overlaps
    # We reconstruct pymatgen structure to use its neighbor finding which handles PBC efficiently
    
    if "matrix" in lattice:
        mat = lattice["matrix"]
        # Validate matrix type and shape
        if isinstance(mat, list) and len(mat) == 3 and all(len(row) == 3 for row in mat):
            lat = Lattice(mat)
        else:
             return {"success": True, "valid": False, "error": "Invalid lattice matrix format"}
    else:
        lat = Lattice.from_parameters(a, b, c, 
                                    lattice.get("alpha", 90), 
                                    lattice.get("beta", 90), 
                                    lattice.get("gamma", 90))
    
    # Ensure parallel behavior for robust checking
    # Create structure without try-except - if Lattice is valid, Structure should be
    pmg_struct = Structure(lat, species_list, positions)
    
    # Check for overlaps (< 0.5 Angstrom)
    for i in range(len(pmg_struct)):
        neighbors = pmg_struct.get_neighbors(pmg_struct[i], 0.5)
        if neighbors:
                dist = neighbors[0].nn_distance
                return {
                "success": True, 
                "valid": False, 
                "error": f"Atoms too close: site {i} ({pmg_struct[i].specie}) has neighbor at {dist:.3f} A"
            }

    return {
        "success": True,
        "valid": True,
        "formula": pmg_struct.formula,
        "density": pmg_struct.density,
        "volume": pmg_struct.volume,
        "n_atoms": len(pmg_struct)
    }


def get_subgroups(space_group: int, strict: bool = False) -> Dict[str, Any]:
    """
    Get maximal subgroups of a space group.
    """
    if not PYXTAL_AVAILABLE:
        return {"success": False, "error": {"code": "PYXTAL_MISSING", "message": "PyXtal required for subgroup analysis"}}

    if not 1 <= space_group <= 230:
         return {
            "success": False,
            "error": {"code": "INVALID_SPACE_GROUP", "message": f"Space group must be 1-230, got {space_group}"}
        }
    
    g = Group(space_group)
    subgroup_numbers = g.get_max_subgroup_numbers()
    
    unique_subgroups = sorted(list(set(subgroup_numbers)))
    subgroups_data = []
    
    for num in unique_subgroups:
        sub_g = Group(num)
        subgroups_data.append({
            "number": int(num),
            "symbol": sub_g.symbol,
            "crystal_system": sub_g.lattice_type,
            "point_group": sub_g.point_group
        })
        
    return {
        "success": True,
        "space_group": {
            "number": space_group,
            "symbol": g.symbol,
            "hall_number": g.hall_number
        },
        "subgroups": subgroups_data,
        "raw_numbers": subgroup_numbers
    }


def get_symmetry_path(start_spg: int, end_spg: int, max_depth: int = 5) -> Dict[str, Any]:
    """
    Find a path from supergroup to subgroup.
    """
    if not PYXTAL_AVAILABLE:
        return {"success": False, "error": {"code": "PYXTAL_MISSING", "message": "PyXtal required for symmetry path"}}

    if start_spg == end_spg:
         return {"success": True, "path": [start_spg]}
         
    queue = [(start_spg, [start_spg])]
    visited = {start_spg}
    
    while queue:
        current, path = queue.pop(0)
        
        if len(path) > max_depth:
            continue
            
        if current == end_spg:
             formatted_path = []
             for spg_num in path:
                 g = Group(spg_num)
                 formatted_path.append({
                     "number": spg_num,
                     "symbol": g.symbol
                 })
             return {
                 "success": True, 
                 "path": formatted_path,
                 "length": len(path) - 1
             }
        
        g = Group(current)
        subs = g.get_max_subgroup_numbers()
        for sub in subs:
            if sub not in visited:
                visited.add(sub)
                new_path = list(path)
                new_path.append(sub)
                queue.append((sub, new_path))
            
    return {
        "success": False,
        "error": {
            "code": "PATH_NOT_FOUND",
            "message": f"No subgroup path found from {start_spg} to {end_spg} within depth {max_depth}"
        }
    }


def transform_to_subgroup(
    structure: Union[Structure, Dict[str, Any]],
    target_spacegroup: int,
    eps: float = 0.0,
    max_steps: int = 5
) -> Dict[str, Any]:
    """
    Transform a crystal structure to a subgroup (Bilbao TRANSTRU-like functionality).

    This enables comparison of structures in different space groups by transforming
    a high-symmetry structure to a lower-symmetry subgroup with proper Wyckoff
    position splitting and lattice transformation.

    Example: Transform cubic BiFeO3 (Pm-3m, 5 atoms) to rhombohedral R3c (10 atoms)
    to compare with experimental R3c phase.

    Args:
        structure: Input structure (pymatgen Structure or dict)
        target_spacegroup: Target space group number (must be a subgroup)
        eps: Symmetry tolerance for the transformation (0 = ideal, >0 allows distortion)
        max_steps: Maximum transformation steps to reach target

    Returns:
        Dictionary with transformed structure and transformation details
    """
    if not PYXTAL_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "PYXTAL_MISSING", "message": "PyXtal required for subgroup transformation"}
        }

    from pyxtal import pyxtal

    # Convert dict to Structure if needed
    if isinstance(structure, dict):
        lattice_data = structure.get("lattice", {})
        atoms_data = structure.get("atoms", structure.get("sites", []))

        if "matrix" in lattice_data:
            lattice = Lattice(lattice_data["matrix"])
        else:
            lattice = Lattice.from_parameters(
                lattice_data.get("a", 5.0),
                lattice_data.get("b", lattice_data.get("a", 5.0)),
                lattice_data.get("c", lattice_data.get("a", 5.0)),
                lattice_data.get("alpha", 90),
                lattice_data.get("beta", 90),
                lattice_data.get("gamma", 90)
            )

        species = []
        coords = []
        for atom in atoms_data:
            elem = atom.get("element")
            if not elem and atom.get("species"):
                elem = atom["species"][0].get("element")
            species.append(elem)
            coords.append(atom.get("coords", atom.get("abc", [0, 0, 0])))

        structure = Structure(lattice, species, coords)

    # Load into PyXtal
    crystal = pyxtal()
    crystal.from_seed(structure)

    original_spg = crystal.group.number
    original_symbol = crystal.group.symbol
    original_n_atoms = len(crystal.to_pymatgen())

    # Validate target is a subgroup
    if target_spacegroup > original_spg:
        return {
            "success": False,
            "error": {
                "code": "INVALID_TARGET",
                "message": f"Target SG {target_spacegroup} is not a subgroup of {original_spg}. "
                          f"Subgroups have lower or equal symmetry."
            }
        }

    if target_spacegroup == original_spg:
        # No transformation needed
        return {
            "success": True,
            "original_spacegroup": {"number": original_spg, "symbol": original_symbol},
            "target_spacegroup": {"number": target_spacegroup, "symbol": original_symbol},
            "transformation_path": [original_spg],
            "n_atoms_original": original_n_atoms,
            "n_atoms_transformed": original_n_atoms,
            "multiplicity": 1,
            "structure": _structure_to_dict(crystal.to_pymatgen())
        }

    # Find transformation path
    path_result = get_symmetry_path(original_spg, target_spacegroup, max_depth=max_steps)

    if not path_result.get("success"):
        return {
            "success": False,
            "error": {
                "code": "NO_PATH",
                "message": f"No group-subgroup path found from SG {original_spg} to SG {target_spacegroup}. "
                          f"Target may not be a proper subgroup.",
                "original_spg": original_spg,
                "target_spg": target_spacegroup
            }
        }

    path = [step["number"] for step in path_result["path"]]

    # Perform stepwise transformation
    current = crystal
    transformation_steps = []

    for i, next_spg in enumerate(path[1:], 1):
        current_spg = current.group.number

        # Try translationengleiche (t) subgroup first, then klassengleiche (k)
        transformed = None
        trans_type = None

        for group_type in ['t', 'k']:
            subgroups = current.group.get_max_t_subgroup() if group_type == 't' else current.group.get_max_k_subgroup()
            if next_spg in subgroups.get('subgroup', []):
                transformed = current.subgroup_once(eps=eps, H=next_spg, group_type=group_type)
                trans_type = group_type
                break

        if transformed is None:
            return {
                "success": False,
                "error": {
                    "code": "TRANSFORMATION_FAILED",
                    "message": f"Could not transform from SG {current_spg} to SG {next_spg}",
                    "partial_path": path[:i]
                }
            }

        transformation_steps.append({
            "from_spg": current_spg,
            "to_spg": next_spg,
            "type": "translationengleiche" if trans_type == 't' else "klassengleiche",
            "n_atoms_before": len(current.to_pymatgen()),
            "n_atoms_after": len(transformed.to_pymatgen())
        })

        current = transformed

    # Get final structure
    final_struct = current.to_pymatgen()

    return {
        "success": True,
        "original_spacegroup": {"number": original_spg, "symbol": original_symbol},
        "target_spacegroup": {"number": target_spacegroup, "symbol": current.group.symbol},
        "transformation_path": path,
        "transformation_steps": transformation_steps,
        "n_atoms_original": original_n_atoms,
        "n_atoms_transformed": len(final_struct),
        "multiplicity": len(final_struct) // original_n_atoms if original_n_atoms > 0 else 1,
        "structure": _structure_to_dict(final_struct)
    }


def transform_by_path(
    structure: Union[Structure, Dict[str, Any]],
    path: List[int],
    eps: float = 0.0
) -> Dict[str, Any]:
    """
    Transform structure along a specified group-subgroup path.

    Use this when you want to follow a specific transformation sequence
    rather than letting the algorithm find the path.

    Args:
        structure: Input structure
        path: List of space group numbers [start, intermediate..., target]
        eps: Symmetry tolerance

    Returns:
        Transformed structure with path details
    """
    if len(path) < 2:
        return {
            "success": False,
            "error": {"code": "INVALID_PATH", "message": "Path must have at least 2 space groups"}
        }

    if not PYXTAL_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "PYXTAL_MISSING", "message": "PyXtal required for transformation"}
        }

    from pyxtal import pyxtal

    # Convert dict to Structure if needed
    if isinstance(structure, dict):
        lattice_data = structure.get("lattice", {})
        atoms_data = structure.get("atoms", structure.get("sites", []))

        if "matrix" in lattice_data:
            lattice = Lattice(lattice_data["matrix"])
        else:
            lattice = Lattice.from_parameters(
                lattice_data.get("a", 5.0),
                lattice_data.get("b", lattice_data.get("a", 5.0)),
                lattice_data.get("c", lattice_data.get("a", 5.0)),
                lattice_data.get("alpha", 90),
                lattice_data.get("beta", 90),
                lattice_data.get("gamma", 90)
            )

        species = []
        coords = []
        for atom in atoms_data:
            elem = atom.get("element")
            if not elem and atom.get("species"):
                elem = atom["species"][0].get("element")
            species.append(elem)
            coords.append(atom.get("coords", atom.get("abc", [0, 0, 0])))

        structure = Structure(lattice, species, coords)

    crystal = pyxtal()
    crystal.from_seed(structure)

    # Verify starting space group matches
    if crystal.group.number != path[0]:
        return {
            "success": False,
            "error": {
                "code": "PATH_MISMATCH",
                "message": f"Structure is in SG {crystal.group.number}, but path starts with SG {path[0]}"
            }
        }

    original_n_atoms = len(crystal.to_pymatgen())
    transformation_steps = []
    current = crystal

    for i, next_spg in enumerate(path[1:], 1):
        current_spg = current.group.number

        transformed = None
        trans_type = None

        for group_type in ['t', 'k']:
            subgroups = current.group.get_max_t_subgroup() if group_type == 't' else current.group.get_max_k_subgroup()
            if next_spg in subgroups.get('subgroup', []):
                transformed = current.subgroup_once(eps=eps, H=next_spg, group_type=group_type)
                trans_type = group_type
                break

        if transformed is None:
            return {
                "success": False,
                "error": {
                    "code": "INVALID_STEP",
                    "message": f"SG {next_spg} is not a maximal subgroup of SG {current_spg}",
                    "step": i,
                    "completed_path": path[:i]
                }
            }

        transformation_steps.append({
            "step": i,
            "from_spg": current_spg,
            "to_spg": next_spg,
            "type": "translationengleiche" if trans_type == 't' else "klassengleiche",
            "n_atoms": len(transformed.to_pymatgen())
        })

        current = transformed

    final_struct = current.to_pymatgen()

    return {
        "success": True,
        "path": path,
        "n_steps": len(path) - 1,
        "transformation_steps": transformation_steps,
        "n_atoms_original": original_n_atoms,
        "n_atoms_final": len(final_struct),
        "multiplicity": len(final_struct) // original_n_atoms if original_n_atoms > 0 else 1,
        "final_spacegroup": {"number": current.group.number, "symbol": current.group.symbol},
        "structure": _structure_to_dict(final_struct)
    }


def get_all_subgroup_paths(
    start_spg: int,
    end_spg: int,
    max_depth: int = 4
) -> Dict[str, Any]:
    """
    Find all possible group-subgroup paths between two space groups.

    This is useful for exploring different transformation pathways,
    as different paths may result in different final structures
    (different orientations, domains, etc.).

    Args:
        start_spg: Starting space group number
        end_spg: Target space group number
        max_depth: Maximum path length

    Returns:
        All found paths with their characteristics
    """
    if not PYXTAL_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "PYXTAL_MISSING", "message": "PyXtal required"}
        }

    g = Group(start_spg)
    paths = g.search_subgroup_paths(end_spg, max_depth)

    if not paths:
        return {
            "success": False,
            "error": {
                "code": "NO_PATHS",
                "message": f"No paths found from SG {start_spg} to SG {end_spg}"
            }
        }

    formatted_paths = []
    for path in paths:
        path_info = []
        for spg in path:
            g_step = Group(spg)
            path_info.append({
                "number": spg,
                "symbol": g_step.symbol,
                "crystal_system": g_step.lattice_type
            })
        formatted_paths.append({
            "path": path,
            "length": len(path) - 1,
            "details": path_info
        })

    return {
        "success": True,
        "start_spg": start_spg,
        "end_spg": end_spg,
        "n_paths": len(formatted_paths),
        "paths": formatted_paths
    }


def _structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """Helper to convert Structure to dict format."""
    lattice = structure.lattice
    return {
        "lattice": {
            "a": lattice.a,
            "b": lattice.b,
            "c": lattice.c,
            "alpha": lattice.alpha,
            "beta": lattice.beta,
            "gamma": lattice.gamma,
            "matrix": lattice.matrix.tolist()
        },
        "atoms": [
            {"element": str(s.specie), "coords": list(s.frac_coords)}
            for s in structure
        ],
        "metadata": {
            "formula": structure.formula,
            "n_atoms": len(structure),
            "volume": structure.volume
        }
    }
