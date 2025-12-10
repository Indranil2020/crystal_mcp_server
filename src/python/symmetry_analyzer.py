#!/usr/bin/env python3
"""
Symmetry Analyzer - Spglib Wrapper for Symmetry Operations

This module provides symmetry analysis capabilities using Spglib.
All functions follow defensive programming without try/except blocks.
"""

from typing import Dict, List, Optional, Any, Tuple
import json
import sys
import numpy as np
from pymatgen.core import Structure, Lattice
import spglib


def validate_structure_for_spglib(structure_dict: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """
    Validate that a structure dictionary has the required fields for spglib.

    Args:
        structure_dict: Structure dictionary to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not structure_dict or not isinstance(structure_dict, dict):
        return False, "Structure must be a dictionary"

    if "lattice" not in structure_dict:
        return False, "Structure missing 'lattice' field"

    if "atoms" not in structure_dict:
        return False, "Structure missing 'atoms' field"

    lattice = structure_dict["lattice"]
    if "matrix" not in lattice:
        return False, "Lattice missing 'matrix' field"

    atoms = structure_dict["atoms"]
    if not atoms or len(atoms) == 0:
        return False, "Structure has no atoms"

    for atom in atoms:
        if "element" not in atom:
            return False, "Atom missing 'element' field"
        if "coords" not in atom:
            return False, "Atom missing 'coords' field"

    return True, None


def dict_to_spglib_cell(structure_dict: Dict[str, Any]) -> Optional[Tuple]:
    """
    Convert structure dictionary to Spglib cell tuple.
    
    Args:
        structure_dict: Structure dictionary
    
    Returns:
        Spglib cell tuple (lattice, positions, numbers) or None
    """
    if not structure_dict or not isinstance(structure_dict, dict):
        return None
    
    if "lattice" not in structure_dict or "atoms" not in structure_dict:
        return None
    
    lattice = structure_dict["lattice"]
    if "matrix" not in lattice:
        return None
    
    lattice_matrix = np.array(lattice["matrix"])
    if lattice_matrix.shape != (3, 3):
        return None
    
    atoms_data = structure_dict["atoms"]
    if not atoms_data or len(atoms_data) == 0:
        return None
    
    positions = []
    numbers = []
    
    # Element to atomic number mapping
    element_to_z = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
        'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Ti': 22, 'V': 23,
        'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36
    }
    
    for atom in atoms_data:
        if "element" not in atom or "coords" not in atom:
            return None
        
        element = atom["element"]
        if element not in element_to_z:
            # Use dummy number for unknown elements
            z = 1
        else:
            z = element_to_z[element]
        
        positions.append(atom["coords"])
        numbers.append(z)
    
    positions_array = np.array(positions)
    numbers_array = np.array(numbers)
    
    return (lattice_matrix, positions_array, numbers_array)


def analyze_symmetry(
    structure_dict: Dict[str, Any],
    symprec: float = 1e-3,
    angle_tolerance: float = 5.0,
    detect_primitive: bool = False,
    standardize: bool = False
) -> Dict[str, Any]:
    """
    Analyze crystal symmetry using Spglib.
    
    Args:
        structure_dict: Structure to analyze
        symprec: Symmetry precision in Angstroms
        angle_tolerance: Angle tolerance in degrees
        detect_primitive: Find primitive cell
        standardize: Standardize to conventional cell
    
    Returns:
        Dictionary with symmetry analysis results
    """
    # Validate inputs
    if symprec <= 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"symprec must be positive, got {symprec}",
                "details": {"symprec": symprec}
            }
        }
    
    if angle_tolerance <= 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"angle_tolerance must be positive, got {angle_tolerance}",
                "details": {"angle_tolerance": angle_tolerance}
            }
        }
    
    # Convert to Spglib cell
    cell = dict_to_spglib_cell(structure_dict)
    if cell is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure to Spglib cell format",
                "details": {}
            }
        }
    
    # Get space group
    spacegroup = spglib.get_spacegroup(cell, symprec=symprec, angle_tolerance=angle_tolerance)
    if spacegroup is None:
        return {
            "success": False,
            "error": {
                "code": "SYMMETRY_DETECTION_FAILED",
                "message": "Failed to detect space group",
                "details": {"symprec": symprec}
            }
        }
    
    # Parse space group string (format: "Fd-3m (227)")
    spg_parts = spacegroup.split()
    symbol = spg_parts[0] if len(spg_parts) > 0 else "Unknown"
    number_str = spg_parts[1] if len(spg_parts) > 1 else "(1)"
    number = int(number_str.strip("()"))
    
    # Get symmetry dataset
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=angle_tolerance)
    if dataset is None:
        return {
            "success": False,
            "error": {
                "code": "SYMMETRY_ANALYSIS_FAILED",
                "message": "Failed to get symmetry dataset",
                "details": {}
            }
        }
    
    # Extract space group information
    space_group_info = {
        "number": int(dataset['number']),
        "symbol": dataset['international'],
        "hall_symbol": dataset['hall'],
        "point_group": dataset.get('pointgroup', 'Unknown'),
        "crystal_system": get_crystal_system(int(dataset['number']))
    }
    
    # Extract symmetry operations
    rotations = dataset['rotations'].tolist()
    translations = dataset['translations'].tolist()
    
    symmetry_operations = {
        "rotations": rotations,
        "translations": translations,
        "n_operations": len(rotations)
    }
    
    # Extract Wyckoff positions
    wyckoff_letters = dataset['wyckoffs']
    equivalent_atoms = dataset['equivalent_atoms']
    
    wyckoff_positions = []
    unique_atoms = {}
    
    lattice_matrix, positions, numbers = cell
    
    for i, (letter, equiv_idx) in enumerate(zip(wyckoff_letters, equivalent_atoms)):
        if equiv_idx not in unique_atoms:
            wyckoff_positions.append({
                "element": structure_dict["atoms"][i]["element"],
                "wyckoff": letter,
                "multiplicity": int(np.sum(equivalent_atoms == equiv_idx)),
                "site_symmetry": "Unknown",  # Spglib doesn't provide this directly
                "coords": positions[i].tolist()
            })
            unique_atoms[equiv_idx] = True
    
    result = {
        "success": True,
        "space_group": space_group_info,
        "symmetry_operations": symmetry_operations,
        "wyckoff_positions": wyckoff_positions,
        "is_standardized": False
    }
    
    # Get primitive cell if requested
    if detect_primitive:
        primitive_cell = spglib.find_primitive(cell, symprec=symprec, angle_tolerance=angle_tolerance)
        if primitive_cell is not None:
            prim_dict = spglib_cell_to_dict(primitive_cell)
            if prim_dict is not None:
                result["primitive_cell"] = prim_dict
    
    # Get conventional cell if requested
    if standardize:
        conventional_cell = spglib.standardize_cell(cell, to_primitive=False, 
                                                     no_idealize=False, symprec=symprec,
                                                     angle_tolerance=angle_tolerance)
        if conventional_cell is not None:
            conv_dict = spglib_cell_to_dict(conventional_cell)
            if conv_dict is not None:
                result["conventional_cell"] = conv_dict
                result["is_standardized"] = True
    
    return result


def spglib_cell_to_dict(cell: Tuple) -> Optional[Dict[str, Any]]:
    """
    Convert Spglib cell to structure dictionary.
    
    Args:
        cell: Spglib cell tuple (lattice, positions, numbers)
    
    Returns:
        Structure dictionary or None
    """
    if cell is None or len(cell) != 3:
        return None
    
    lattice_matrix, positions, numbers = cell
    
    # Calculate lattice parameters
    a = float(np.linalg.norm(lattice_matrix[0]))
    b = float(np.linalg.norm(lattice_matrix[1]))
    c = float(np.linalg.norm(lattice_matrix[2]))
    
    alpha = float(np.arccos(np.dot(lattice_matrix[1], lattice_matrix[2]) / (b * c)) * 180 / np.pi)
    beta = float(np.arccos(np.dot(lattice_matrix[0], lattice_matrix[2]) / (a * c)) * 180 / np.pi)
    gamma = float(np.arccos(np.dot(lattice_matrix[0], lattice_matrix[1]) / (a * b)) * 180 / np.pi)
    
    volume = float(np.abs(np.dot(lattice_matrix[0], np.cross(lattice_matrix[1], lattice_matrix[2]))))
    
    # Atomic number to element mapping (simplified)
    z_to_element = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 22: 'Ti', 23: 'V',
        24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
        31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr'
    }
    
    atoms = []
    for pos, z in zip(positions, numbers):
        element = z_to_element.get(int(z), f"X{z}")
        cartesian = np.dot(pos, lattice_matrix)
        
        atoms.append({
            "element": element,
            "coords": pos.tolist(),
            "cartesian": cartesian.tolist()
        })
    
    return {
        "lattice": {
            "a": a,
            "b": b,
            "c": c,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "matrix": lattice_matrix.tolist(),
            "volume": volume
        },
        "atoms": atoms,
        "metadata": {
            "natoms": len(atoms)
        }
    }


def get_crystal_system(space_group_number: int) -> str:
    """
    Get crystal system from space group number.
    
    Args:
        space_group_number: Space group number (1-230)
    
    Returns:
        Crystal system name
    """
    if 1 <= space_group_number <= 2:
        return "triclinic"
    elif 3 <= space_group_number <= 15:
        return "monoclinic"
    elif 16 <= space_group_number <= 74:
        return "orthorhombic"
    elif 75 <= space_group_number <= 142:
        return "tetragonal"
    elif 143 <= space_group_number <= 167:
        return "trigonal"
    elif 168 <= space_group_number <= 194:
        return "hexagonal"
    elif 195 <= space_group_number <= 230:
        return "cubic"
    else:
        return "unknown"


def main():
    """
    Main entry point for command-line execution.
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python symmetry_analyzer.py <input_json_file>"
            }
        }))
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input file
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    # Analyze symmetry
    result = analyze_symmetry(**params)
    
    # Output result
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
