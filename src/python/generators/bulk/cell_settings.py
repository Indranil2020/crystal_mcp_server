"""
bulk/cell_settings.py - Cell Setting Transformations

Transforms structures between different cell conventions:
- Primitive cell (smallest unit cell)
- Conventional cell (standard crystallographic setting)
- Niggli reduced cell (unique reduced form)
- Standardized cell (ITA standard setting)

Scientific References:
- International Tables for Crystallography, Vol. A, 6th Ed. (2016)
- Niggli, P. "Krystallographische und strukturtheoretische Grundbegriffe" (1928)
- Grosse-Kunstleve et al. Acta Cryst. A60, 1-6 (2004)
"""

from typing import Dict, Any, List, Optional, Union, Tuple
import importlib.util
import numpy as np

Structure = None
Lattice = None
SpacegroupAnalyzer = None
PYMATGEN_AVAILABLE = False
if (
    importlib.util.find_spec("pymatgen.core") is not None
    and importlib.util.find_spec("pymatgen.symmetry.analyzer") is not None
):
    from pymatgen.core import Structure, Lattice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    PYMATGEN_AVAILABLE = True

SPGLIB_AVAILABLE = importlib.util.find_spec("spglib") is not None
if SPGLIB_AVAILABLE:
    import spglib


# Cell setting types
CELL_SETTINGS = {
    "primitive": {
        "description": "Smallest unit cell with one lattice point",
        "use_case": "DFT calculations (fewer atoms)"
    },
    "conventional": {
        "description": "Standard crystallographic cell (may have multiple lattice points)",
        "use_case": "Visualization, comparison with literature"
    },
    "niggli": {
        "description": "Unique reduced cell form (Niggli reduction)",
        "use_case": "Structure comparison, database matching"
    },
    "standardized": {
        "description": "ITA standard setting with origin choice",
        "use_case": "Crystallographic databases, publication"
    }
}


def _structure_to_spglib(structure: Dict[str, Any]) -> Tuple:
    """Convert structure dict to spglib format (cell, positions, numbers)."""
    if "structure" in structure:
        struct = structure["structure"]
    else:
        struct = structure

    lattice = struct.get("lattice", {})
    atoms = struct.get("atoms", struct.get("sites", []))

    # Get lattice matrix
    if "matrix" in lattice:
        cell = np.array(lattice["matrix"])
    else:
        a = lattice.get("a", 5.0)
        b = lattice.get("b", a)
        c = lattice.get("c", a)
        alpha = np.radians(lattice.get("alpha", 90))
        beta = np.radians(lattice.get("beta", 90))
        gamma = np.radians(lattice.get("gamma", 90))

        # Crystallographic to Cartesian
        cos_alpha = np.cos(alpha)
        cos_beta = np.cos(beta)
        cos_gamma = np.cos(gamma)
        sin_gamma = np.sin(gamma)

        val = (cos_alpha - cos_beta * cos_gamma) / sin_gamma

        cell = np.array([
            [a, 0, 0],
            [b * cos_gamma, b * sin_gamma, 0],
            [c * cos_beta, c * val, c * np.sqrt(1 - cos_beta**2 - val**2)]
        ])

    # Get positions and atomic numbers
    positions = []
    numbers = []

    # Complete element to atomic number map (Z=1 to Z=118)
    ELEMENT_Z = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
        "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
        "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
        "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
        "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
        "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
        "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99,
        "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
        "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111,
        "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
    }

    for atom in atoms:
        elem = atom.get("element", "X")
        coords = atom.get("coords", [0, 0, 0])
        positions.append(coords)
        numbers.append(ELEMENT_Z.get(elem, 0))

    return (cell, positions, numbers)


def _spglib_to_structure(cell: np.ndarray, positions: List, numbers: List) -> Dict[str, Any]:
    """Convert spglib format back to structure dict."""
    # Complete atomic number to element map (Z=1 to Z=118)
    Z_ELEMENT = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
        9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
        16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
        23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",
        30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
        37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
        44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
        51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La",
        58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
        65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
        72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt",
        79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
        86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
        93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es",
        100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db",
        106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg",
        112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
    }

    # Calculate lattice parameters from matrix
    a_vec = cell[0]
    b_vec = cell[1]
    c_vec = cell[2]

    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)

    alpha = np.degrees(np.arccos(np.dot(b_vec, c_vec) / (b * c)))
    beta = np.degrees(np.arccos(np.dot(a_vec, c_vec) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(a_vec, b_vec) / (a * b)))

    atoms = []
    for pos, z in zip(positions, numbers):
        atoms.append({
            "element": Z_ELEMENT.get(z, "X"),
            "coords": list(pos)
        })

    return {
        "lattice": {
            "a": float(a), "b": float(b), "c": float(c),
            "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
            "matrix": cell.tolist()
        },
        "atoms": atoms,
        "sites": atoms
    }


def get_primitive_cell(
    structure: Dict[str, Any],
    symprec: float = 1e-5
) -> Dict[str, Any]:
    """
    Convert structure to primitive cell.

    The primitive cell is the smallest unit cell that contains exactly
    one lattice point and can generate the entire crystal by translation.

    Args:
        structure: Input structure dictionary
        symprec: Symmetry tolerance for spglib

    Returns:
        Structure in primitive cell setting

    Reference:
        International Tables for Crystallography, Vol. A, Chapter 5
    """
    if not SPGLIB_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "SPGLIB_NOT_AVAILABLE",
                     "message": "spglib required for cell transformations"}
        }

    cell, positions, numbers = _structure_to_spglib(structure)

    # Find primitive cell
    primitive = spglib.find_primitive((cell, positions, numbers), symprec=symprec)

    if primitive is None:
        return {
            "success": False,
            "error": {"code": "PRIMITIVE_FAILED",
                     "message": "Could not find primitive cell"}
        }

    prim_cell, prim_pos, prim_numbers = primitive

    # Get symmetry info
    sym_dataset = spglib.get_symmetry_dataset((prim_cell, prim_pos, prim_numbers), symprec=symprec)

    result_structure = _spglib_to_structure(prim_cell, prim_pos, prim_numbers)

    return {
        "success": True,
        "cell_type": "primitive",
        "n_atoms_original": len(numbers),
        "n_atoms_primitive": len(prim_numbers),
        "reduction_factor": len(numbers) / len(prim_numbers) if len(prim_numbers) > 0 else 1,
        "spacegroup_number": sym_dataset["number"] if sym_dataset else None,
        "spacegroup_symbol": sym_dataset["international"] if sym_dataset else None,
        "structure": result_structure
    }


def get_conventional_cell(
    structure: Dict[str, Any],
    symprec: float = 1e-5
) -> Dict[str, Any]:
    """
    Convert structure to conventional cell.

    The conventional cell is the standard crystallographic unit cell
    that may contain multiple lattice points but follows crystallographic
    conventions for each crystal system.

    Args:
        structure: Input structure dictionary
        symprec: Symmetry tolerance

    Returns:
        Structure in conventional cell setting

    Reference:
        International Tables for Crystallography, Vol. A, Chapter 2
    """
    if not SPGLIB_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "SPGLIB_NOT_AVAILABLE",
                     "message": "spglib required for cell transformations"}
        }

    cell, positions, numbers = _structure_to_spglib(structure)

    # Get standardized conventional cell
    std = spglib.standardize_cell((cell, positions, numbers),
                                   to_primitive=False,
                                   no_idealize=False,
                                   symprec=symprec)

    if std is None:
        return {
            "success": False,
            "error": {"code": "CONVENTIONAL_FAILED",
                     "message": "Could not find conventional cell"}
        }

    conv_cell, conv_pos, conv_numbers = std

    # Get symmetry info
    sym_dataset = spglib.get_symmetry_dataset((conv_cell, conv_pos, conv_numbers), symprec=symprec)

    result_structure = _spglib_to_structure(conv_cell, conv_pos, conv_numbers)

    return {
        "success": True,
        "cell_type": "conventional",
        "n_atoms_original": len(numbers),
        "n_atoms_conventional": len(conv_numbers),
        "expansion_factor": len(conv_numbers) / len(numbers) if len(numbers) > 0 else 1,
        "spacegroup_number": sym_dataset["number"] if sym_dataset else None,
        "spacegroup_symbol": sym_dataset["international"] if sym_dataset else None,
        "hall_symbol": sym_dataset["hall"] if sym_dataset else None,
        "structure": result_structure
    }


def get_niggli_cell(
    structure: Dict[str, Any],
    eps: float = 1e-5
) -> Dict[str, Any]:
    """
    Convert structure to Niggli reduced cell.

    The Niggli cell is a unique reduced cell form that allows
    unambiguous comparison of crystal structures. It satisfies
    specific conditions on the metric tensor.

    Args:
        structure: Input structure dictionary
        eps: Tolerance for Niggli reduction

    Returns:
        Structure in Niggli reduced cell

    Reference:
        Niggli, P. (1928). Krystallographische und strukturtheoretische Grundbegriffe.
        Grosse-Kunstleve, R.W. et al. Acta Cryst. A60, 1-6 (2004).
    """
    if not SPGLIB_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "SPGLIB_NOT_AVAILABLE",
                     "message": "spglib required for cell transformations"}
        }

    cell, positions, numbers = _structure_to_spglib(structure)

    # Apply Niggli reduction
    niggli_cell = spglib.niggli_reduce(cell, eps=eps)

    if niggli_cell is None:
        return {
            "success": False,
            "error": {"code": "NIGGLI_FAILED",
                     "message": "Niggli reduction failed"}
        }

    # Transform atomic positions to new cell
    # Find transformation matrix
    inv_original = np.linalg.inv(cell)
    transform = niggli_cell @ inv_original

    # Transform fractional coordinates
    new_positions = []
    for pos in positions:
        cart = np.array(pos) @ cell
        new_frac = cart @ np.linalg.inv(niggli_cell)
        # Wrap to [0, 1)
        new_frac = new_frac % 1.0
        new_positions.append(new_frac.tolist())

    result_structure = _spglib_to_structure(niggli_cell, new_positions, numbers)

    # Calculate Niggli parameters for verification
    a_vec = niggli_cell[0]
    b_vec = niggli_cell[1]
    c_vec = niggli_cell[2]

    A = np.dot(a_vec, a_vec)
    B = np.dot(b_vec, b_vec)
    C = np.dot(c_vec, c_vec)
    xi = 2 * np.dot(b_vec, c_vec)
    eta = 2 * np.dot(a_vec, c_vec)
    zeta = 2 * np.dot(a_vec, b_vec)

    return {
        "success": True,
        "cell_type": "niggli",
        "n_atoms": len(numbers),
        "niggli_parameters": {
            "A": float(A), "B": float(B), "C": float(C),
            "xi": float(xi), "eta": float(eta), "zeta": float(zeta)
        },
        "is_type_I": xi * eta * zeta > 0,  # Type I or Type II Niggli cell
        "structure": result_structure
    }


def get_standardized_cell(
    structure: Dict[str, Any],
    to_primitive: bool = False,
    no_idealize: bool = False,
    symprec: float = 1e-5
) -> Dict[str, Any]:
    """
    Convert structure to ITA standard setting.

    Standardizes the cell according to International Tables for
    Crystallography conventions, including proper origin choice
    and axis orientation.

    Args:
        structure: Input structure dictionary
        to_primitive: If True, return primitive cell instead of conventional
        no_idealize: If True, don't idealize lattice parameters
        symprec: Symmetry tolerance

    Returns:
        Structure in standardized setting

    Reference:
        International Tables for Crystallography, Vol. A, 6th Ed. (2016)
    """
    if not SPGLIB_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "SPGLIB_NOT_AVAILABLE",
                     "message": "spglib required for cell transformations"}
        }

    cell, positions, numbers = _structure_to_spglib(structure)

    # Get symmetry dataset first
    sym_dataset = spglib.get_symmetry_dataset((cell, positions, numbers), symprec=symprec)

    if sym_dataset is None:
        return {
            "success": False,
            "error": {"code": "SYMMETRY_DETECTION_FAILED",
                     "message": "Could not detect symmetry"}
        }

    # Standardize
    std = spglib.standardize_cell((cell, positions, numbers),
                                   to_primitive=to_primitive,
                                   no_idealize=no_idealize,
                                   symprec=symprec)

    if std is None:
        return {
            "success": False,
            "error": {"code": "STANDARDIZATION_FAILED",
                     "message": "Cell standardization failed"}
        }

    std_cell, std_pos, std_numbers = std

    result_structure = _spglib_to_structure(std_cell, std_pos, std_numbers)

    # Determine crystal system
    sg_number = sym_dataset["number"]
    if sg_number <= 2:
        crystal_system = "triclinic"
    elif sg_number <= 15:
        crystal_system = "monoclinic"
    elif sg_number <= 74:
        crystal_system = "orthorhombic"
    elif sg_number <= 142:
        crystal_system = "tetragonal"
    elif sg_number <= 167:
        crystal_system = "trigonal"
    elif sg_number <= 194:
        crystal_system = "hexagonal"
    else:
        crystal_system = "cubic"

    return {
        "success": True,
        "cell_type": "standardized",
        "to_primitive": to_primitive,
        "idealized": not no_idealize,
        "n_atoms": len(std_numbers),
        "spacegroup_number": sym_dataset["number"],
        "spacegroup_symbol": sym_dataset["international"],
        "hall_symbol": sym_dataset["hall"],
        "hall_number": sym_dataset["hall_number"],
        "crystal_system": crystal_system,
        "pointgroup": sym_dataset["pointgroup"],
        "origin_shift": sym_dataset["origin_shift"].tolist(),
        "transformation_matrix": sym_dataset["transformation_matrix"].tolist(),
        "wyckoff_letters": list(sym_dataset["wyckoffs"]) if sym_dataset["wyckoffs"] is not None else [],
        "structure": result_structure
    }


def transform_cell(
    structure: Dict[str, Any],
    setting: str = "conventional",
    symprec: float = 1e-5
) -> Dict[str, Any]:
    """
    Transform structure to specified cell setting.

    Convenience function that wraps all cell transformation functions.

    Args:
        structure: Input structure dictionary
        setting: Target cell type: "primitive", "conventional", "niggli", "standardized"
        symprec: Symmetry tolerance

    Returns:
        Structure in requested cell setting
    """
    transformers = {
        "primitive": lambda s: get_primitive_cell(s, symprec=symprec),
        "conventional": lambda s: get_conventional_cell(s, symprec=symprec),
        "niggli": lambda s: get_niggli_cell(s, eps=symprec),
        "standardized": lambda s: get_standardized_cell(s, symprec=symprec)
    }

    if setting not in transformers:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SETTING",
                "message": f"Unknown cell setting: {setting}. Valid: {list(transformers.keys())}"
            }
        }

    return transformers[setting](structure)


def compare_cell_settings(
    structure: Dict[str, Any],
    symprec: float = 1e-5
) -> Dict[str, Any]:
    """
    Compare all cell settings for a structure.

    Useful for understanding the relationship between different
    cell representations of the same crystal structure.

    Args:
        structure: Input structure dictionary
        symprec: Symmetry tolerance

    Returns:
        Comparison of all cell settings
    """
    results = {}

    for setting in ["primitive", "conventional", "niggli", "standardized"]:
        result = transform_cell(structure, setting, symprec)
        if result["success"]:
            struct = result["structure"]
            lattice = struct["lattice"]
            results[setting] = {
                "n_atoms": len(struct.get("atoms", [])),
                "a": lattice["a"],
                "b": lattice["b"],
                "c": lattice["c"],
                "alpha": lattice["alpha"],
                "beta": lattice["beta"],
                "gamma": lattice["gamma"],
                "volume": abs(np.linalg.det(np.array(lattice["matrix"])))
            }
        else:
            results[setting] = {"error": result.get("error", {}).get("message", "Failed")}

    return {
        "success": True,
        "comparison": results,
        "description": CELL_SETTINGS
    }
