"""
Workflow Schema Adapter

Converts between workflow editing schema and core structure schema.

Workflow schema:
    {
        "species": ["Si", "Si", "O", "O"],
        "positions": [[x, y, z], ...],  # Cartesian coordinates
        "cell": [[ax, ay, az], [bx, by, bz], [cx, cy, cz]],  # 3x3 matrix
        "n_atoms": 4
    }

Core schema:
    {
        "atoms": [
            {"element": "Si", "coords": [fx, fy, fz]},  # Fractional coordinates
            ...
        ],
        "lattice": {
            "matrix": [[ax, ay, az], [bx, by, bz], [cx, cy, cz]],
            "a": float, "b": float, "c": float,
            "alpha": float, "beta": float, "gamma": float,
            "volume": float
        },
        "metadata": {...},
        "space_group": {...}
    }
"""

from typing import Dict, Any, List
import numpy as np


def workflow_to_core(workflow_structure: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert workflow schema to core schema.

    Args:
        workflow_structure: Structure in workflow format (species/positions/cell)

    Returns:
        Structure in core format (atoms/lattice)
    """
    species = workflow_structure.get("species", [])
    positions_cart = np.array(workflow_structure.get("positions", []))
    cell = np.array(workflow_structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))

    # Calculate lattice parameters from cell matrix
    lattice_params = _cell_matrix_to_params(cell)

    # Convert Cartesian to fractional coordinates
    if len(positions_cart) > 0:
        cell_matrix = np.array(cell)
        if np.linalg.det(cell_matrix) == 0:
            raise ValueError("Cell matrix is singular (determinant = 0)")

        inv_cell = np.linalg.inv(cell_matrix)
        positions_frac = positions_cart @ inv_cell
    else:
        positions_frac = []

    # Build atoms list
    atoms = []
    for element, frac_coords in zip(species, positions_frac):
        atoms.append({
            "element": element,
            "coords": frac_coords.tolist() if isinstance(frac_coords, np.ndarray) else list(frac_coords)
        })

    # Build core structure
    core_structure = {
        "atoms": atoms,
        "lattice": {
            "matrix": cell.tolist() if isinstance(cell, np.ndarray) else cell,
            "a": lattice_params["a"],
            "b": lattice_params["b"],
            "c": lattice_params["c"],
            "alpha": lattice_params["alpha"],
            "beta": lattice_params["beta"],
            "gamma": lattice_params["gamma"],
            "volume": lattice_params["volume"]
        },
        "metadata": {
            "formula": _get_formula(species),
            "n_atoms": len(atoms),
            "source": "workflow_editing"
        },
        "space_group": {
            "number": 1,
            "symbol": "P1"
        }
    }

    return core_structure


def core_to_workflow(core_structure: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert core schema to workflow schema.

    Args:
        core_structure: Structure in core format (atoms/lattice)

    Returns:
        Structure in workflow format (species/positions/cell)
    """
    atoms = core_structure.get("atoms", [])
    lattice = core_structure.get("lattice", {})
    cell = np.array(lattice.get("matrix", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))

    # Extract species and fractional coords
    species = []
    positions_frac = []

    for atom in atoms:
        species.append(atom.get("element", "X"))
        coords = atom.get("coords", [0, 0, 0])
        positions_frac.append(coords)

    # Convert fractional to Cartesian coordinates
    if len(positions_frac) > 0:
        positions_frac = np.array(positions_frac)
        positions_cart = positions_frac @ cell
    else:
        positions_cart = []

    # Build workflow structure
    workflow_structure = {
        "species": species,
        "positions": positions_cart.tolist() if isinstance(positions_cart, np.ndarray) else positions_cart,
        "cell": cell.tolist() if isinstance(cell, np.ndarray) else cell,
        "n_atoms": len(species)
    }

    return workflow_structure


def _cell_matrix_to_params(cell: np.ndarray) -> Dict[str, float]:
    """
    Calculate lattice parameters from cell matrix.

    Cell matrix format:
        [[ax, ay, az],
         [bx, by, bz],
         [cx, cy, cz]]

    Returns:
        Dictionary with a, b, c, alpha, beta, gamma, volume
    """
    if cell.shape != (3, 3):
        raise ValueError(f"Cell matrix must be 3x3, got {cell.shape}")

    # Lattice vector lengths
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])

    # Angles (in degrees)
    # alpha: angle between b and c
    # beta: angle between a and c
    # gamma: angle between a and b

    def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
        """Calculate angle between two vectors in degrees."""
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        # Clamp to avoid numerical issues with arccos
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))

    alpha = angle_between(cell[1], cell[2])
    beta = angle_between(cell[0], cell[2])
    gamma = angle_between(cell[0], cell[1])

    # Volume (absolute value of determinant)
    volume = abs(np.linalg.det(cell))

    return {
        "a": float(a),
        "b": float(b),
        "c": float(c),
        "alpha": float(alpha),
        "beta": float(beta),
        "gamma": float(gamma),
        "volume": float(volume)
    }


def _get_formula(species: List[str]) -> str:
    """
    Generate chemical formula from species list.

    Args:
        species: List of element symbols

    Returns:
        Chemical formula string (e.g., "Si2O4")
    """
    if not species:
        return ""

    from collections import Counter
    composition = Counter(species)

    # Sort by element symbol
    formula_parts = []
    for element in sorted(composition.keys()):
        count = composition[element]
        if count == 1:
            formula_parts.append(element)
        else:
            formula_parts.append(f"{element}{count}")

    return "".join(formula_parts)


def validate_workflow_structure(structure: Dict[str, Any]) -> tuple[bool, str]:
    """
    Validate workflow structure format.

    Args:
        structure: Structure to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(structure, dict):
        return False, "Structure must be a dictionary"

    # Check required fields
    if "species" not in structure:
        return False, "Missing required field: species"

    if "positions" not in structure:
        return False, "Missing required field: positions"

    if "cell" not in structure:
        return False, "Missing required field: cell"

    species = structure["species"]
    positions = structure["positions"]
    cell = structure["cell"]

    # Validate types
    if not isinstance(species, list):
        return False, "Field 'species' must be a list"

    if not isinstance(positions, list):
        return False, "Field 'positions' must be a list"

    if not isinstance(cell, list):
        return False, "Field 'cell' must be a list"

    # Validate lengths
    if len(species) != len(positions):
        return False, f"Length mismatch: {len(species)} species but {len(positions)} positions"

    # Validate cell matrix
    if len(cell) != 3:
        return False, f"Cell must have 3 rows, got {len(cell)}"

    for i, row in enumerate(cell):
        if not isinstance(row, list) or len(row) != 3:
            return False, f"Cell row {i} must be a list of 3 numbers"

    # Check cell is not singular
    cell_matrix = np.array(cell)
    if abs(np.linalg.det(cell_matrix)) < 1e-10:
        return False, "Cell matrix is singular (determinant near zero)"

    return True, ""


def validate_core_structure(structure: Dict[str, Any]) -> tuple[bool, str]:
    """
    Validate core structure format.

    Args:
        structure: Structure to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(structure, dict):
        return False, "Structure must be a dictionary"

    # Check required fields
    if "atoms" not in structure:
        return False, "Missing required field: atoms"

    if "lattice" not in structure:
        return False, "Missing required field: lattice"

    atoms = structure["atoms"]
    lattice = structure["lattice"]

    # Validate atoms
    if not isinstance(atoms, list):
        return False, "Field 'atoms' must be a list"

    for i, atom in enumerate(atoms):
        if not isinstance(atom, dict):
            return False, f"Atom {i} must be a dictionary"

        if "element" not in atom:
            return False, f"Atom {i} missing required field: element"

        if "coords" not in atom:
            return False, f"Atom {i} missing required field: coords"

        coords = atom["coords"]
        if not isinstance(coords, list) or len(coords) != 3:
            return False, f"Atom {i} coords must be a list of 3 numbers"

    # Validate lattice
    if not isinstance(lattice, dict):
        return False, "Field 'lattice' must be a dictionary"

    if "matrix" not in lattice:
        return False, "Lattice missing required field: matrix"

    matrix = lattice["matrix"]
    if not isinstance(matrix, list) or len(matrix) != 3:
        return False, "Lattice matrix must be 3x3"

    for i, row in enumerate(matrix):
        if not isinstance(row, list) or len(row) != 3:
            return False, f"Lattice matrix row {i} must be a list of 3 numbers"

    return True, ""
