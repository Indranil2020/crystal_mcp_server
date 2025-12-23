"""
Iterative Structure Editing Module

Provides capabilities for iteratively editing structures without starting from scratch.
Supports:
- Chained operations with history tracking
- In-place modifications
- Undo/redo capability
- Structure comparison
- Edit validation

Example usage:
    # Create editor instance
    editor = StructureEditor(initial_structure)

    # Chain modifications
    editor.add_adsorbate("CO", site="ontop")
    editor.create_vacancy(site_index=5)
    editor.apply_strain(0.02, "biaxial")

    # Get result
    final_structure = editor.get_structure()
    history = editor.get_history()
"""

from typing import Dict, Any, List, Optional, Tuple, Union, Callable
from copy import deepcopy
import numpy as np
from dataclasses import dataclass, field
from datetime import datetime


@dataclass
class EditOperation:
    """Record of a single edit operation."""
    name: str
    params: Dict[str, Any]
    timestamp: str
    structure_before_hash: str
    structure_after_hash: str
    reversible: bool = True
    undo_params: Optional[Dict[str, Any]] = None


class StructureEditor:
    """
    Interactive structure editor with history tracking.

    Allows iterative modifications to structures with undo/redo capability
    and operation chaining.
    """

    def __init__(self, structure: Optional[Dict[str, Any]] = None, track_history: bool = True):
        """
        Initialize the editor.

        Parameters
        ----------
        structure : dict, optional
            Initial structure to edit
        track_history : bool
            Whether to track edit history (default: True)
        """
        self._structure = deepcopy(structure) if structure else None
        self._track_history = track_history
        self._history: List[Dict[str, Any]] = []
        self._undo_stack: List[Dict[str, Any]] = []
        self._redo_stack: List[Dict[str, Any]] = []

    def set_structure(self, structure: Dict[str, Any]) -> 'StructureEditor':
        """Set the current structure."""
        self._structure = deepcopy(structure)
        if self._track_history:
            self._history.append({
                "operation": "set_structure",
                "timestamp": datetime.now().isoformat(),
                "n_atoms": structure.get("n_atoms", len(structure.get("species", [])))
            })
        return self

    def get_structure(self) -> Dict[str, Any]:
        """Get the current structure."""
        return deepcopy(self._structure)

    def get_history(self) -> List[Dict[str, Any]]:
        """Get the edit history."""
        return deepcopy(self._history)

    def undo(self) -> bool:
        """Undo the last operation."""
        if not self._undo_stack:
            return False

        state = self._undo_stack.pop()
        self._redo_stack.append(deepcopy(self._structure))
        self._structure = state["structure"]

        if self._track_history:
            self._history.append({
                "operation": "undo",
                "undone_operation": state["operation"],
                "timestamp": datetime.now().isoformat()
            })

        return True

    def redo(self) -> bool:
        """Redo the last undone operation."""
        if not self._redo_stack:
            return False

        state = self._redo_stack.pop()
        self._undo_stack.append({
            "structure": deepcopy(self._structure),
            "operation": "redo"
        })
        self._structure = state

        if self._track_history:
            self._history.append({
                "operation": "redo",
                "timestamp": datetime.now().isoformat()
            })

        return True

    def _save_undo_state(self, operation: str):
        """Save current state for undo."""
        self._undo_stack.append({
            "structure": deepcopy(self._structure),
            "operation": operation
        })
        self._redo_stack.clear()

    # =========================================================================
    # Atom Manipulation Operations
    # =========================================================================

    def add_atom(
        self,
        element: str,
        position: List[float],
        validate: bool = True
    ) -> 'StructureEditor':
        """
        Add a single atom to the structure.

        Parameters
        ----------
        element : str
            Element symbol (e.g., "C", "O", "Fe")
        position : list
            [x, y, z] position in Angstroms
        validate : bool
            Check for overlaps (default: True)

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set. Use set_structure() first.")

        self._save_undo_state("add_atom")

        species = list(self._structure.get("species", []))
        positions = list(self._structure.get("positions", []))

        if validate:
            # Check for overlaps
            min_dist = self._get_min_distance(position, positions)
            if min_dist < 0.5:  # 0.5 Å minimum
                raise ValueError(f"Atom too close to existing atom (d={min_dist:.2f} Å)")

        species.append(element)
        positions.append(position)

        self._structure["species"] = species
        self._structure["positions"] = positions
        self._structure["n_atoms"] = len(species)

        if self._track_history:
            self._history.append({
                "operation": "add_atom",
                "element": element,
                "position": position,
                "timestamp": datetime.now().isoformat()
            })

        return self

    def remove_atom(self, index: int) -> 'StructureEditor':
        """
        Remove an atom by index.

        Parameters
        ----------
        index : int
            Index of atom to remove (0-based)

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("remove_atom")

        species = list(self._structure.get("species", []))
        positions = list(self._structure.get("positions", []))

        if index < 0 or index >= len(species):
            raise IndexError(f"Atom index {index} out of range")

        removed_element = species.pop(index)
        removed_pos = positions.pop(index)

        self._structure["species"] = species
        self._structure["positions"] = positions
        self._structure["n_atoms"] = len(species)

        if self._track_history:
            self._history.append({
                "operation": "remove_atom",
                "index": index,
                "removed_element": removed_element,
                "removed_position": removed_pos,
                "timestamp": datetime.now().isoformat()
            })

        return self

    def move_atom(self, index: int, new_position: List[float]) -> 'StructureEditor':
        """
        Move an atom to a new position.

        Parameters
        ----------
        index : int
            Index of atom to move
        new_position : list
            New [x, y, z] position

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("move_atom")

        positions = list(self._structure.get("positions", []))
        old_position = positions[index]
        positions[index] = new_position
        self._structure["positions"] = positions

        if self._track_history:
            self._history.append({
                "operation": "move_atom",
                "index": index,
                "old_position": old_position,
                "new_position": new_position,
                "timestamp": datetime.now().isoformat()
            })

        return self

    def substitute_atom(self, index: int, new_element: str) -> 'StructureEditor':
        """
        Substitute an atom with a different element.

        Parameters
        ----------
        index : int
            Index of atom to substitute
        new_element : str
            New element symbol

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("substitute_atom")

        species = list(self._structure.get("species", []))
        old_element = species[index]
        species[index] = new_element
        self._structure["species"] = species

        if self._track_history:
            self._history.append({
                "operation": "substitute_atom",
                "index": index,
                "old_element": old_element,
                "new_element": new_element,
                "timestamp": datetime.now().isoformat()
            })

        return self

    def translate_atoms(
        self,
        translation: List[float],
        indices: Optional[List[int]] = None
    ) -> 'StructureEditor':
        """
        Translate atoms by a vector.

        Parameters
        ----------
        translation : list
            [dx, dy, dz] translation vector
        indices : list, optional
            Specific atom indices to translate. If None, translates all.

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("translate_atoms")

        positions = np.array(self._structure.get("positions", []))
        translation = np.array(translation)

        if indices is None:
            positions += translation
        else:
            for i in indices:
                positions[i] += translation

        self._structure["positions"] = positions.tolist()

        if self._track_history:
            self._history.append({
                "operation": "translate_atoms",
                "translation": list(translation),
                "indices": indices,
                "timestamp": datetime.now().isoformat()
            })

        return self

    def rotate_atoms(
        self,
        angle: float,
        axis: str = "z",
        center: Optional[List[float]] = None,
        indices: Optional[List[int]] = None
    ) -> 'StructureEditor':
        """
        Rotate atoms around an axis.

        Parameters
        ----------
        angle : float
            Rotation angle in degrees
        axis : str
            Rotation axis ("x", "y", "z")
        center : list, optional
            Center of rotation. If None, uses center of mass.
        indices : list, optional
            Specific atom indices to rotate. If None, rotates all.

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("rotate_atoms")

        from scipy.spatial.transform import Rotation as R

        positions = np.array(self._structure.get("positions", []))

        # Determine center
        if center is None:
            if indices is None:
                center = positions.mean(axis=0)
            else:
                center = positions[indices].mean(axis=0)
        else:
            center = np.array(center)

        # Create rotation
        axis_map = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
        axis_vec = axis_map.get(axis.lower(), [0, 0, 1])
        rot = R.from_rotvec(np.radians(angle) * np.array(axis_vec))

        # Apply rotation
        if indices is None:
            positions = rot.apply(positions - center) + center
        else:
            for i in indices:
                positions[i] = rot.apply(positions[i] - center) + center

        self._structure["positions"] = positions.tolist()

        if self._track_history:
            self._history.append({
                "operation": "rotate_atoms",
                "angle": angle,
                "axis": axis,
                "center": list(center),
                "indices": indices,
                "timestamp": datetime.now().isoformat()
            })

        return self

    # =========================================================================
    # Defect Operations
    # =========================================================================

    def create_vacancy(self, site_index: int) -> 'StructureEditor':
        """
        Create a vacancy by removing an atom.

        Parameters
        ----------
        site_index : int
            Index of atom to remove

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        return self.remove_atom(site_index)

    def create_interstitial(
        self,
        element: str,
        position: Optional[List[float]] = None,
        interstitial_type: str = "octahedral"
    ) -> 'StructureEditor':
        """
        Create an interstitial defect.

        Parameters
        ----------
        element : str
            Element to insert
        position : list, optional
            Position for interstitial. If None, finds best site.
        interstitial_type : str
            Type of interstitial site ("octahedral", "tetrahedral")

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if position is None:
            position = self._find_interstitial_site(interstitial_type)

        return self.add_atom(element, position)

    def create_substitution(
        self,
        site_index: int,
        dopant: str
    ) -> 'StructureEditor':
        """
        Create a substitutional defect.

        Parameters
        ----------
        site_index : int
            Index of atom to substitute
        dopant : str
            Element to substitute with

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        return self.substitute_atom(site_index, dopant)

    # =========================================================================
    # Cell Operations
    # =========================================================================

    def scale_cell(self, factor: Union[float, List[float]]) -> 'StructureEditor':
        """
        Scale the unit cell.

        Parameters
        ----------
        factor : float or list
            Scaling factor(s). Single value for uniform scaling,
            [fx, fy, fz] for anisotropic.

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("scale_cell")

        cell = np.array(self._structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
        positions = np.array(self._structure.get("positions", []))

        if isinstance(factor, (int, float)):
            factor = [factor, factor, factor]

        factor = np.array(factor)

        # Scale cell
        new_cell = cell * factor.reshape(3, 1)

        # Scale positions proportionally
        if len(positions) > 0:
            # Convert to fractional, then back to Cartesian with new cell
            inv_cell = np.linalg.inv(cell)
            frac_positions = positions @ inv_cell
            new_positions = frac_positions @ new_cell
            self._structure["positions"] = new_positions.tolist()

        self._structure["cell"] = new_cell.tolist()

        if self._track_history:
            self._history.append({
                "operation": "scale_cell",
                "factor": list(factor),
                "timestamp": datetime.now().isoformat()
            })

        return self

    def apply_strain(
        self,
        strain: float,
        strain_type: str = "biaxial"
    ) -> 'StructureEditor':
        """
        Apply strain to the structure.

        Parameters
        ----------
        strain : float
            Strain magnitude (e.g., 0.02 for 2%)
        strain_type : str
            Type of strain ("biaxial", "uniaxial_x", "uniaxial_y", "uniaxial_z", "hydrostatic")

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if strain_type == "biaxial":
            factor = [1 + strain, 1 + strain, 1.0]
        elif strain_type == "uniaxial_x":
            factor = [1 + strain, 1.0, 1.0]
        elif strain_type == "uniaxial_y":
            factor = [1.0, 1 + strain, 1.0]
        elif strain_type == "uniaxial_z":
            factor = [1.0, 1.0, 1 + strain]
        elif strain_type == "hydrostatic":
            factor = [1 + strain, 1 + strain, 1 + strain]
        else:
            raise ValueError(f"Unknown strain type: {strain_type}")

        return self.scale_cell(factor)

    def make_supercell(self, scaling: Union[int, List[int]]) -> 'StructureEditor':
        """
        Create a supercell.

        Parameters
        ----------
        scaling : int or list
            Scaling factors. Single value or [nx, ny, nz].

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("make_supercell")

        if isinstance(scaling, int):
            scaling = [scaling, scaling, scaling]

        cell = np.array(self._structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
        species = self._structure.get("species", [])
        positions = np.array(self._structure.get("positions", []))

        new_species = []
        new_positions = []

        for i in range(scaling[0]):
            for j in range(scaling[1]):
                for k in range(scaling[2]):
                    translation = i * cell[0] + j * cell[1] + k * cell[2]
                    for sym, pos in zip(species, positions):
                        new_species.append(sym)
                        new_positions.append(pos + translation)

        new_cell = cell * np.array(scaling).reshape(3, 1)

        self._structure["species"] = new_species
        self._structure["positions"] = [list(p) for p in new_positions]
        self._structure["cell"] = new_cell.tolist()
        self._structure["n_atoms"] = len(new_species)

        if self._track_history:
            self._history.append({
                "operation": "make_supercell",
                "scaling": scaling,
                "new_n_atoms": len(new_species),
                "timestamp": datetime.now().isoformat()
            })

        return self

    # =========================================================================
    # Adsorbate Operations
    # =========================================================================

    def add_adsorbate(
        self,
        adsorbate: str,
        site: str = "ontop",
        site_index: Optional[int] = None,
        height: float = 2.0,
        rotation: float = 0.0
    ) -> 'StructureEditor':
        """
        Add an adsorbate molecule to the surface.

        Parameters
        ----------
        adsorbate : str
            Molecule name (e.g., "CO", "H2O", "NH3")
        site : str
            Adsorption site type ("ontop", "bridge", "hollow")
        site_index : int, optional
            Specific surface atom index
        height : float
            Height above surface in Angstroms
        rotation : float
            Rotation angle in degrees

        Returns
        -------
        StructureEditor
            Self for chaining
        """
        if self._structure is None:
            raise ValueError("No structure set")

        self._save_undo_state("add_adsorbate")

        from ase.build import molecule as ase_molecule

        # Get adsorbate
        try:
            mol = ase_molecule(adsorbate)
        except:
            mol = self._build_simple_molecule(adsorbate)

        # Find adsorption site
        positions = np.array(self._structure.get("positions", []))
        z_max = positions[:, 2].max()

        if site_index is not None:
            site_pos = positions[site_index, :2]
        else:
            site_pos = self._find_site_position(positions, site)

        # Position adsorbate
        mol_center = mol.get_center_of_mass()
        mol.translate([
            site_pos[0] - mol_center[0],
            site_pos[1] - mol_center[1],
            z_max + height - mol_center[2]
        ])

        if rotation != 0:
            from scipy.spatial.transform import Rotation as R
            rot = R.from_euler('z', rotation, degrees=True)
            mol.positions = rot.apply(mol.positions - mol.get_center_of_mass()) + mol.get_center_of_mass()

        # Add to structure
        species = list(self._structure.get("species", []))
        struct_positions = list(self._structure.get("positions", []))

        for sym, pos in zip(mol.get_chemical_symbols(), mol.positions):
            species.append(sym)
            struct_positions.append(list(pos))

        self._structure["species"] = species
        self._structure["positions"] = struct_positions
        self._structure["n_atoms"] = len(species)

        if self._track_history:
            self._history.append({
                "operation": "add_adsorbate",
                "adsorbate": adsorbate,
                "site": site,
                "height": height,
                "rotation": rotation,
                "timestamp": datetime.now().isoformat()
            })

        return self

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _get_min_distance(self, pos: List[float], positions: List[List[float]]) -> float:
        """Get minimum distance from pos to any position in list."""
        if not positions:
            return float('inf')
        pos = np.array(pos)
        positions = np.array(positions)
        distances = np.linalg.norm(positions - pos, axis=1)
        return distances.min()

    def _find_interstitial_site(self, site_type: str) -> List[float]:
        """Find an interstitial site position."""
        positions = np.array(self._structure.get("positions", []))
        cell = np.array(self._structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))

        # Simple approach: find void in structure
        center = positions.mean(axis=0)
        return list(center)

    def _find_site_position(self, positions: np.ndarray, site: str) -> np.ndarray:
        """Find adsorption site position."""
        z_max = positions[:, 2].max()
        top_mask = positions[:, 2] > z_max - 0.5
        top_positions = positions[top_mask]

        if site == "ontop":
            return top_positions[0, :2]
        elif site == "bridge":
            if len(top_positions) >= 2:
                return (top_positions[0, :2] + top_positions[1, :2]) / 2
        elif site == "hollow":
            if len(top_positions) >= 3:
                return np.mean(top_positions[:3, :2], axis=0)
        elif site == "center":
            cell = np.array(self._structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
            return np.array([cell[0, 0]/2, cell[1, 1]/2])

        return top_positions[0, :2]

    def _build_simple_molecule(self, name: str):
        """Build simple molecule by name."""
        from ase import Atoms

        molecules = {
            "CO": Atoms("CO", positions=[[0, 0, 0], [0, 0, 1.13]]),
            "H2O": Atoms("H2O", positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]]),
            "NH3": Atoms("NH4", positions=[[0, 0, 0], [0, 0.94, 0.38], [0.81, -0.47, 0.38], [-0.81, -0.47, 0.38]]),
        }

        return molecules.get(name, Atoms("H", positions=[[0, 0, 0]]))


# =============================================================================
# Functional API for single operations
# =============================================================================

def edit_structure(
    structure: Dict[str, Any],
    operations: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """
    Apply a sequence of edit operations to a structure.

    Parameters
    ----------
    structure : dict
        Input structure
    operations : list
        List of operation dicts, each with:
        - "operation": operation name
        - other params specific to the operation

    Returns
    -------
    dict
        Result with edited structure and history

    Example
    -------
    >>> result = edit_structure(
    ...     structure=si_structure,
    ...     operations=[
    ...         {"operation": "make_supercell", "scaling": [2, 2, 2]},
    ...         {"operation": "create_vacancy", "site_index": 0},
    ...         {"operation": "apply_strain", "strain": 0.02, "strain_type": "biaxial"}
    ...     ]
    ... )
    """
    editor = StructureEditor(structure, track_history=True)

    for op in operations:
        op_name = op.pop("operation")
        method = getattr(editor, op_name, None)

        if method is None:
            return {
                "success": False,
                "error": f"Unknown operation: {op_name}"
            }

        try:
            method(**op)
        except Exception as e:
            return {
                "success": False,
                "error": f"Operation '{op_name}' failed: {str(e)}",
                "partial_structure": editor.get_structure(),
                "completed_operations": editor.get_history()[:-1]
            }

    return {
        "success": True,
        "structure": editor.get_structure(),
        "history": editor.get_history(),
        "n_operations": len(operations)
    }


def compare_structures(
    structure1: Dict[str, Any],
    structure2: Dict[str, Any],
    tolerance: float = 0.01
) -> Dict[str, Any]:
    """
    Compare two structures and report differences.

    Parameters
    ----------
    structure1 : dict
        First structure
    structure2 : dict
        Second structure
    tolerance : float
        Tolerance for position comparison in Angstroms

    Returns
    -------
    dict
        Comparison results
    """
    s1 = structure1
    s2 = structure2

    species1 = s1.get("species", [])
    species2 = s2.get("species", [])
    pos1 = np.array(s1.get("positions", []))
    pos2 = np.array(s2.get("positions", []))
    cell1 = np.array(s1.get("cell", [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    cell2 = np.array(s2.get("cell", [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))

    # Compare composition
    from collections import Counter
    comp1 = Counter(species1)
    comp2 = Counter(species2)

    # Compare cell
    cell_diff = np.abs(cell1 - cell2).max()

    # Position comparison (if same composition)
    position_rmsd = None
    if len(pos1) == len(pos2) and comp1 == comp2:
        # Simple RMSD (doesn't account for atom reordering)
        position_rmsd = np.sqrt(np.mean((pos1 - pos2)**2))

    return {
        "success": True,
        "identical": comp1 == comp2 and (position_rmsd is not None and position_rmsd < tolerance),
        "n_atoms_1": len(species1),
        "n_atoms_2": len(species2),
        "composition_1": dict(comp1),
        "composition_2": dict(comp2),
        "compositions_match": comp1 == comp2,
        "cell_difference_max": float(cell_diff),
        "position_rmsd": float(position_rmsd) if position_rmsd is not None else None,
        "tolerance": tolerance
    }


def get_available_edit_operations() -> Dict[str, Any]:
    """Get list of available edit operations."""
    return {
        "atom_operations": {
            "add_atom": "Add a single atom",
            "remove_atom": "Remove an atom by index",
            "move_atom": "Move an atom to new position",
            "substitute_atom": "Replace atom with different element",
            "translate_atoms": "Translate atoms by vector",
            "rotate_atoms": "Rotate atoms around axis"
        },
        "defect_operations": {
            "create_vacancy": "Remove atom to create vacancy",
            "create_interstitial": "Add interstitial atom",
            "create_substitution": "Substitutional doping"
        },
        "cell_operations": {
            "scale_cell": "Scale unit cell",
            "apply_strain": "Apply strain tensor",
            "make_supercell": "Create supercell"
        },
        "adsorbate_operations": {
            "add_adsorbate": "Add molecule to surface"
        },
        "utilities": {
            "edit_structure": "Apply sequence of operations",
            "compare_structures": "Compare two structures"
        }
    }
