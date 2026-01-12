"""
Molecular Editor Module

Provides the main editing engine that orchestrates molecular modifications:
- Applies edit operations to molecules
- Manages RWMol state
- Handles coordinate regeneration
- Coordinates with validation pipeline
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union
import logging
import importlib.util
import copy

# Check module availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors

from .annotated_molecule import AnnotatedMolecule
from .edit_operations import (
    EditOperation,
    EditOperationType,
    EditResult,
    execute_operation,
    resolve_target,
)

logger = logging.getLogger(__name__)


@dataclass
class EditorState:
    """Tracks the state of the molecular editor."""
    original_smiles: str
    current_smiles: str
    operations_applied: List[EditResult] = field(default_factory=list)
    is_modified: bool = False
    is_valid: bool = True
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


class MolecularEditor:
    """
    RWMol-based molecular editing engine.

    This class manages the editing of molecules through structured operations.
    It maintains the RWMol state and coordinates with the AnnotatedMolecule
    for target resolution and annotation preservation.

    Usage:
        # Create from AnnotatedMolecule
        editor = MolecularEditor(annotated_mol)

        # Apply operations
        result = editor.apply_operation(operation)

        # Or apply multiple operations
        results = editor.apply_operations([op1, op2, op3])

        # Finalize and get edited molecule
        edited_mol = editor.finalize(optimize=True)
    """

    def __init__(self, molecule: AnnotatedMolecule):
        """
        Initialize editor with an AnnotatedMolecule.

        Args:
            molecule: The molecule to edit
        """
        self.original_molecule = molecule
        self.state = EditorState(
            original_smiles=molecule.smiles,
            current_smiles=molecule.smiles,
        )

        # Get RWMol for editing
        self._rwmol = molecule.get_rwmol()
        if self._rwmol is None and RDKIT_AVAILABLE:
            # Try to create from SMILES
            mol = Chem.MolFromSmiles(molecule.smiles)
            if mol is not None:
                self._rwmol = Chem.RWMol(mol)

        # Working copy of annotations
        self._annotations = copy.deepcopy({
            "stereocenters": [s.to_dict() for s in molecule.stereocenters],
            "rings": [r.to_dict() for r in molecule.rings],
            "functional_groups": [f.to_dict() for f in molecule.functional_groups],
            "atom_labels": molecule.atom_labels.copy(),
            "edit_history": molecule.edit_history.copy(),
        })

    @property
    def rwmol(self):
        """Get the RWMol object."""
        return self._rwmol

    @property
    def is_valid(self) -> bool:
        """Check if current molecule state is valid."""
        if self._rwmol is None:
            return False
        if not RDKIT_AVAILABLE:
            return False

        # Check if molecule can be sanitized
        mol_copy = Chem.RWMol(self._rwmol)
        result = Chem.SanitizeMol(mol_copy, catchErrors=True)
        return result == Chem.SanitizeFlags.SANITIZE_NONE

    def get_current_smiles(self) -> str:
        """Get current SMILES string."""
        if self._rwmol is None or not RDKIT_AVAILABLE:
            return self.state.current_smiles
        return Chem.MolToSmiles(self._rwmol)

    def apply_operation(self, operation: EditOperation) -> EditResult:
        """
        Apply a single edit operation to the molecule.

        Args:
            operation: The operation to apply

        Returns:
            EditResult with success/failure information
        """
        if self._rwmol is None:
            return EditResult(
                success=False,
                operation_type=operation.operation_type,
                target=operation.target,
                message="No valid molecule to edit",
                error_code="NO_MOLECULE",
            )

        # Create working copy for annotations lookup
        working_mol = AnnotatedMolecule.from_dict({
            **self.original_molecule.to_dict(),
            "smiles": self.get_current_smiles(),
        })
        if working_mol is not None:
            working_mol._mol = Chem.Mol(self._rwmol)

        # Execute the operation
        result = execute_operation(self._rwmol, operation, working_mol)

        # Record the operation
        self.state.operations_applied.append(result)

        if result.success:
            self.state.is_modified = True
            self.state.current_smiles = self.get_current_smiles()

            # Update edit history
            self._annotations["edit_history"].append({
                "operation": operation.operation_type.value,
                "target": operation.target.to_dict(),
                "params": operation.params,
                "success": True,
                "message": result.message,
            })

        else:
            self.state.errors.append(result.message)
            self._annotations["edit_history"].append({
                "operation": operation.operation_type.value,
                "target": operation.target.to_dict(),
                "params": operation.params,
                "success": False,
                "error": result.error_code,
                "message": result.message,
            })

        # Collect warnings
        self.state.warnings.extend(result.warnings)

        return result

    def apply_operations(self, operations: List[EditOperation]) -> List[EditResult]:
        """
        Apply multiple edit operations in sequence.

        Operations are applied in order. If an operation fails, subsequent
        operations are still attempted unless stop_on_error is set.

        Args:
            operations: List of operations to apply

        Returns:
            List of EditResult objects
        """
        results = []
        for operation in operations:
            result = self.apply_operation(operation)
            results.append(result)

            # Log each operation
            if result.success:
                logger.debug(f"Applied {operation.operation_type.value}: {result.message}")
            else:
                logger.warning(f"Failed {operation.operation_type.value}: {result.message}")

        return results

    def finalize(
        self,
        optimize: bool = True,
        add_hydrogens: bool = True,
        generate_3d: bool = True
    ) -> Optional[AnnotatedMolecule]:
        """
        Finalize editing and return the modified AnnotatedMolecule.

        This method:
        1. Adds explicit hydrogens if requested
        2. Generates 3D coordinates if requested
        3. Optimizes geometry if requested
        4. Creates new AnnotatedMolecule with updated annotations

        Args:
            optimize: Whether to optimize geometry with force field
            add_hydrogens: Whether to add explicit hydrogens
            generate_3d: Whether to generate 3D coordinates

        Returns:
            AnnotatedMolecule or None if finalization fails
        """
        if self._rwmol is None:
            logger.error("No valid molecule to finalize")
            return None

        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available for finalization")
            return None

        # Get final mol
        mol = self._rwmol.GetMol()

        # Sanitize
        sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
        if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
            logger.error(f"Sanitization failed during finalization: {sanitize_result}")
            self.state.is_valid = False
            self.state.errors.append("Sanitization failed")
            # Continue anyway to return what we have

        # Add hydrogens
        if add_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        coordinates = None
        has_3d = False

        if generate_3d:
            embed_result = self._generate_3d_coordinates(mol, optimize)
            if embed_result is not None:
                mol, coordinates = embed_result
                has_3d = True

        # Create new AnnotatedMolecule
        smiles = Chem.MolToSmiles(mol, canonical=True)
        canonical = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)

        result = AnnotatedMolecule(
            smiles=smiles,
            canonical_smiles=canonical,
            formula=formula,
            molecular_weight=mw,
            coordinates=coordinates,
            has_3d=has_3d,
            edit_history=self._annotations["edit_history"],
            _mol=mol,
        )

        # Re-annotate
        result._annotate_atoms()
        result._annotate_stereocenters()
        result._annotate_double_bond_stereo()
        result._annotate_rings()
        result._annotate_functional_groups()
        result._generate_atom_labels()

        return result

    def _generate_3d_coordinates(
        self,
        mol,
        optimize: bool = True
    ) -> Optional[Tuple[Any, List[Tuple[float, float, float]]]]:
        """
        Generate 3D coordinates for the molecule.

        Returns:
            Tuple of (mol with conformer, coordinates list) or None on failure
        """
        if not RDKIT_AVAILABLE:
            return None

        # Check if ETKDGv3 is available
        etkdg_params = None
        if hasattr(AllChem, 'ETKDGv3'):
            etkdg_params = AllChem.ETKDGv3()
        elif hasattr(AllChem, 'ETKDG'):
            etkdg_params = AllChem.ETKDG()

        if etkdg_params is not None:
            etkdg_params.randomSeed = 42

        # Embed molecule
        embed_result = AllChem.EmbedMolecule(mol, etkdg_params) if etkdg_params else AllChem.EmbedMolecule(mol)

        if embed_result == -1:
            # Try with random coordinates
            embed_result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

        if embed_result == -1:
            logger.warning("Failed to generate 3D coordinates")
            return None

        # Optimize if requested
        if optimize:
            # Try MMFF first
            mmff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            if mmff_result != 0:
                # Fall back to UFF
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        # Extract coordinates
        conf = mol.GetConformer()
        coordinates = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coordinates.append((pos.x, pos.y, pos.z))

        return (mol, coordinates)

    def rollback(self) -> bool:
        """
        Rollback to the original molecule state.

        Returns:
            True if rollback succeeded
        """
        if not RDKIT_AVAILABLE:
            return False

        mol = Chem.MolFromSmiles(self.original_molecule.smiles)
        if mol is None:
            return False

        self._rwmol = Chem.RWMol(mol)
        self.state = EditorState(
            original_smiles=self.original_molecule.smiles,
            current_smiles=self.original_molecule.smiles,
        )
        self._annotations = copy.deepcopy({
            "stereocenters": [s.to_dict() for s in self.original_molecule.stereocenters],
            "rings": [r.to_dict() for r in self.original_molecule.rings],
            "functional_groups": [f.to_dict() for f in self.original_molecule.functional_groups],
            "atom_labels": self.original_molecule.atom_labels.copy(),
            "edit_history": self.original_molecule.edit_history.copy(),
        })

        return True

    def get_state(self) -> EditorState:
        """Get current editor state."""
        return self.state

    def get_operation_history(self) -> List[Dict[str, Any]]:
        """Get list of applied operations."""
        return [r.to_dict() for r in self.state.operations_applied]


def edit_molecule(
    smiles: str,
    operations: List[Union[EditOperation, Dict[str, Any]]],
    optimize: bool = True,
    generate_3d: bool = True
) -> Tuple[Optional[AnnotatedMolecule], List[EditResult]]:
    """
    Convenience function to edit a molecule from SMILES.

    Args:
        smiles: SMILES string of molecule to edit
        operations: List of EditOperation objects or operation dicts
        optimize: Whether to optimize geometry
        generate_3d: Whether to generate 3D coordinates

    Returns:
        Tuple of (edited AnnotatedMolecule, list of EditResults)
        Returns (None, results) if molecule creation fails
    """
    # Create AnnotatedMolecule
    mol = AnnotatedMolecule.from_smiles(smiles)
    if mol is None:
        return (None, [EditResult(
            success=False,
            operation_type=EditOperationType.CUSTOM,
            target=None,
            message=f"Failed to parse SMILES: {smiles}",
            error_code="INVALID_SMILES",
        )])

    # Create editor
    editor = MolecularEditor(mol)

    # Parse operations if needed
    parsed_ops = []
    for op in operations:
        if isinstance(op, EditOperation):
            parsed_ops.append(op)
        elif isinstance(op, dict):
            parsed_ops.append(EditOperation.from_dict(op))

    # Apply operations
    results = editor.apply_operations(parsed_ops)

    # Finalize
    edited = editor.finalize(optimize=optimize, generate_3d=generate_3d)

    return (edited, results)


def create_editor_from_smiles(smiles: str) -> Optional[MolecularEditor]:
    """
    Create a MolecularEditor from a SMILES string.

    Args:
        smiles: SMILES string

    Returns:
        MolecularEditor or None if SMILES is invalid
    """
    mol = AnnotatedMolecule.from_smiles(smiles)
    if mol is None:
        return None
    return MolecularEditor(mol)
