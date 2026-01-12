"""
Edit Operations Module

Defines structured edit operations for molecular modification:
- Operation types (hydroxylate, acetylate, epimerize, etc.)
- Target specifications (by label, SMARTS, atom index, ring)
- Operation executors
- Result objects
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union, Callable
from enum import Enum
import re
import logging
import importlib.util

# Check module availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

logger = logging.getLogger(__name__)


class EditOperationType(Enum):
    """All supported edit operation types."""
    # Functional group additions
    HYDROXYLATE = "hydroxylate"
    ACETYLATE = "acetylate"
    METHYLATE = "methylate"
    HALOGENATE = "halogenate"
    AMINATE = "aminate"
    CARBOXYLATE = "carboxylate"
    PHOSPHORYLATE = "phosphorylate"
    SULFONATE = "sulfonate"
    NITRATE = "nitrate"

    # Functional group removals
    DEHYDROXYLATE = "dehydroxylate"
    DEACETYLATE = "deacetylate"
    DEMETHYLATE = "demethylate"
    DEHALOGENATE = "dehalogenate"

    # Stereochemistry
    EPIMERIZE = "epimerize"
    INVERT_STEREO = "invert_stereo"
    SET_STEREO = "set_stereo"
    SET_E_Z = "set_e_z"

    # Atom operations
    ADD_ATOM = "add_atom"
    DELETE_ATOM = "delete_atom"
    CHANGE_ELEMENT = "change_element"
    SET_CHARGE = "set_charge"

    # Bond operations
    ADD_BOND = "add_bond"
    DELETE_BOND = "delete_bond"
    CHANGE_BOND_ORDER = "change_bond_order"

    # Ring operations
    ADD_RING = "add_ring"
    DELETE_RING = "delete_ring"
    FUSE_RING = "fuse_ring"
    OPEN_RING = "open_ring"
    EXPAND_RING = "expand_ring"
    CONTRACT_RING = "contract_ring"

    # Substitution
    SUBSTITUTE = "substitute"
    REPLACE_GROUP = "replace_group"

    # Protection
    PROTECT = "protect"
    DEPROTECT = "deprotect"

    # Oxidation/Reduction
    OXIDIZE = "oxidize"
    REDUCE = "reduce"

    # Custom
    CUSTOM = "custom"


class TargetType(Enum):
    """Types of target specifications."""
    ATOM_INDEX = "atom_idx"
    ATOM_LABEL = "label"
    SMARTS = "smarts"
    RING_INDEX = "ring_idx"
    RING_TYPE = "ring_type"
    FUNCTIONAL_GROUP = "group"
    ATOM_PAIR = "atom_pair"
    ALL = "all"


@dataclass
class TargetSpecification:
    """
    Specification for targeting atoms/bonds/groups in a molecule.

    Supports multiple targeting methods:
    - By atom index: target_type=ATOM_INDEX, value=5
    - By label: target_type=ATOM_LABEL, value="C-7"
    - By SMARTS: target_type=SMARTS, value="[OH]"
    - By ring: target_type=RING_INDEX, value=0 or RING_TYPE, value="oxetane"
    - By functional group: target_type=FUNCTIONAL_GROUP, value="hydroxyl:0"
    """
    target_type: TargetType
    value: Union[int, str, Tuple[int, int]]
    match_index: int = 0  # Which match to use if multiple (0 = first)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target_type": self.target_type.value,
            "value": self.value if not isinstance(self.value, tuple) else list(self.value),
            "match_index": self.match_index,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "TargetSpecification":
        target_type_value = data.get("target_type", "atom_idx")
        target_type = TargetType(target_type_value) if target_type_value in [e.value for e in TargetType] else TargetType.ATOM_INDEX
        value = data.get("value", 0)
        if isinstance(value, list):
            value = tuple(value)
        return cls(
            target_type=target_type,
            value=value,
            match_index=data.get("match_index", 0),
        )

    @classmethod
    def from_string(cls, target_str: str) -> "TargetSpecification":
        """
        Parse target specification from string format.

        Examples:
            "label:C-7" -> TargetSpecification(ATOM_LABEL, "C-7")
            "smarts:[OH]" -> TargetSpecification(SMARTS, "[OH]")
            "atom_idx:5" -> TargetSpecification(ATOM_INDEX, 5)
            "ring:oxetane:0" -> TargetSpecification(RING_TYPE, "oxetane", 0)
            "group:hydroxyl:1" -> TargetSpecification(FUNCTIONAL_GROUP, "hydroxyl", 1)
        """
        if ":" not in target_str:
            # Default to SMARTS if no prefix
            return cls(target_type=TargetType.SMARTS, value=target_str)

        parts = target_str.split(":", maxsplit=2)
        prefix = parts[0].lower()

        if prefix == "label":
            return cls(target_type=TargetType.ATOM_LABEL, value=parts[1])
        elif prefix == "smarts":
            return cls(target_type=TargetType.SMARTS, value=parts[1])
        elif prefix in ("atom_idx", "idx", "atom"):
            return cls(target_type=TargetType.ATOM_INDEX, value=int(parts[1]))
        elif prefix == "ring":
            match_idx = int(parts[2]) if len(parts) > 2 else 0
            # Check if value is numeric (ring index) or string (ring type)
            if parts[1].isdigit():
                return cls(target_type=TargetType.RING_INDEX, value=int(parts[1]), match_index=match_idx)
            return cls(target_type=TargetType.RING_TYPE, value=parts[1], match_index=match_idx)
        elif prefix == "group":
            match_idx = int(parts[2]) if len(parts) > 2 else 0
            return cls(target_type=TargetType.FUNCTIONAL_GROUP, value=parts[1], match_index=match_idx)
        elif prefix in ("pair", "bond"):
            # Format: pair:5,7 or bond:5,7
            atoms = parts[1].split(",")
            return cls(target_type=TargetType.ATOM_PAIR, value=(int(atoms[0]), int(atoms[1])))
        else:
            # Default to SMARTS
            return cls(target_type=TargetType.SMARTS, value=target_str)


@dataclass
class EditResult:
    """Result of an edit operation."""
    success: bool
    operation_type: EditOperationType
    target: TargetSpecification
    message: str = ""
    affected_atoms: List[int] = field(default_factory=list)
    new_atoms: List[int] = field(default_factory=list)
    removed_atoms: List[int] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    error_code: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "operation_type": self.operation_type.value,
            "target": self.target.to_dict(),
            "message": self.message,
            "affected_atoms": self.affected_atoms,
            "new_atoms": self.new_atoms,
            "removed_atoms": self.removed_atoms,
            "warnings": self.warnings,
            "error_code": self.error_code,
        }


@dataclass
class EditOperation:
    """
    A structured molecular edit operation.

    Contains:
    - Operation type (what to do)
    - Target specification (where to apply)
    - Parameters (how to do it)
    """
    operation_type: EditOperationType
    target: TargetSpecification
    params: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.operation_type.value,
            "target": self.target.to_dict() if isinstance(self.target, TargetSpecification) else self.target,
            "params": self.params,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EditOperation":
        """Parse from dictionary (as received from MCP)."""
        # Parse operation type
        op_type_str = data.get("type", "custom")
        op_type = EditOperationType.CUSTOM
        for e in EditOperationType:
            if e.value == op_type_str.lower():
                op_type = e
                break

        # Parse target
        target_data = data.get("target", {})
        if isinstance(target_data, str):
            target = TargetSpecification.from_string(target_data)
        elif isinstance(target_data, dict):
            target = TargetSpecification.from_dict(target_data)
        else:
            target = TargetSpecification(TargetType.ATOM_INDEX, 0)

        return cls(
            operation_type=op_type,
            target=target,
            params=data.get("params", {}),
        )


def create_operation(
    op_type: Union[str, EditOperationType],
    target: Union[str, TargetSpecification],
    **params
) -> EditOperation:
    """
    Convenience function to create an EditOperation.

    Args:
        op_type: Operation type (string or enum)
        target: Target specification (string like "label:C-7" or TargetSpecification)
        **params: Operation-specific parameters

    Returns:
        EditOperation instance
    """
    # Parse operation type
    if isinstance(op_type, str):
        parsed_type = EditOperationType.CUSTOM
        for e in EditOperationType:
            if e.value == op_type.lower():
                parsed_type = e
                break
    else:
        parsed_type = op_type

    # Parse target
    if isinstance(target, str):
        parsed_target = TargetSpecification.from_string(target)
    else:
        parsed_target = target

    return EditOperation(
        operation_type=parsed_type,
        target=parsed_target,
        params=params,
    )


# ============================================================================
# Operation Execution Functions
# ============================================================================
# These functions implement the actual molecular modifications using RDKit.
# Each returns an EditResult object.

def resolve_target(mol, target: TargetSpecification, annotated_mol=None) -> List[int]:
    """
    Resolve a target specification to a list of atom indices.

    Args:
        mol: RDKit molecule object
        target: Target specification
        annotated_mol: Optional AnnotatedMolecule for label/group resolution

    Returns:
        List of atom indices matching the target
    """
    if not RDKIT_AVAILABLE or mol is None:
        return []

    if target.target_type == TargetType.ATOM_INDEX:
        idx = int(target.value)
        if 0 <= idx < mol.GetNumAtoms():
            return [idx]
        return []

    elif target.target_type == TargetType.ATOM_LABEL:
        if annotated_mol is not None:
            idx = annotated_mol.get_atom_by_label(str(target.value))
            if idx is not None:
                return [idx]
        return []

    elif target.target_type == TargetType.SMARTS:
        pattern = Chem.MolFromSmarts(str(target.value))
        if pattern is None:
            return []
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) > target.match_index:
            return list(matches[target.match_index])
        elif len(matches) > 0:
            return list(matches[0])
        return []

    elif target.target_type == TargetType.RING_TYPE:
        if annotated_mol is not None:
            from .annotated_molecule import RingType
            rings = annotated_mol.get_ring_by_type(str(target.value))
            if len(rings) > target.match_index:
                return rings[target.match_index].atom_indices
        return []

    elif target.target_type == TargetType.RING_INDEX:
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        idx = int(target.value)
        if 0 <= idx < len(rings):
            return list(rings[idx])
        return []

    elif target.target_type == TargetType.FUNCTIONAL_GROUP:
        if annotated_mol is not None:
            from .annotated_molecule import FunctionalGroupType
            groups = annotated_mol.get_functional_group_by_type(str(target.value))
            if len(groups) > target.match_index:
                return groups[target.match_index].atom_indices
        return []

    elif target.target_type == TargetType.ATOM_PAIR:
        pair = target.value
        if isinstance(pair, (list, tuple)) and len(pair) == 2:
            return list(pair)
        return []

    elif target.target_type == TargetType.ALL:
        return list(range(mol.GetNumAtoms()))

    return []


def execute_hydroxylate(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Add hydroxyl group (-OH) to target atom.

    Params:
        orientation: "axial" or "equatorial" (for ring atoms, affects 3D placement)
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HYDROXYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, target_atoms[0] if target_atoms else 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HYDROXYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    # Validate target atom can accept a substituent
    atom = rwmol.GetAtomWithIdx(target_idx)
    if atom is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HYDROXYLATE,
            target=target_spec,
            message=f"Atom index {target_idx} not found",
            error_code="ATOM_NOT_FOUND",
        )

    # Check if atom has available valence
    explicit_valence = atom.GetExplicitValence()
    default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
    if isinstance(default_valence, tuple):
        default_valence = max(default_valence)

    # For carbon, check if we can add an -OH (need to replace an H)
    num_hs = atom.GetTotalNumHs()
    if num_hs < 1:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HYDROXYLATE,
            target=target_spec,
            message=f"Atom {target_idx} has no hydrogen to replace",
            error_code="NO_HYDROGEN",
        )

    # Add oxygen atom
    o_idx = rwmol.AddAtom(Chem.Atom(8))  # Oxygen

    # Add bond between target and oxygen
    rwmol.AddBond(target_idx, o_idx, Chem.BondType.SINGLE)

    # Update hydrogen counts
    atom.SetNumExplicitHs(num_hs - 1)
    oxygen = rwmol.GetAtomWithIdx(o_idx)
    oxygen.SetNumExplicitHs(1)

    # Sanitize to ensure valid molecule
    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HYDROXYLATE,
            target=target_spec,
            message=f"Sanitization failed after hydroxylation",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.HYDROXYLATE,
        target=target_spec,
        message=f"Added hydroxyl group to atom {target_idx}",
        affected_atoms=[target_idx],
        new_atoms=[o_idx],
    )


def execute_methylate(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """Add methyl group (-CH3) to target atom."""
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.METHYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.METHYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    atom = rwmol.GetAtomWithIdx(target_idx)
    if atom is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.METHYLATE,
            target=target_spec,
            message=f"Atom index {target_idx} not found",
            error_code="ATOM_NOT_FOUND",
        )

    num_hs = atom.GetTotalNumHs()
    if num_hs < 1:
        return EditResult(
            success=False,
            operation_type=EditOperationType.METHYLATE,
            target=target_spec,
            message=f"Atom {target_idx} has no hydrogen to replace",
            error_code="NO_HYDROGEN",
        )

    # Add carbon (methyl group)
    c_idx = rwmol.AddAtom(Chem.Atom(6))  # Carbon
    rwmol.AddBond(target_idx, c_idx, Chem.BondType.SINGLE)

    # Update hydrogen counts
    atom.SetNumExplicitHs(num_hs - 1)
    methyl = rwmol.GetAtomWithIdx(c_idx)
    methyl.SetNumExplicitHs(3)

    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.METHYLATE,
            target=target_spec,
            message="Sanitization failed after methylation",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.METHYLATE,
        target=target_spec,
        message=f"Added methyl group to atom {target_idx}",
        affected_atoms=[target_idx],
        new_atoms=[c_idx],
    )


def execute_halogenate(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Add halogen to target atom.

    Params:
        element: "F", "Cl", "Br", or "I" (default: "Cl")
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    element = params.get("element", "Cl").upper()
    element_map = {"F": 9, "CL": 17, "BR": 35, "I": 53}

    if element not in element_map:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, target_atoms[0]),
            message=f"Invalid halogen: {element}",
            error_code="INVALID_ELEMENT",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    atom = rwmol.GetAtomWithIdx(target_idx)
    if atom is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=target_spec,
            message=f"Atom index {target_idx} not found",
            error_code="ATOM_NOT_FOUND",
        )

    num_hs = atom.GetTotalNumHs()
    if num_hs < 1:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=target_spec,
            message=f"Atom {target_idx} has no hydrogen to replace",
            error_code="NO_HYDROGEN",
        )

    # Add halogen
    hal_idx = rwmol.AddAtom(Chem.Atom(element_map[element]))
    rwmol.AddBond(target_idx, hal_idx, Chem.BondType.SINGLE)

    atom.SetNumExplicitHs(num_hs - 1)

    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.HALOGENATE,
            target=target_spec,
            message="Sanitization failed after halogenation",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.HALOGENATE,
        target=target_spec,
        message=f"Added {element} to atom {target_idx}",
        affected_atoms=[target_idx],
        new_atoms=[hal_idx],
    )


def execute_acetylate(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Add acetyl group (-COCH3) to target atom.

    The acetyl group replaces a hydrogen on a heteroatom (O, N, S).
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    atom = rwmol.GetAtomWithIdx(target_idx)
    if atom is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=target_spec,
            message=f"Atom index {target_idx} not found",
            error_code="ATOM_NOT_FOUND",
        )

    # Acetylation typically on O, N, or S
    symbol = atom.GetSymbol()
    if symbol not in ["O", "N", "S"]:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=target_spec,
            message=f"Acetylation requires O, N, or S atom, got {symbol}",
            error_code="INVALID_ATOM_TYPE",
            warnings=["Consider using a different operation for carbon atoms"],
        )

    num_hs = atom.GetTotalNumHs()
    if num_hs < 1:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=target_spec,
            message=f"Atom {target_idx} has no hydrogen to replace",
            error_code="NO_HYDROGEN",
        )

    # Add acetyl group: C(=O)CH3
    # First, add carbonyl carbon
    co_idx = rwmol.AddAtom(Chem.Atom(6))
    # Add carbonyl oxygen
    o_idx = rwmol.AddAtom(Chem.Atom(8))
    # Add methyl carbon
    ch3_idx = rwmol.AddAtom(Chem.Atom(6))

    # Bond carbonyl carbon to target
    rwmol.AddBond(target_idx, co_idx, Chem.BondType.SINGLE)
    # Double bond to oxygen
    rwmol.AddBond(co_idx, o_idx, Chem.BondType.DOUBLE)
    # Single bond to methyl
    rwmol.AddBond(co_idx, ch3_idx, Chem.BondType.SINGLE)

    # Update hydrogens
    atom.SetNumExplicitHs(num_hs - 1)
    methyl = rwmol.GetAtomWithIdx(ch3_idx)
    methyl.SetNumExplicitHs(3)

    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ACETYLATE,
            target=target_spec,
            message="Sanitization failed after acetylation",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.ACETYLATE,
        target=target_spec,
        message=f"Added acetyl group to atom {target_idx}",
        affected_atoms=[target_idx],
        new_atoms=[co_idx, o_idx, ch3_idx],
    )


def execute_epimerize(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Invert stereochemistry at a chiral center.

    Params:
        from_config: "R" or "S" (optional, for validation)
        to_config: "R" or "S" (optional, for validation)
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.EPIMERIZE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.EPIMERIZE,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    atom = rwmol.GetAtomWithIdx(target_idx)
    if atom is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.EPIMERIZE,
            target=target_spec,
            message=f"Atom index {target_idx} not found",
            error_code="ATOM_NOT_FOUND",
        )

    # Get current chiral tag
    current_tag = atom.GetChiralTag()
    from rdkit.Chem.rdchem import ChiralType

    # Map chiral tags
    if current_tag == ChiralType.CHI_TETRAHEDRAL_CW:
        new_tag = ChiralType.CHI_TETRAHEDRAL_CCW
    elif current_tag == ChiralType.CHI_TETRAHEDRAL_CCW:
        new_tag = ChiralType.CHI_TETRAHEDRAL_CW
    elif current_tag == ChiralType.CHI_UNSPECIFIED:
        return EditResult(
            success=False,
            operation_type=EditOperationType.EPIMERIZE,
            target=target_spec,
            message=f"Atom {target_idx} is not a chiral center",
            error_code="NOT_CHIRAL",
        )
    else:
        return EditResult(
            success=False,
            operation_type=EditOperationType.EPIMERIZE,
            target=target_spec,
            message=f"Unknown chirality type at atom {target_idx}",
            error_code="UNKNOWN_CHIRALITY",
        )

    # Invert
    atom.SetChiralTag(new_tag)

    return EditResult(
        success=True,
        operation_type=EditOperationType.EPIMERIZE,
        target=target_spec,
        message=f"Inverted stereochemistry at atom {target_idx}",
        affected_atoms=[target_idx],
    )


def execute_delete_atom(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Delete an atom and its bonds.

    The atom is removed and neighboring atoms get their hydrogen counts adjusted.
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_ATOM,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if not target_atoms:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_ATOM,
            target=TargetSpecification(TargetType.ATOM_INDEX, 0),
            message="No target atoms specified",
            error_code="NO_TARGET",
        )

    target_idx = target_atoms[0]
    target_spec = TargetSpecification(TargetType.ATOM_INDEX, target_idx)

    if target_idx >= rwmol.GetNumAtoms():
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_ATOM,
            target=target_spec,
            message=f"Atom index {target_idx} out of range",
            error_code="ATOM_NOT_FOUND",
        )

    # Get neighbors before removal
    atom = rwmol.GetAtomWithIdx(target_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]

    # Remove atom
    rwmol.RemoveAtom(target_idx)

    # Sanitize
    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_ATOM,
            target=target_spec,
            message="Sanitization failed after atom removal",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.DELETE_ATOM,
        target=target_spec,
        message=f"Removed atom {target_idx}",
        affected_atoms=neighbors,
        removed_atoms=[target_idx],
    )


def execute_add_bond(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Add a bond between two atoms.

    Params:
        order: 1, 2, or 3 (default: 1)
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ADD_BOND,
            target=TargetSpecification(TargetType.ATOM_PAIR, (0, 0)),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if len(target_atoms) < 2:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ADD_BOND,
            target=TargetSpecification(TargetType.ATOM_PAIR, (0, 0)),
            message="Need two atoms to add a bond",
            error_code="INSUFFICIENT_ATOMS",
        )

    atom1_idx, atom2_idx = target_atoms[0], target_atoms[1]
    target_spec = TargetSpecification(TargetType.ATOM_PAIR, (atom1_idx, atom2_idx))

    # Check atoms exist
    if atom1_idx >= rwmol.GetNumAtoms() or atom2_idx >= rwmol.GetNumAtoms():
        return EditResult(
            success=False,
            operation_type=EditOperationType.ADD_BOND,
            target=target_spec,
            message="Atom index out of range",
            error_code="ATOM_NOT_FOUND",
        )

    # Check if bond already exists
    existing_bond = rwmol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
    if existing_bond is not None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ADD_BOND,
            target=target_spec,
            message=f"Bond already exists between atoms {atom1_idx} and {atom2_idx}",
            error_code="BOND_EXISTS",
        )

    # Determine bond order
    order = params.get("order", 1)
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    bond_type = bond_type_map.get(order, Chem.BondType.SINGLE)

    # Add bond
    rwmol.AddBond(atom1_idx, atom2_idx, bond_type)

    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.ADD_BOND,
            target=target_spec,
            message="Sanitization failed after adding bond",
            error_code="SANITIZE_FAILED",
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.ADD_BOND,
        target=target_spec,
        message=f"Added bond (order {order}) between atoms {atom1_idx} and {atom2_idx}",
        affected_atoms=[atom1_idx, atom2_idx],
    )


def execute_delete_bond(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """Delete a bond between two atoms."""
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_BOND,
            target=TargetSpecification(TargetType.ATOM_PAIR, (0, 0)),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    if len(target_atoms) < 2:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_BOND,
            target=TargetSpecification(TargetType.ATOM_PAIR, (0, 0)),
            message="Need two atoms to delete a bond",
            error_code="INSUFFICIENT_ATOMS",
        )

    atom1_idx, atom2_idx = target_atoms[0], target_atoms[1]
    target_spec = TargetSpecification(TargetType.ATOM_PAIR, (atom1_idx, atom2_idx))

    # Check bond exists
    bond = rwmol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
    if bond is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_BOND,
            target=target_spec,
            message=f"No bond between atoms {atom1_idx} and {atom2_idx}",
            error_code="BOND_NOT_FOUND",
        )

    # Remove bond
    rwmol.RemoveBond(atom1_idx, atom2_idx)

    sanitize_result = Chem.SanitizeMol(rwmol, catchErrors=True)
    if sanitize_result != Chem.SanitizeFlags.SANITIZE_NONE:
        return EditResult(
            success=False,
            operation_type=EditOperationType.DELETE_BOND,
            target=target_spec,
            message="Sanitization failed after removing bond",
            error_code="SANITIZE_FAILED",
            warnings=["Molecule may now have multiple fragments"],
        )

    return EditResult(
        success=True,
        operation_type=EditOperationType.DELETE_BOND,
        target=target_spec,
        message=f"Removed bond between atoms {atom1_idx} and {atom2_idx}",
        affected_atoms=[atom1_idx, atom2_idx],
    )


def execute_substitute(rwmol, target_atoms: List[int], params: Dict[str, Any]) -> EditResult:
    """
    Substitute a substructure with another.

    Params:
        replacement: SMARTS pattern for replacement group
        preserve_stereo: Whether to preserve stereochemistry (default: True)
    """
    if not RDKIT_AVAILABLE or rwmol is None:
        return EditResult(
            success=False,
            operation_type=EditOperationType.SUBSTITUTE,
            target=TargetSpecification(TargetType.SMARTS, ""),
            message="RDKit not available",
            error_code="RDKIT_NOT_AVAILABLE",
        )

    replacement_smarts = params.get("replacement", "")
    if not replacement_smarts:
        return EditResult(
            success=False,
            operation_type=EditOperationType.SUBSTITUTE,
            target=TargetSpecification(TargetType.SMARTS, ""),
            message="No replacement pattern specified",
            error_code="NO_REPLACEMENT",
        )

    # This is a complex operation - for now return a placeholder
    # Full implementation would use ReplaceSubstructs
    return EditResult(
        success=False,
        operation_type=EditOperationType.SUBSTITUTE,
        target=TargetSpecification(TargetType.SMARTS, ""),
        message="Substitution operation not yet fully implemented",
        error_code="NOT_IMPLEMENTED",
    )


# Operation executor registry
OPERATION_EXECUTORS: Dict[EditOperationType, Callable] = {
    EditOperationType.HYDROXYLATE: execute_hydroxylate,
    EditOperationType.METHYLATE: execute_methylate,
    EditOperationType.HALOGENATE: execute_halogenate,
    EditOperationType.ACETYLATE: execute_acetylate,
    EditOperationType.EPIMERIZE: execute_epimerize,
    EditOperationType.DELETE_ATOM: execute_delete_atom,
    EditOperationType.ADD_BOND: execute_add_bond,
    EditOperationType.DELETE_BOND: execute_delete_bond,
    EditOperationType.SUBSTITUTE: execute_substitute,
}


def execute_operation(
    rwmol,
    operation: EditOperation,
    annotated_mol=None
) -> EditResult:
    """
    Execute a single edit operation on an RWMol.

    Args:
        rwmol: RDKit RWMol object
        operation: EditOperation to execute
        annotated_mol: Optional AnnotatedMolecule for target resolution

    Returns:
        EditResult with success/failure information
    """
    # Resolve target atoms
    target_atoms = resolve_target(rwmol, operation.target, annotated_mol)

    if not target_atoms and operation.operation_type != EditOperationType.CUSTOM:
        return EditResult(
            success=False,
            operation_type=operation.operation_type,
            target=operation.target,
            message=f"Could not resolve target: {operation.target.to_dict()}",
            error_code="TARGET_NOT_FOUND",
        )

    # Get executor for this operation type
    executor = OPERATION_EXECUTORS.get(operation.operation_type)

    if executor is None:
        return EditResult(
            success=False,
            operation_type=operation.operation_type,
            target=operation.target,
            message=f"Operation {operation.operation_type.value} not implemented",
            error_code="NOT_IMPLEMENTED",
        )

    # Execute
    return executor(rwmol, target_atoms, operation.params)
