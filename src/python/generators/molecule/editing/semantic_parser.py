"""
Semantic Parser Module

Converts natural language molecular editing instructions into structured
EditOperation objects. Uses pattern matching with a chemical knowledge base.

Supports instructions like:
- "Hydroxylate at C-7 with axial orientation"
- "Methylate the nitrogen"
- "Epimerize C-3 from R to S"
- "Add a chlorine to the benzene ring"
- "Remove the methyl group at C-4"
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any, Pattern
import re
import logging

from .edit_operations import (
    EditOperation,
    EditOperationType,
    TargetSpecification,
    TargetType,
    create_operation,
)
from .annotated_molecule import AnnotatedMolecule, FunctionalGroupType, RingType

logger = logging.getLogger(__name__)


# ============================================================================
# Chemical Knowledge Base
# ============================================================================

# Verb patterns mapping to operation types
VERB_PATTERNS: Dict[str, EditOperationType] = {
    # Additions
    "hydroxylate": EditOperationType.HYDROXYLATE,
    "add hydroxyl": EditOperationType.HYDROXYLATE,
    "add oh": EditOperationType.HYDROXYLATE,
    "add -oh": EditOperationType.HYDROXYLATE,

    "methylate": EditOperationType.METHYLATE,
    "add methyl": EditOperationType.METHYLATE,
    "add ch3": EditOperationType.METHYLATE,
    "add -ch3": EditOperationType.METHYLATE,

    "acetylate": EditOperationType.ACETYLATE,
    "add acetyl": EditOperationType.ACETYLATE,
    "add ac": EditOperationType.ACETYLATE,

    "halogenate": EditOperationType.HALOGENATE,
    "add halogen": EditOperationType.HALOGENATE,
    "chlorinate": EditOperationType.HALOGENATE,
    "add chlorine": EditOperationType.HALOGENATE,
    "add cl": EditOperationType.HALOGENATE,
    "brominate": EditOperationType.HALOGENATE,
    "add bromine": EditOperationType.HALOGENATE,
    "add br": EditOperationType.HALOGENATE,
    "fluorinate": EditOperationType.HALOGENATE,
    "add fluorine": EditOperationType.HALOGENATE,
    "add f": EditOperationType.HALOGENATE,
    "iodinate": EditOperationType.HALOGENATE,
    "add iodine": EditOperationType.HALOGENATE,
    "add i": EditOperationType.HALOGENATE,

    "aminate": EditOperationType.AMINATE,
    "add amino": EditOperationType.AMINATE,
    "add amine": EditOperationType.AMINATE,
    "add nh2": EditOperationType.AMINATE,
    "add -nh2": EditOperationType.AMINATE,

    "carboxylate": EditOperationType.CARBOXYLATE,
    "add carboxyl": EditOperationType.CARBOXYLATE,
    "add cooh": EditOperationType.CARBOXYLATE,
    "add -cooh": EditOperationType.CARBOXYLATE,

    "phosphorylate": EditOperationType.PHOSPHORYLATE,
    "add phosphate": EditOperationType.PHOSPHORYLATE,
    "add po4": EditOperationType.PHOSPHORYLATE,

    "sulfonate": EditOperationType.SULFONATE,
    "add sulfonate": EditOperationType.SULFONATE,
    "add so3h": EditOperationType.SULFONATE,

    "nitrate": EditOperationType.NITRATE,
    "add nitro": EditOperationType.NITRATE,
    "add no2": EditOperationType.NITRATE,

    # Removals
    "dehydroxylate": EditOperationType.DEHYDROXYLATE,
    "remove hydroxyl": EditOperationType.DEHYDROXYLATE,
    "remove oh": EditOperationType.DEHYDROXYLATE,
    "remove -oh": EditOperationType.DEHYDROXYLATE,

    "demethylate": EditOperationType.DEMETHYLATE,
    "remove methyl": EditOperationType.DEMETHYLATE,
    "remove ch3": EditOperationType.DEMETHYLATE,
    "remove -ch3": EditOperationType.DEMETHYLATE,

    "deacetylate": EditOperationType.DEACETYLATE,
    "remove acetyl": EditOperationType.DEACETYLATE,
    "remove ac": EditOperationType.DEACETYLATE,

    "dehalogenate": EditOperationType.DEHALOGENATE,
    "remove halogen": EditOperationType.DEHALOGENATE,
    "remove chlorine": EditOperationType.DEHALOGENATE,
    "remove cl": EditOperationType.DEHALOGENATE,
    "remove bromine": EditOperationType.DEHALOGENATE,
    "remove br": EditOperationType.DEHALOGENATE,

    # Stereochemistry
    "epimerize": EditOperationType.EPIMERIZE,
    "invert": EditOperationType.INVERT_STEREO,
    "invert stereochemistry": EditOperationType.INVERT_STEREO,
    "invert stereo": EditOperationType.INVERT_STEREO,
    "flip stereochemistry": EditOperationType.INVERT_STEREO,
    "flip stereo": EditOperationType.INVERT_STEREO,
    "set stereo": EditOperationType.SET_STEREO,
    "set stereochemistry": EditOperationType.SET_STEREO,
    "set e/z": EditOperationType.SET_E_Z,
    "set ez": EditOperationType.SET_E_Z,
    "set cis": EditOperationType.SET_E_Z,
    "set trans": EditOperationType.SET_E_Z,

    # Atom operations
    "add atom": EditOperationType.ADD_ATOM,
    "delete atom": EditOperationType.DELETE_ATOM,
    "remove atom": EditOperationType.DELETE_ATOM,
    "change element": EditOperationType.CHANGE_ELEMENT,
    "replace element": EditOperationType.CHANGE_ELEMENT,
    "set charge": EditOperationType.SET_CHARGE,
    "change charge": EditOperationType.SET_CHARGE,

    # Bond operations
    "add bond": EditOperationType.ADD_BOND,
    "create bond": EditOperationType.ADD_BOND,
    "form bond": EditOperationType.ADD_BOND,
    "delete bond": EditOperationType.DELETE_BOND,
    "remove bond": EditOperationType.DELETE_BOND,
    "break bond": EditOperationType.DELETE_BOND,
    "change bond order": EditOperationType.CHANGE_BOND_ORDER,
    "make double bond": EditOperationType.CHANGE_BOND_ORDER,
    "make triple bond": EditOperationType.CHANGE_BOND_ORDER,
    "make single bond": EditOperationType.CHANGE_BOND_ORDER,

    # Ring operations
    "add ring": EditOperationType.ADD_RING,
    "delete ring": EditOperationType.DELETE_RING,
    "remove ring": EditOperationType.DELETE_RING,
    "fuse ring": EditOperationType.FUSE_RING,
    "open ring": EditOperationType.OPEN_RING,
    "break ring": EditOperationType.OPEN_RING,
    "expand ring": EditOperationType.EXPAND_RING,
    "contract ring": EditOperationType.CONTRACT_RING,

    # Substitution
    "substitute": EditOperationType.SUBSTITUTE,
    "replace": EditOperationType.REPLACE_GROUP,
    "replace group": EditOperationType.REPLACE_GROUP,

    # Protection
    "protect": EditOperationType.PROTECT,
    "deprotect": EditOperationType.DEPROTECT,
    "remove protecting group": EditOperationType.DEPROTECT,

    # Oxidation/Reduction
    "oxidize": EditOperationType.OXIDIZE,
    "reduce": EditOperationType.REDUCE,
}

# Halogen element patterns
HALOGEN_PATTERNS = {
    "chlorine": "Cl",
    "chloro": "Cl",
    "cl": "Cl",
    "bromine": "Br",
    "bromo": "Br",
    "br": "Br",
    "fluorine": "F",
    "fluoro": "F",
    "f": "F",
    "iodine": "I",
    "iodo": "I",
    "i": "I",
}

# Orientation patterns
ORIENTATION_PATTERNS = {
    "axial": "axial",
    "equatorial": "equatorial",
    "ax": "axial",
    "eq": "equatorial",
    "alpha": "alpha",
    "beta": "beta",
    "α": "alpha",
    "β": "beta",
}

# Stereo configuration patterns
STEREO_PATTERNS = {
    "r": "R",
    "s": "S",
    "(r)": "R",
    "(s)": "S",
    "r-config": "R",
    "s-config": "S",
    "r configuration": "R",
    "s configuration": "S",
    "e": "E",
    "z": "Z",
    "(e)": "E",
    "(z)": "Z",
    "cis": "Z",
    "trans": "E",
}

# Ring type patterns
RING_NAME_PATTERNS = {
    "benzene": RingType.BENZENE,
    "phenyl": RingType.BENZENE,
    "aromatic ring": RingType.BENZENE,
    "naphthalene": RingType.NAPHTHALENE,
    "naphthyl": RingType.NAPHTHALENE,
    "cyclohexane": RingType.CYCLOHEXANE,
    "cyclopentane": RingType.CYCLOPENTANE,
    "furan": RingType.FURAN,
    "furyl": RingType.FURAN,
    "thiophene": RingType.THIOPHENE,
    "thienyl": RingType.THIOPHENE,
    "pyrrole": RingType.PYRROLE,
    "pyrrolyl": RingType.PYRROLE,
    "pyridine": RingType.PYRIDINE,
    "pyridyl": RingType.PYRIDINE,
    "imidazole": RingType.IMIDAZOLE,
    "imidazolyl": RingType.IMIDAZOLE,
    "oxetane": RingType.OXETANE,
    "epoxide": RingType.EPOXIDE,
    "3-membered ring": RingType.EPOXIDE,
    "4-membered ring": RingType.OXETANE,
}

# Functional group patterns
GROUP_NAME_PATTERNS = {
    "hydroxyl": FunctionalGroupType.HYDROXYL,
    "oh": FunctionalGroupType.HYDROXYL,
    "alcohol": FunctionalGroupType.HYDROXYL,
    "carboxyl": FunctionalGroupType.CARBOXYL,
    "cooh": FunctionalGroupType.CARBOXYL,
    "carboxylic acid": FunctionalGroupType.CARBOXYL,
    "amino": FunctionalGroupType.AMINO,
    "amine": FunctionalGroupType.AMINO,
    "nitrogen": FunctionalGroupType.AMINO,  # "the nitrogen" usually refers to amine
    "nh2": FunctionalGroupType.AMINO,
    "carbonyl": FunctionalGroupType.CARBONYL,
    "ketone": FunctionalGroupType.CARBONYL,
    "aldehyde": FunctionalGroupType.CARBONYL,
    "ester": FunctionalGroupType.ESTER,
    "ether": FunctionalGroupType.ETHER,
    "amide": FunctionalGroupType.AMIDE,
    "nitro": FunctionalGroupType.NITRO,
    "nitrile": FunctionalGroupType.NITRILE,
    "cyano": FunctionalGroupType.NITRILE,
    "cn": FunctionalGroupType.NITRILE,
    "sulfhydryl": FunctionalGroupType.SULFHYDRYL,
    "thiol": FunctionalGroupType.SULFHYDRYL,
    "sh": FunctionalGroupType.SULFHYDRYL,
    "halide": FunctionalGroupType.HALIDE,
    "halogen": FunctionalGroupType.HALIDE,
    "methyl": FunctionalGroupType.METHYL,
    "ch3": FunctionalGroupType.METHYL,
    "ethyl": FunctionalGroupType.ETHYL,
    "phenyl": FunctionalGroupType.PHENYL,
    "acetyl": FunctionalGroupType.ACETYL,
    "ac": FunctionalGroupType.ACETYL,
}


@dataclass
class ParseResult:
    """Result of parsing a natural language instruction."""
    success: bool
    operations: List[EditOperation]
    message: str = ""
    confidence: float = 0.0
    warnings: List[str] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


class SemanticParser:
    """
    Parser for natural language molecular editing instructions.

    Uses pattern matching with a chemical knowledge base to extract:
    - Operation type (what to do)
    - Target specification (where to apply)
    - Parameters (how to do it)
    """

    def __init__(self, molecule: Optional[AnnotatedMolecule] = None):
        """
        Initialize parser with optional molecule context.

        Args:
            molecule: AnnotatedMolecule for context-aware parsing
        """
        self.molecule = molecule

    def parse(self, instruction: str) -> ParseResult:
        """
        Parse a natural language instruction into EditOperation(s).

        Args:
            instruction: Natural language instruction

        Returns:
            ParseResult with operations and metadata
        """
        instruction = instruction.strip().lower()

        if not instruction:
            return ParseResult(
                success=False,
                operations=[],
                message="Empty instruction",
            )

        # Split on common conjunction patterns to handle multiple operations
        sub_instructions = self._split_instructions(instruction)

        operations = []
        warnings = []
        total_confidence = 0.0

        for sub_instr in sub_instructions:
            result = self._parse_single_instruction(sub_instr)
            if result.success:
                operations.extend(result.operations)
                total_confidence += result.confidence
            else:
                warnings.append(f"Could not parse: '{sub_instr}'")

        if not operations:
            return ParseResult(
                success=False,
                operations=[],
                message="Could not parse any operations from instruction",
                warnings=warnings,
            )

        avg_confidence = total_confidence / len(sub_instructions) if sub_instructions else 0.0

        return ParseResult(
            success=True,
            operations=operations,
            message=f"Parsed {len(operations)} operation(s)",
            confidence=avg_confidence,
            warnings=warnings,
        )

    def _split_instructions(self, instruction: str) -> List[str]:
        """Split instruction into sub-instructions on conjunctions."""
        # Split on common patterns
        patterns = [
            r'\s+and\s+',
            r'\s*,\s+then\s+',
            r'\s*;\s*',
            r'\s+then\s+',
            r'\s+also\s+',
        ]

        result = [instruction]
        for pattern in patterns:
            new_result = []
            for item in result:
                new_result.extend(re.split(pattern, item))
            result = new_result

        return [s.strip() for s in result if s.strip()]

    def _parse_single_instruction(self, instruction: str) -> ParseResult:
        """Parse a single instruction (no conjunctions)."""
        # Step 1: Extract operation type
        op_type, op_confidence = self._extract_operation_type(instruction)

        if op_type is None:
            return ParseResult(
                success=False,
                operations=[],
                message=f"Unknown operation in: {instruction}",
            )

        # Step 2: Extract target
        target, target_confidence = self._extract_target(instruction)

        if target is None:
            # Use a generic target for operations that don't require one
            target = TargetSpecification(TargetType.ALL, "all")
            target_confidence = 0.5

        # Step 3: Extract parameters
        params = self._extract_parameters(instruction, op_type)

        # Create operation
        operation = EditOperation(
            operation_type=op_type,
            target=target,
            params=params,
        )

        confidence = (op_confidence + target_confidence) / 2

        return ParseResult(
            success=True,
            operations=[operation],
            confidence=confidence,
        )

    def _extract_operation_type(self, instruction: str) -> Tuple[Optional[EditOperationType], float]:
        """
        Extract the operation type from instruction.

        Returns:
            Tuple of (operation type, confidence)
        """
        instruction_lower = instruction.lower()

        # Try exact matches first (highest confidence)
        for pattern, op_type in VERB_PATTERNS.items():
            if pattern in instruction_lower:
                return (op_type, 0.95)

        # Try fuzzy matching on verbs
        words = instruction_lower.split()
        for word in words:
            for pattern, op_type in VERB_PATTERNS.items():
                # Check if word starts with pattern verb
                pattern_verb = pattern.split()[0]
                if word.startswith(pattern_verb[:4]) and len(pattern_verb) > 3:
                    return (op_type, 0.7)

        return (None, 0.0)

    def _extract_target(self, instruction: str) -> Tuple[Optional[TargetSpecification], float]:
        """
        Extract target specification from instruction.

        Returns:
            Tuple of (target specification, confidence)
        """
        instruction_lower = instruction.lower()

        # Pattern 1: Atom label like "C-7", "N-3", "O-1"
        label_match = re.search(r'\b([A-Z][a-z]?)-(\d+)\b', instruction, re.IGNORECASE)
        if label_match:
            label = f"{label_match.group(1).upper()}-{label_match.group(2)}"
            return (TargetSpecification(TargetType.ATOM_LABEL, label), 0.95)

        # Pattern 2: Position references like "at position 7", "position 3"
        position_match = re.search(r'(?:at\s+)?position\s+(\d+)', instruction_lower)
        if position_match:
            pos = int(position_match.group(1))
            # Check if molecule context available to resolve
            if self.molecule is not None:
                # Try to find atom at this position
                for idx, label in self.molecule.atom_labels.items():
                    if label.endswith(f"-{pos}"):
                        return (TargetSpecification(TargetType.ATOM_LABEL, label), 0.85)
            return (TargetSpecification(TargetType.ATOM_INDEX, pos - 1), 0.7)

        # Pattern 3: Ring type like "benzene ring", "oxetane"
        for ring_name, ring_type in RING_NAME_PATTERNS.items():
            if ring_name in instruction_lower:
                # Check for index like "first oxetane", "second benzene"
                ordinal_match = re.search(rf'(first|second|third|1st|2nd|3rd|\d+)(?:st|nd|rd|th)?\s+{ring_name}', instruction_lower)
                idx = 0
                if ordinal_match:
                    ordinal = ordinal_match.group(1)
                    ordinal_map = {"first": 0, "1st": 0, "second": 1, "2nd": 1, "third": 2, "3rd": 2}
                    if ordinal in ordinal_map:
                        idx = ordinal_map[ordinal]
                    elif ordinal.isdigit():
                        idx = int(ordinal) - 1
                return (TargetSpecification(TargetType.RING_TYPE, ring_type.value, idx), 0.85)

        # Pattern 4: Functional group like "the hydroxyl", "primary amine"
        for group_name, group_type in GROUP_NAME_PATTERNS.items():
            if group_name in instruction_lower:
                ordinal_match = re.search(rf'(first|second|third|primary|secondary|tertiary|\d+)(?:st|nd|rd|th)?\s+{group_name}', instruction_lower)
                idx = 0
                if ordinal_match:
                    ordinal = ordinal_match.group(1)
                    ordinal_map = {
                        "first": 0, "primary": 0, "1st": 0,
                        "second": 1, "secondary": 1, "2nd": 1,
                        "third": 2, "tertiary": 2, "3rd": 2
                    }
                    if ordinal in ordinal_map:
                        idx = ordinal_map[ordinal]
                    elif ordinal.isdigit():
                        idx = int(ordinal) - 1
                return (TargetSpecification(TargetType.FUNCTIONAL_GROUP, group_type.value, idx), 0.8)

        # Pattern 5: Element reference like "the nitrogen", "carbon atom"
        element_patterns = [
            (r'\bthe\s+(carbon|nitrogen|oxygen|sulfur|phosphorus)\b', {"carbon": "C", "nitrogen": "N", "oxygen": "O", "sulfur": "S", "phosphorus": "P"}),
            (r'\b(carbon|nitrogen|oxygen|sulfur|phosphorus)\s+atom\b', {"carbon": "C", "nitrogen": "N", "oxygen": "O", "sulfur": "S", "phosphorus": "P"}),
        ]
        for pattern, elem_map in element_patterns:
            match = re.search(pattern, instruction_lower)
            if match:
                elem = elem_map.get(match.group(1), "C")
                # Use SMARTS to find this element
                return (TargetSpecification(TargetType.SMARTS, f"[{elem}]"), 0.6)

        # Pattern 6: SMARTS-like pattern in brackets
        smarts_match = re.search(r'\[([^\]]+)\]', instruction)
        if smarts_match:
            smarts = f"[{smarts_match.group(1)}]"
            return (TargetSpecification(TargetType.SMARTS, smarts), 0.9)

        # Pattern 7: Generic "the ring", "a ring"
        if re.search(r'\b(?:the|a)\s+ring\b', instruction_lower):
            return (TargetSpecification(TargetType.RING_INDEX, 0), 0.5)

        return (None, 0.0)

    def _extract_parameters(self, instruction: str, op_type: EditOperationType) -> Dict[str, Any]:
        """
        Extract operation-specific parameters from instruction.

        Returns:
            Dictionary of parameters
        """
        params = {}
        instruction_lower = instruction.lower()

        # Orientation parameters
        for pattern, orientation in ORIENTATION_PATTERNS.items():
            if pattern in instruction_lower:
                params["orientation"] = orientation
                break

        # Halogen element
        if op_type == EditOperationType.HALOGENATE:
            for pattern, element in HALOGEN_PATTERNS.items():
                if pattern in instruction_lower:
                    params["element"] = element
                    break
            # Default to Cl if no specific halogen mentioned
            if "element" not in params:
                params["element"] = "Cl"

        # Stereo configuration
        if op_type in (EditOperationType.EPIMERIZE, EditOperationType.SET_STEREO):
            # Look for "from R to S" or "R to S" patterns
            from_to_match = re.search(r'from\s+([rs])\s+to\s+([rs])', instruction_lower)
            if from_to_match:
                params["from_config"] = from_to_match.group(1).upper()
                params["to_config"] = from_to_match.group(2).upper()
            else:
                # Just look for target config
                for pattern, config in STEREO_PATTERNS.items():
                    if pattern in instruction_lower:
                        params["to_config"] = config
                        break

        # E/Z configuration
        if op_type == EditOperationType.SET_E_Z:
            if "cis" in instruction_lower or "(z)" in instruction_lower or " z " in instruction_lower:
                params["configuration"] = "Z"
            elif "trans" in instruction_lower or "(e)" in instruction_lower or " e " in instruction_lower:
                params["configuration"] = "E"

        # Bond order
        if op_type in (EditOperationType.ADD_BOND, EditOperationType.CHANGE_BOND_ORDER):
            if "double" in instruction_lower:
                params["order"] = 2
            elif "triple" in instruction_lower:
                params["order"] = 3
            elif "single" in instruction_lower:
                params["order"] = 1

        # Ring size
        if op_type in (EditOperationType.ADD_RING, EditOperationType.EXPAND_RING, EditOperationType.CONTRACT_RING):
            size_match = re.search(r'(\d+)[- ]membered', instruction_lower)
            if size_match:
                params["ring_size"] = int(size_match.group(1))

        # Substitution replacement
        if op_type in (EditOperationType.SUBSTITUTE, EditOperationType.REPLACE_GROUP):
            with_match = re.search(r'with\s+(\S+)', instruction_lower)
            if with_match:
                replacement = with_match.group(1)
                # Check if it's a known group name
                for name, group_type in GROUP_NAME_PATTERNS.items():
                    if name in replacement:
                        params["replacement"] = group_type.value
                        break
                if "replacement" not in params:
                    params["replacement"] = replacement

        return params


def parse_nl_instruction(
    instruction: str,
    molecule: Optional[AnnotatedMolecule] = None
) -> List[EditOperation]:
    """
    Parse natural language instruction into edit operations.

    This is the main entry point for NL parsing.

    Args:
        instruction: Natural language instruction
        molecule: Optional molecule for context

    Returns:
        List of EditOperation objects
    """
    parser = SemanticParser(molecule)
    result = parser.parse(instruction)

    if result.success:
        return result.operations

    logger.warning(f"Failed to parse instruction: {instruction}")
    logger.warning(f"Reason: {result.message}")

    return []


def get_supported_operations() -> List[str]:
    """Get list of supported operation verbs."""
    return list(VERB_PATTERNS.keys())


def get_supported_ring_types() -> List[str]:
    """Get list of supported ring type names."""
    return list(RING_NAME_PATTERNS.keys())


def get_supported_functional_groups() -> List[str]:
    """Get list of supported functional group names."""
    return list(GROUP_NAME_PATTERNS.keys())
