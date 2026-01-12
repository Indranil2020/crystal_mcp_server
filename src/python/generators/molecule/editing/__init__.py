"""
Molecular Editing Module

Provides tools for modifying molecular structures through natural language
or structured operations. Supports stereochemistry, ring systems, functional
groups, and atom-level modifications.

Key Components:
- AnnotatedMolecule: Rich molecular representation with annotations
- MolecularEditor: RWMol-based editing engine
- EditOperation: Structured operation definitions
- SemanticParser: NL to operations conversion
- ValidationPipeline: Post-edit validation
"""

from .annotated_molecule import (
    AnnotatedMolecule,
    StereocenterInfo,
    RingInfo,
    FunctionalGroupInfo,
    AtomInfo,
)

from .edit_operations import (
    EditOperation,
    EditOperationType,
    TargetSpecification,
    EditResult,
    create_operation,
)

from .molecular_editor import (
    MolecularEditor,
)

from .semantic_parser import (
    parse_nl_instruction,
    SemanticParser,
)

from .validation import (
    ValidationResult,
    validate_molecule,
    ValidationPipeline,
)

__all__ = [
    # Data structures
    "AnnotatedMolecule",
    "StereocenterInfo",
    "RingInfo",
    "FunctionalGroupInfo",
    "AtomInfo",
    # Operations
    "EditOperation",
    "EditOperationType",
    "TargetSpecification",
    "EditResult",
    "create_operation",
    # Editor
    "MolecularEditor",
    # Parser
    "parse_nl_instruction",
    "SemanticParser",
    # Validation
    "ValidationResult",
    "validate_molecule",
    "ValidationPipeline",
]
