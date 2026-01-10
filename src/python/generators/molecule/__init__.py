"""
molecule/ - Molecular Structures

Covers Category 5 of structure_catalogue.md.
"""

from .small_molecules import (
    generate_molecule,
    get_molecule_database as get_small_molecule_database,  # Legacy name
    MOLECULE_DATABASE,
)

from .universal_molecule import (
    generate_molecule_universal,
    smiles_to_3d_structure,
    resolve_to_smiles,
    get_supported_molecules,
    MOLECULE_ALIASES,
)

from .molecule_database import (
    MoleculeDatabase,
    get_molecule_database,
    lookup_molecule,
)

from .molecular_cluster import (
    generate_molecular_cluster,
    create_dimer,
    create_hetero_dimer,
    create_stack,
    StackingType,
    STACKING_DEFAULTS,
)

# Unified arrangement adapter (delegates to molecular_arrangement.py)
from .arrangement_adapter import (
    generate_molecular_cluster as generate_molecular_cluster_unified,
    create_arrangement_from_formula,
    create_arrangement_with_constraints,
)

from .frameworks import (
    generate_mof,
    generate_cof,
)

from .cages import (
    generate_fullerene,
    generate_cage_compound,
)


__all__ = [
    # Small molecules (local DB)
    "generate_molecule", "get_small_molecule_database", "MOLECULE_DATABASE",
    # Universal molecule generation (any identifier)
    "generate_molecule_universal", "smiles_to_3d_structure", "resolve_to_smiles",
    "get_supported_molecules", "MOLECULE_ALIASES",
    # Unified molecule database (SQLite)
    "MoleculeDatabase", "get_molecule_database", "lookup_molecule",
    # Molecular clusters
    "generate_molecular_cluster", "create_dimer", "create_hetero_dimer", "create_stack",
    "StackingType", "STACKING_DEFAULTS",
    # Unified arrangement (formula-based, constraint-based)
    "generate_molecular_cluster_unified", "create_arrangement_from_formula",
    "create_arrangement_with_constraints",
    # Frameworks
    "generate_mof", "generate_cof",
    # Cages
    "generate_fullerene", "generate_cage_compound",
]
