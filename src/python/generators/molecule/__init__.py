"""
molecule/ - Molecular Structures

Covers Category 5 of structure_catalogue.md.
"""

from .small_molecules import (
    generate_molecule,
    get_molecule_database,
    MOLECULE_DATABASE,
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
    # Small molecules
    "generate_molecule", "get_molecule_database", "MOLECULE_DATABASE",
    # Frameworks
    "generate_mof", "generate_cof",
    # Cages
    "generate_fullerene", "generate_cage_compound",
]
