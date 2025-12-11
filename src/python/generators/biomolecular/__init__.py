"""
biomolecular/ - Biomolecular and Soft Matter Structures

Comprehensive module for biomolecular structures.
Covers Category 12 of structure_catalogue.md (6 items).
"""

from .lipids import (
    generate_lipid_bilayer,
    generate_membrane_with_protein,
    LIPID_DATABASE,
)

from .interfaces import (
    generate_protein_surface,
    generate_dna_surface,
    generate_peptide_surface,
    SURFACE_PROTEINS,
)


__all__ = [
    "generate_lipid_bilayer", "generate_membrane_with_protein", "LIPID_DATABASE",
    "generate_protein_surface", "generate_dna_surface", "generate_peptide_surface",
    "SURFACE_PROTEINS",
]
