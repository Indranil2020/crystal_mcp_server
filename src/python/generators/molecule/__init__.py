"""
molecule/ - Molecular Structures

Comprehensive module for molecular structures.
Covers Category 5 of structure_catalogue.md (15 items).

Submodules:
- small_molecules: Common molecules
- frameworks: MOFs, COFs
- cages: Fullerenes, cages
- organometallics: Metallocenes, carbenes
- porphyrins: Porphyrins, phthalocyanines
- carbon_nanostructures: Nanohorns, peapods
- biomolecules: DNA, peptides
"""

# Small molecules
from .small_molecules import generate_molecule, COMMON_MOLECULES

# Frameworks
from .frameworks import generate_mof, generate_cof, MOF_DATABASE

# Cages and fullerenes
from .cages import (
    generate_fullerene,
    generate_cage,
    FULLERENE_DATABASE,
)

# Need CAGE_DATABASE - add it
CAGE_DATABASE = {
    "M12L24": {"n_metal": 12, "n_ligand": 24},
    "M24L48": {"n_metal": 24, "n_ligand": 48},
    "Fe4L6": {"n_metal": 4, "n_ligand": 6},
    "Pd6L4": {"n_metal": 6, "n_ligand": 4},
}

# Organometallics
from .organometallics import (
    generate_metallocene,
    generate_piano_stool,
    generate_nhc_carbene,
    METALLOCENE_DATABASE,
    PIANO_STOOL_DATABASE,
    CARBENE_CATALYSTS,
)

# Porphyrins
from .porphyrins import (
    generate_porphyrin,
    generate_phthalocyanine,
    generate_corrole,
    generate_subphthalocyanine,
    PORPHYRIN_DATABASE,
    PHTHALOCYANINE_DATABASE,
    CORROLE_DATABASE,
)

# Carbon nanostructures
from .carbon_nanostructures import (
    generate_nanohorn,
    generate_peapod,
    generate_nanobud,
    generate_graphene_quantum_dot,
    NANOHORN_DATABASE,
    PEAPOD_DATABASE,
)

# Biomolecules
from .biomolecules import (
    generate_peptide,
    generate_dna_strand,
    generate_g_quadruplex,
    generate_dipeptide_library,
    generate_tripeptide_library,
    AMINO_ACIDS,
    NUCLEOTIDES,
    HELIX_PARAMS,
)


__all__ = [
    # Small molecules
    "generate_molecule", "COMMON_MOLECULES",
    # Frameworks
    "generate_mof", "generate_cof", "MOF_DATABASE",
    # Cages
    "generate_fullerene", "generate_cage",
    "FULLERENE_DATABASE", "CAGE_DATABASE",
    # Organometallics
    "generate_metallocene", "generate_piano_stool", "generate_nhc_carbene",
    "METALLOCENE_DATABASE", "PIANO_STOOL_DATABASE", "CARBENE_CATALYSTS",
    # Porphyrins
    "generate_porphyrin", "generate_phthalocyanine", "generate_corrole",
    "generate_subphthalocyanine", "PORPHYRIN_DATABASE", "PHTHALOCYANINE_DATABASE",
    "CORROLE_DATABASE",
    # Carbon nanostructures
    "generate_nanohorn", "generate_peapod", "generate_nanobud",
    "generate_graphene_quantum_dot", "NANOHORN_DATABASE", "PEAPOD_DATABASE",
    # Biomolecules
    "generate_peptide", "generate_dna_strand", "generate_g_quadruplex",
    "generate_dipeptide_library", "generate_tripeptide_library",
    "AMINO_ACIDS", "NUCLEOTIDES", "HELIX_PARAMS",
]
