"""
twist/ - Twistronics and Multilayer Stacks

Comprehensive module for twisted multilayer structures.
Covers all 6 structure types from Category 3 of structure_catalogue.md.

Submodules:
- moire: Moiré superlattice generation
- bilayers: Twisted homo- and hetero-bilayers
- multilayers: Trilayers to pentalayers
- intercalation: Li, Na, K intercalated bilayers
- stacking: ABA, ABC, rhombohedral stacking
- ferroelectric: Layer-locked ferroelectrics
"""

# Base utilities
from .base import TwistStructure, calculate_moire_angle, find_commensurate_cell

# Moiré superlattices
from .moire import (
    generate_moire_superlattice,
    calculate_moire_periodicity,
    MAGIC_ANGLES,
)

# Twisted bilayers
from .bilayers import (
    generate_twisted_bilayer,
    generate_twisted_heterobilayer,
)

# Multilayers
from .multilayers import (
    generate_twisted_trilayer,
    generate_twisted_multilayer,
)

# Intercalation
from .intercalation import (
    generate_intercalated_bilayer,
    INTERCALANTS,
)

# Stacking
from .stacking import (
    generate_stacked_layers,
    STACKING_SEQUENCES,
)

# Ferroelectric
from .ferroelectric import (
    generate_ferroelectric_2d,
    FERROELECTRIC_2D_MATERIALS,
)


__all__ = [
    # Base
    "TwistStructure", "calculate_moire_angle", "find_commensurate_cell",
    # Moiré
    "generate_moire_superlattice", "calculate_moire_periodicity", "MAGIC_ANGLES",
    # Bilayers
    "generate_twisted_bilayer", "generate_twisted_heterobilayer",
    # Multilayers
    "generate_twisted_trilayer", "generate_twisted_multilayer",
    # Intercalation
    "generate_intercalated_bilayer", "INTERCALANTS",
    # Stacking
    "generate_stacked_layers", "STACKING_SEQUENCES",
    # Ferroelectric
    "generate_ferroelectric_2d", "FERROELECTRIC_2D_MATERIALS",
]
