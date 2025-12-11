"""
two_d/ - 2D Materials and Exfoliable Slabs

Comprehensive module for 2D material structure generation.
Covers all 11 structure types from Category 2 of structure_catalogue.md.

Submodules:
- layergroups: 80 layer groups, slab representations
- xenes: Graphene, silicene, germanene, phosphorene, etc.
- tmds: MoS2, WS2, with 1T, 2H, 1T' phases
- mxenes: Ti3C2Tx with terminations
- janus: Broken inversion (MoSSe, WSSe)
- magnetic_2d: CrI3, CrGeTe3, Fe3GeTe2
- he_2d: High-entropy 2D materials
- ribbons: Armchair, zigzag, chiral nanoribbons
- porous: g-C3N4, holey graphene, biphenylene
- flakes: Nanomeshes, quantum dots
- grain_boundaries: 2D grain boundaries
"""

# Base utilities
from .base import TwoD_Structure, structure_to_dict, atoms_to_dict

# X-enes
from .xenes import (
    generate_xene,
    XENE_PARAMS,
)

# TMDs
from .tmds import (
    generate_tmd,
    TMD_DATABASE,
)

# MXenes
from .mxenes import (
    generate_mxene,
    MXENE_DATABASE,
)

# Janus
from .janus import (
    generate_janus_2d,
    JANUS_MATERIALS,
)

# Magnetic 2D
from .magnetic_2d import (
    generate_magnetic_2d,
    MAGNETIC_2D_MATERIALS,
)

# High-entropy 2D
from .he_2d import (
    generate_he_2d,
)

# Nanoribbons
from .ribbons import (
    generate_nanoribbon,
    EDGE_TYPES,
)

# Porous 2D
from .porous import (
    generate_porous_2d,
    POROUS_2D_MATERIALS,
)

# Flakes and QDs
from .flakes import (
    generate_nanoflake,
    generate_quantum_dot,
)

# Grain boundaries
from .grain_boundaries import (
    generate_2d_grain_boundary,
)


__all__ = [
    # Base
    "TwoD_Structure", "structure_to_dict", "atoms_to_dict",
    # X-enes
    "generate_xene", "XENE_PARAMS",
    # TMDs
    "generate_tmd", "TMD_DATABASE",
    # MXenes
    "generate_mxene", "MXENE_DATABASE",
    # Janus
    "generate_janus_2d", "JANUS_MATERIALS",
    # Magnetic 2D
    "generate_magnetic_2d", "MAGNETIC_2D_MATERIALS",
    # High-entropy 2D
    "generate_he_2d",
    # Ribbons
    "generate_nanoribbon", "EDGE_TYPES",
    # Porous
    "generate_porous_2d", "POROUS_2D_MATERIALS",
    # Flakes
    "generate_nanoflake", "generate_quantum_dot",
    # Grain boundaries
    "generate_2d_grain_boundary",
]
