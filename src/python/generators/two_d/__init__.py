"""
two_d/ - Two-Dimensional Materials

Covers Category 3 of structure_catalogue.md.
"""

try:
    from .he_2d import (
        generate_vertical_heterostructure,
        generate_lateral_heterostructure,
        generate_encapsulated_structure,
        MATERIAL_2D_DATABASE,
        HETEROSTRUCTURE_SYSTEMS,
        INTERLAYER_DISTANCES,
    )
except ImportError:
    pass

try:
    from .janus import generate_janus_tmd, JANUS_MATERIALS
except ImportError:
    pass

try:
    from .magnetic_2d import generate_magnetic_2d, MAGNETIC_2D_DATABASE
except ImportError:
    pass

try:
    from .mxenes import generate_mxene, MXENE_DATABASE
except ImportError:
    pass

try:
    from .base import TWOD_MATERIALS
except ImportError:
    TWOD_MATERIALS = {}


__all__ = [
    # Heterostructures
    "generate_vertical_heterostructure", "generate_lateral_heterostructure",
    "generate_encapsulated_structure", "MATERIAL_2D_DATABASE",
    "HETEROSTRUCTURE_SYSTEMS", "INTERLAYER_DISTANCES",
    # Janus
    "generate_janus_tmd", "JANUS_MATERIALS",
    # Magnetic 2D
    "generate_magnetic_2d", "MAGNETIC_2D_DATABASE",
    # MXenes
    "generate_mxene", "MXENE_DATABASE",
    # Base
    "TWOD_MATERIALS",
]
