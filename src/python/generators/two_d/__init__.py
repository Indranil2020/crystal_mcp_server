"""
two_d/ - Two-Dimensional Materials

Covers Category 3 of structure_catalogue.md.
"""

import importlib.util

def _module_available(module_path: str) -> bool:
    return importlib.util.find_spec(module_path) is not None

# Heterostructures
if _module_available("generators.two_d.he_2d"):
    from .he_2d import (
        generate_vertical_heterostructure,
        generate_lateral_heterostructure,
        generate_encapsulated_structure,
        MATERIAL_2D_DATABASE,
        HETEROSTRUCTURE_SYSTEMS,
        INTERLAYER_DISTANCES,
    )
else:
    generate_vertical_heterostructure = None
    generate_lateral_heterostructure = None
    generate_encapsulated_structure = None
    MATERIAL_2D_DATABASE = {}
    HETEROSTRUCTURE_SYSTEMS = {}
    INTERLAYER_DISTANCES = {}

# Janus
if _module_available("generators.two_d.janus"):
    from .janus import generate_janus_tmd, JANUS_MATERIALS
else:
    generate_janus_tmd = None
    JANUS_MATERIALS = {}

# Magnetic 2D
if _module_available("generators.two_d.magnetic_2d"):
    from .magnetic_2d import generate_magnetic_2d, MAGNETIC_2D_DATABASE
else:
    generate_magnetic_2d = None
    MAGNETIC_2D_DATABASE = {}

# MXenes
if _module_available("generators.two_d.mxenes"):
    from .mxenes import generate_mxene, MXENE_DATABASE
else:
    generate_mxene = None
    MXENE_DATABASE = {}

# Base
if _module_available("generators.two_d.base"):
    from .base import TWOD_MATERIALS
else:
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
