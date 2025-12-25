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
    from .janus import generate_janus_2d, list_janus_materials, JANUS_MATERIALS
    # Alias for backwards compatibility
    generate_janus_tmd = generate_janus_2d
else:
    generate_janus_2d = None
    generate_janus_tmd = None
    list_janus_materials = None
    JANUS_MATERIALS = {}

# Magnetic 2D
if _module_available("generators.two_d.magnetic_2d"):
    from .magnetic_2d import generate_magnetic_2d, MAGNETIC_2D_MATERIALS
    MAGNETIC_2D_DATABASE = MAGNETIC_2D_MATERIALS  # Alias
else:
    generate_magnetic_2d = None
    MAGNETIC_2D_DATABASE = {}
    MAGNETIC_2D_MATERIALS = {}

# MXenes
if _module_available("generators.two_d.mxenes"):
    from .mxenes import generate_mxene, MXENE_DATABASE
else:
    generate_mxene = None
    MXENE_DATABASE = {}

# Base utilities are always available
from .base import structure_to_dict, atoms_to_dict, add_vacuum, center_in_vacuum

# Placeholder for 2D materials database
TWOD_MATERIALS = {}


__all__ = [
    # Heterostructures
    "generate_vertical_heterostructure", "generate_lateral_heterostructure",
    "generate_encapsulated_structure", "MATERIAL_2D_DATABASE",
    "HETEROSTRUCTURE_SYSTEMS", "INTERLAYER_DISTANCES",
    # Janus
    "generate_janus_2d", "generate_janus_tmd", "list_janus_materials", "JANUS_MATERIALS",
    # Magnetic 2D
    "generate_magnetic_2d", "MAGNETIC_2D_DATABASE", "MAGNETIC_2D_MATERIALS",
    # MXenes
    "generate_mxene", "MXENE_DATABASE",
    # Base utilities
    "structure_to_dict", "atoms_to_dict", "add_vacuum", "center_in_vacuum",
    "TWOD_MATERIALS",
]
