"""
external_fields/ - External Field and Driven Structures

Comprehensive module for field-driven structures.
Covers Category 14 of structure_catalogue.md (5 items).
"""

from .magnetic_textures import (
    generate_skyrmion,
    generate_antiskyrmion,
    generate_magnetic_bobber,
    SKYRMION_DATABASE,
)

from .strain_fields import (
    generate_strained_structure,
    generate_rippled_structure,
    apply_strain_field,
    STRAIN_PROFILES,
)


__all__ = [
    "generate_skyrmion", "generate_antiskyrmion", "generate_magnetic_bobber",
    "SKYRMION_DATABASE",
    "generate_strained_structure", "generate_rippled_structure", "apply_strain_field",
    "STRAIN_PROFILES",
]
