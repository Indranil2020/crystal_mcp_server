"""
surface/ - Surface and Slab Structures

Covers Category 4 from structure_catalogue.md.
"""

from .slabs import (
    generate_slab,
    generate_vicinal_surface,
    generate_kinked_surface,
    generate_surface_with_adterrace,
    SURFACE_ENERGY_DATABASE,
    SLAB_CONFIGURATIONS,
)

from .reconstructions import (
    generate_reconstruction,
)


__all__ = [
    "generate_slab", "generate_vicinal_surface", "generate_kinked_surface",
    "generate_surface_with_adterrace", "SURFACE_ENERGY_DATABASE", "SLAB_CONFIGURATIONS",
    "generate_reconstruction",
]
