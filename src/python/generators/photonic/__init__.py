"""
photonic/ - Photonic and Metamaterial Structures

Comprehensive module for photonic structures.
Covers Category 13 of new_structure_catalogue.md.

Submodules:
- metamaterials: Metasurfaces, photonic crystals
- crystals: Additional photonic crystal types
"""

from .metamaterials import (
    generate_metasurface_unit_cell,
    generate_photonic_crystal,
    generate_broken_symmetry_photonic,
    METAMATERIAL_DATABASE,
    PHOTONIC_CRYSTAL_DATABASE,
)


__all__ = [
    "generate_metasurface_unit_cell", "generate_photonic_crystal",
    "generate_broken_symmetry_photonic",
    "METAMATERIAL_DATABASE", "PHOTONIC_CRYSTAL_DATABASE",
]
