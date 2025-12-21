"""
meta_structures - Meta-structure and Advanced Lattice Generation

This module provides generators for:
- Quasicrystals (Penrose tilings, Fibonacci approximants)
- Artificial lattices (Kagome, Lieb, checkerboard)
- Moiré superlattices
- Photonic crystals
"""

from .quasicrystals import (
    generate_quasicrystal,
    generate_penrose_tiling,
    generate_fibonacci_approximant,
    generate_icosahedral_quasicrystal
)

from .artificial_lattices import (
    generate_kagome_lattice,
    generate_lieb_lattice,
    generate_checkerboard_lattice,
    generate_dice_lattice,
    generate_honeycomb_variants
)

__all__ = [
    'generate_quasicrystal',
    'generate_penrose_tiling',
    'generate_fibonacci_approximant',
    'generate_icosahedral_quasicrystal',
    'generate_kagome_lattice',
    'generate_lieb_lattice',
    'generate_checkerboard_lattice',
    'generate_dice_lattice',
    'generate_honeycomb_variants'
]
