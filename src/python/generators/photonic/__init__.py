"""
photonic/ - Photonic Material Generation

Category 11: Photonic crystals, metamaterials.
"""

from .crystals import generate_photonic_crystal, PHOTONIC_LATTICES
from .metamaterials import generate_metamaterial


__all__ = ["generate_photonic_crystal", "PHOTONIC_LATTICES", "generate_metamaterial"]
