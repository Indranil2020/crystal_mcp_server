"""
magnetic/ - Magnetic Materials Generation

Category 14: Ferromagnets, antiferromagnets, spintronic materials.
"""

from .materials import generate_magnetic_material, MAGNETIC_DATABASE
from .heusler import generate_heusler, HEUSLER_DATABASE


__all__ = ["generate_magnetic_material", "MAGNETIC_DATABASE", "generate_heusler", "HEUSLER_DATABASE"]
