"""
magnetic/ - Magnetic Materials

Covers magnetic structure generation from structure_catalogue.md.
"""

from .heusler import (
    generate_full_heusler,
    generate_half_heusler,
    generate_inverse_heusler,
    generate_quaternary_heusler,
    FULL_HEUSLER_DATABASE,
    HALF_HEUSLER_DATABASE,
    QUATERNARY_HEUSLER_DATABASE,
)


__all__ = [
    "generate_full_heusler", "generate_half_heusler",
    "generate_inverse_heusler", "generate_quaternary_heusler",
    "FULL_HEUSLER_DATABASE", "HALF_HEUSLER_DATABASE", "QUATERNARY_HEUSLER_DATABASE",
]
