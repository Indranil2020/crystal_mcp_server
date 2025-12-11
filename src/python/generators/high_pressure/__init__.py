"""
high_pressure/ - High Pressure and Extreme Conditions

Comprehensive module for high-pressure structures.
Covers Category 11 of structure_catalogue.md (5 items).
"""

from .phases import (
    generate_high_pressure_phase,
    generate_polymeric_nitrogen,
    generate_super_ionic_phase,
    generate_superhydride,
    HP_PHASE_DATABASE,
    SUPERHYDRIDE_DATABASE,
)


__all__ = [
    "generate_high_pressure_phase",
    "generate_polymeric_nitrogen",
    "generate_super_ionic_phase",
    "generate_superhydride",
    "HP_PHASE_DATABASE",
    "SUPERHYDRIDE_DATABASE",
]
