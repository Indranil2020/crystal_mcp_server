"""
quality_control/ - Quality Control and Metadata

Comprehensive module for structure validation and analysis.
Covers Category 16 of structure_catalogue.md (7 items).
"""

from .symmetry import (
    analyze_symmetry,
    find_space_group,
    check_tolerance,
)

from .validators import (
    check_stoichiometry,
    check_bond_valence,
    check_charge_neutrality,
    remove_duplicates,
)

from .kpath import (
    generate_kpath,
    get_high_symmetry_points,
    calculate_kmesh,
    HIGH_SYMMETRY_POINTS,
)


__all__ = [
    "analyze_symmetry", "find_space_group", "check_tolerance",
    "check_stoichiometry", "check_bond_valence", "check_charge_neutrality", "remove_duplicates",
    "generate_kpath", "get_high_symmetry_points", "calculate_kmesh", "HIGH_SYMMETRY_POINTS",
]
