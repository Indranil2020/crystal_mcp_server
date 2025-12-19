"""
defect/ - Crystal Defects

Covers Category 7 of structure_catalogue.md.
"""

from .point_defects import (
    generate_vacancy,
    generate_interstitial,
    generate_substitution,
    generate_antisite,
    generate_f_center,
    DEFECT_FORMATION_ENERGIES,
    DOPANT_DATABASE,
)

from .extended_defects import (
    generate_dislocation,
    generate_stacking_fault,
    generate_grain_boundary,
    DISLOCATION_DATABASE,
    GRAIN_BOUNDARY_DATABASE,
    STACKING_FAULT_ENERGIES,
)

from .amorphous import (
    generate_amorphous_structure,
    generate_melt_quench,
    generate_amorphous_interface,
    AMORPHOUS_DATABASE,
    COORDINATION_NUMBERS,
)

from .volume_defects import (
    generate_void,
    generate_precipitate,
)

from .sqs import (
    generate_sqs,
)


__all__ = [
    # Point defects
    "generate_vacancy", "generate_interstitial", "generate_substitution",
    "generate_antisite", "generate_f_center",
    "DEFECT_FORMATION_ENERGIES", "DOPANT_DATABASE",
    # Extended defects
    "generate_dislocation", "generate_stacking_fault", "generate_grain_boundary",
    "DISLOCATION_DATABASE", "GRAIN_BOUNDARY_DATABASE", "STACKING_FAULT_ENERGIES",
    # Amorphous
    "generate_amorphous_structure", "generate_melt_quench", "generate_amorphous_interface",
    "AMORPHOUS_DATABASE", "COORDINATION_NUMBERS",
    # Volume defects
    "generate_void", "generate_precipitate",
    # SQS
    "generate_sqs",
]
