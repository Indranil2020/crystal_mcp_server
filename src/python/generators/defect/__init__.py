"""
defect/ - Defect and Disorder Structures

Comprehensive module for defect structures.
Covers Category 7 of structure_catalogue.md (8 items).

Submodules:
- point_defects: Vacancies, interstitials, substitutions
- extended_defects: Dislocations, grain boundaries
- amorphous: Amorphous/glassy structures
- volume_defects: Voids, precipitates, cracks
- sqs: Special Quasirandom Structures, chemical disorder
"""

# Point defects
from .point_defects import (
    generate_vacancy,
    generate_interstitial,
    generate_substitution,
)

# Extended defects
from .extended_defects import (
    generate_dislocation,
    generate_grain_boundary,
)

# Amorphous structures
from .amorphous import (
    generate_amorphous_structure,
    generate_melt_quenched_amorphous,
    AMORPHOUS_DATABASE,
)

# Volume defects
from .volume_defects import (
    generate_void,
    generate_precipitate,
    generate_crack_tip,
    generate_radiation_cascade,
)

# SQS and chemical disorder
from .sqs import (
    generate_sqs,
    generate_sro_structure,
    generate_hea_sqs,
    SQS_ALLOY_DATABASE,
)


__all__ = [
    # Point defects
    "generate_vacancy", "generate_interstitial", "generate_substitution",
    # Extended defects
    "generate_dislocation", "generate_grain_boundary",
    # Amorphous
    "generate_amorphous_structure", "generate_melt_quenched_amorphous", "AMORPHOUS_DATABASE",
    # Volume defects
    "generate_void", "generate_precipitate", "generate_crack_tip", "generate_radiation_cascade",
    # SQS
    "generate_sqs", "generate_sro_structure", "generate_hea_sqs", "SQS_ALLOY_DATABASE",
]
