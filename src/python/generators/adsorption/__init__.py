"""
adsorption/ - Adsorption and Interface Motifs

Comprehensive module for adsorption structures.
Covers Category 6 of structure_catalogue.md (9 items).

Submodules:
- sac: Single Atom Catalysts
- adsorbates: Simple molecular adsorbates
- fragments: Dissociated fragments (OH*, CHx*, etc.)
- water_layers: Water bilayers and ice
- sams: Self-Assembled Monolayers
- ionic_liquids: IL monolayers
"""

# Single Atom Catalysts
from .sac import generate_sac, SAC_SUPPORTS

# Simple adsorbates
from .adsorbates import add_adsorbate, ADSORPTION_SITES

# Dissociated fragments
from .fragments import (
    generate_fragment_adsorbate,
    generate_coverage_pattern,
    FRAGMENT_DATABASE,
)

# Water and ice layers
from .water_layers import (
    generate_water_bilayer,
    generate_ice_layer,
    WATER_LAYER_DATABASE,
)

# Self-Assembled Monolayers
from .sams import (
    generate_sam,
    generate_mixed_sam,
    SAM_DATABASE,
)

# Ionic liquid monolayers
from .ionic_liquids import (
    generate_il_monolayer,
    IL_DATABASE,
)


__all__ = [
    # SAC
    "generate_sac", "SAC_SUPPORTS",
    # Adsorbates
    "add_adsorbate", "ADSORPTION_SITES",
    # Fragments
    "generate_fragment_adsorbate", "generate_coverage_pattern", "FRAGMENT_DATABASE",
    # Water
    "generate_water_bilayer", "generate_ice_layer", "WATER_LAYER_DATABASE",
    # SAMs
    "generate_sam", "generate_mixed_sam", "SAM_DATABASE",
    # Ionic liquids
    "generate_il_monolayer", "IL_DATABASE",
]
