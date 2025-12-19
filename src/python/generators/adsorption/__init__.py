"""
adsorption/ - Adsorption and Interface Structures

Covers Category 6 of structure_catalogue.md.
"""

from .adsorbates import (
    generate_adsorbate_on_surface,
    generate_coadsorption,
    get_adsorbate_molecule,
    ADSORBATE_DATABASE,
    COVERAGE_PATTERNS,
    ADSORPTION_SITES,
)

from .sac import (
    generate_sac_graphene,
    generate_sac_oxide,
    generate_sac_tmd,
    SAC_METALS,
    SAC_SUPPORTS,
    COORDINATION_MOTIFS,
)


__all__ = [
    # Adsorbates
    "generate_adsorbate_on_surface", "generate_coadsorption", "get_adsorbate_molecule",
    "ADSORBATE_DATABASE", "COVERAGE_PATTERNS", "ADSORPTION_SITES",
    # SAC
    "generate_sac_graphene", "generate_sac_oxide", "generate_sac_tmd",
    "SAC_METALS", "SAC_SUPPORTS", "COORDINATION_MOTIFS",
]
