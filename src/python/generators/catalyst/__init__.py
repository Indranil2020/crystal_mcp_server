"""
catalyst/ - Catalysis & Electrocatalysis

Comprehensive module for catalyst structures.
Covers Category 10 of new_structure_catalogue.md (8 items).

Submodules:
- saa: Single Atom Alloys, HEA surfaces, 2D SAC, strained catalysts
- clusters: Oxide-supported clusters, zeolite-confined, photocatalysts
"""

from .saa import (
    generate_saa_surface,
    generate_hea_surface,
    generate_sac_2d_support,
    generate_strained_catalyst,
    SAA_DATABASE,
    SAC_COORDINATION_MOTIFS,
    HEA_SURFACES,
)

from .clusters import (
    generate_cluster,
    generate_supported_cluster,
    generate_zeolite_cluster,
    generate_photocatalyst_heterojunction,
    CLUSTER_DATABASE,
    OXIDE_SUPPORTS,
    ZEOLITE_FRAMEWORKS,
    HETEROJUNCTIONS,
)


__all__ = [
    # SAA module
    "generate_saa_surface", "generate_hea_surface", 
    "generate_sac_2d_support", "generate_strained_catalyst",
    "SAA_DATABASE", "SAC_COORDINATION_MOTIFS", "HEA_SURFACES",
    # Clusters module
    "generate_cluster", "generate_supported_cluster",
    "generate_zeolite_cluster", "generate_photocatalyst_heterojunction",
    "CLUSTER_DATABASE", "OXIDE_SUPPORTS", "ZEOLITE_FRAMEWORKS", "HETEROJUNCTIONS",
]
