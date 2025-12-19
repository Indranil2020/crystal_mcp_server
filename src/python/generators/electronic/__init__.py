"""
electronic/ - Electronic Materials

Comprehensive module for electronic materials.
Covers Category 8 (partial) of new_structure_catalogue.md.

Submodules:
- topological: TIs, Weyl/Dirac semimetals, Majorana platforms
- semiconductors: Semiconductor heterostructures
- superconductors: Superconducting materials
"""

from .topological import (
    generate_topological_insulator,
    generate_weyl_semimetal,
    generate_majorana_nanowire,
    TOPOLOGICAL_INSULATOR_DATABASE,
    SEMIMETAL_DATABASE,
    MAJORANA_DATABASE,
)


__all__ = [
    # Topological
    "generate_topological_insulator", "generate_weyl_semimetal", "generate_majorana_nanowire",
    "TOPOLOGICAL_INSULATOR_DATABASE", "SEMIMETAL_DATABASE", "MAJORANA_DATABASE",
]
