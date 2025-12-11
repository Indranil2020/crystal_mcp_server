"""
battery/ - Battery Materials

Comprehensive module for battery materials.
Covers Category 9 of structure_catalogue.md (8 items).

Submodules:
- cathodes: Layered, olivine, spinel, DRX cathodes
- electrolytes: Solid state electrolytes (garnet, sulfide, NASICON)
- anodes: Graphite, Si, Li metal anodes
"""

from .cathodes import (
    generate_cathode,
    generate_lithiated_series,
    generate_surface_coated_cathode,
    CATHODE_DATABASE,
)

from .electrolytes import (
    generate_electrolyte,
    generate_cathode_electrolyte_interface,
    compare_electrolytes,
    ELECTROLYTE_DATABASE,
)


__all__ = [
    # Cathodes
    "generate_cathode", "generate_lithiated_series", "generate_surface_coated_cathode",
    "CATHODE_DATABASE",
    # Electrolytes
    "generate_electrolyte", "generate_cathode_electrolyte_interface", "compare_electrolytes",
    "ELECTROLYTE_DATABASE",
]
