"""
bulk/ - Bulk Structure Generation

Covers Category 1 and related items in structure_catalogue.md.
"""

# Import from submodules that exist and have the functions

try:
    from .supercells import (
        make_supercell,
        generate_random_alloy_supercell,
        generate_graded_supercell,
        calculate_supercell_for_size,
        generate_slab_supercell,
        COMMON_SUPERCELLS,
        SPECIAL_CELLS,
    )
except ImportError:
    pass

try:
    from .spacegroups import (
        generate_from_spacegroup,
        generate_prototype,
        get_spacegroup_info,
        SPACEGROUP_DATABASE,
        PROTOTYPE_STRUCTURES,
    )
except ImportError:
    pass

try:
    from .zeolites import (
        generate_zeolite,
        generate_zeolite_with_guest,
        ion_exchange_zeolite,
        ZEOLITE_DATABASE,
    )
except ImportError:
    pass

try:
    from .clathrates import (
        generate_clathrate,
        generate_empty_clathrate,
        CLATHRATE_DATABASE,
        CAGE_TYPES,
    )
except ImportError:
    pass

try:
    from .strain import apply_strain, generate_strained_structure
except ImportError:
    pass

try:
    from .polytypes import generate_polytype, POLYTYPE_DATABASE
except ImportError:
    pass

try:
    from .cuprates import generate_cuprate, CUPRATE_DATABASE
except ImportError:
    pass

try:
    from .magnetic import generate_magnetic_bulk
except ImportError:
    pass

try:
    from .base import LATTICE_DATABASE
except ImportError:
    LATTICE_DATABASE = {}


__all__ = [
    # Supercells
    "make_supercell", "generate_random_alloy_supercell", "generate_graded_supercell",
    "calculate_supercell_for_size", "generate_slab_supercell",
    "COMMON_SUPERCELLS", "SPECIAL_CELLS",
    # Space groups
    "generate_from_spacegroup", "generate_prototype", "get_spacegroup_info",
    "SPACEGROUP_DATABASE", "PROTOTYPE_STRUCTURES",
    # Zeolites
    "generate_zeolite", "generate_zeolite_with_guest", "ion_exchange_zeolite",
    "ZEOLITE_DATABASE",
    # Clathrates
    "generate_clathrate", "generate_empty_clathrate",
    "CLATHRATE_DATABASE", "CAGE_TYPES",
    # Other
    "apply_strain", "generate_strained_structure",
    "generate_polytype", "POLYTYPE_DATABASE",
    "generate_cuprate", "CUPRATE_DATABASE",
    "generate_magnetic_bulk",
    "LATTICE_DATABASE",
]
