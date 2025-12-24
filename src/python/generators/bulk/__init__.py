"""
bulk/ - Bulk Structure Generation

Covers Category 1 and related items in structure_catalogue.md.
"""

# Import from submodules that exist and have the functions

import importlib.util


def _module_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


_BASE_MODULE = __name__

if _module_available(f"{_BASE_MODULE}.supercells"):
    from .supercells import (
        make_supercell,
        generate_random_alloy_supercell,
        generate_graded_supercell,
        calculate_supercell_for_size,
        generate_slab_supercell,
        COMMON_SUPERCELLS,
        SPECIAL_CELLS,
    )
else:
    make_supercell = None
    generate_random_alloy_supercell = None
    generate_graded_supercell = None
    calculate_supercell_for_size = None
    generate_slab_supercell = None
    COMMON_SUPERCELLS = {}
    SPECIAL_CELLS = {}

if _module_available(f"{_BASE_MODULE}.spacegroups"):
    from .spacegroups import (
        generate_from_spacegroup,
        generate_prototype,
        get_spacegroup_info,
        SPACEGROUP_DATABASE,
        PROTOTYPE_STRUCTURES,
    )
else:
    generate_from_spacegroup = None
    generate_prototype = None
    get_spacegroup_info = None
    SPACEGROUP_DATABASE = {}
    PROTOTYPE_STRUCTURES = {}

if _module_available(f"{_BASE_MODULE}.zeolites"):
    from .zeolites import (
        generate_zeolite,
        generate_zeolite_with_guest,
        ion_exchange_zeolite,
        ZEOLITE_DATABASE,
    )
else:
    generate_zeolite = None
    generate_zeolite_with_guest = None
    ion_exchange_zeolite = None
    ZEOLITE_DATABASE = {}

if _module_available(f"{_BASE_MODULE}.clathrates"):
    from .clathrates import (
        generate_clathrate,
        generate_empty_clathrate,
        CLATHRATE_DATABASE,
        CAGE_TYPES,
    )
else:
    generate_clathrate = None
    generate_empty_clathrate = None
    CLATHRATE_DATABASE = {}
    CAGE_TYPES = {}

if _module_available(f"{_BASE_MODULE}.strain"):
    from .strain import apply_strain, generate_strained_structure
else:
    apply_strain = None
    generate_strained_structure = None

if _module_available(f"{_BASE_MODULE}.polytypes"):
    from .polytypes import generate_polytype, POLYTYPE_DATABASE
else:
    generate_polytype = None
    POLYTYPE_DATABASE = {}

if _module_available(f"{_BASE_MODULE}.cuprates"):
    from .cuprates import generate_cuprate, CUPRATE_DATABASE
else:
    generate_cuprate = None
    CUPRATE_DATABASE = {}

if _module_available(f"{_BASE_MODULE}.magnetic"):
    from .magnetic import generate_magnetic_bulk
else:
    generate_magnetic_bulk = None

if _module_available(f"{_BASE_MODULE}.base"):
    from .base import LATTICE_DATABASE
else:
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
