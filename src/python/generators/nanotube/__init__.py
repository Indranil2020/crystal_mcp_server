"""
nanotube/ - Nanotube and Nanowire Structures

Comprehensive module for 1D nanostructures.
Covers Category 8 of structure_catalogue.md (9 items).

Submodules:
- cnt: Carbon nanotubes (SWCNTs, MWCNTs, functionalized)
- nanowires: Semiconductor and metal nanowires
- other_tubes: BN nanotubes, TMD nanotubes
"""

from .cnt import (
    generate_cnt,
    generate_mwcnt,
    generate_functionalized_cnt,
    generate_cnt_bundle,
    generate_doped_cnt,
    calculate_cnt_properties,
    CNT_DATABASE,
    FUNCTIONAL_GROUPS,
)

from .nanowires import (
    generate_nanowire,
    generate_coreshell_nanowire,
    generate_axial_heterostructure,
    NANOWIRE_DATABASE,
    CORESHELL_DATABASE,
)


__all__ = [
    # CNT
    "generate_cnt", "generate_mwcnt", "generate_functionalized_cnt",
    "generate_cnt_bundle", "generate_doped_cnt", "calculate_cnt_properties",
    "CNT_DATABASE", "FUNCTIONAL_GROUPS",
    # Nanowires
    "generate_nanowire", "generate_coreshell_nanowire", "generate_axial_heterostructure",
    "NANOWIRE_DATABASE", "CORESHELL_DATABASE",
]
