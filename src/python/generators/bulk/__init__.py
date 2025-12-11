"""
bulk/ - Crystalline Bulk Structure Generators

Comprehensive module for 3D periodic crystal structure generation.
Covers all 22 structure types from Category 1 of structure_catalogue.md.

Submodules:
- spacegroups: 230 space groups, Wyckoff positions
- supercells: Integer/rational matrices, commensurate cells
- polytypes: 2H, 4H, 6H, 15R stacking sequences
- magnetic: FM, AFM, FiM, spiral orderings
- perovskites: Simple, double, layered perovskites
- clathrates: Type I, II, III cage structures
- disorder: SQS, random alloys
- boundaries: Twin, APB, stacking faults
- defects: Vacancy, interstitial, antisite
- zeolites: IZA framework topologies
- cuprates: YBCO, LSCO, Bi-2212
- strain: Bain path, epitaxial strain
- special_phases: Superhydrides, NTE, Prussian blue
"""

# Base utilities
from .base import BulkStructure, structure_to_dict, atoms_to_dict

# Spacegroup generation
from .spacegroups import (
    generate_crystal,
    generate_all_polymorphs,
    get_wyckoff_positions,
)

# Supercells
from .supercells import (
    make_supercell,
    make_rational_supercell,
    make_commensurate_cell,
    SUPERCELL_PRESETS,
)

# Polytypes
from .polytypes import (
    generate_polytype,
    generate_custom_stacking,
    POLYTYPE_SEQUENCES,
)

# Magnetic orderings
from .magnetic import (
    generate_magnetic_ordering,
    get_magnetic_space_group,
    MAGNETIC_ORDERINGS,
)

# Perovskites
from .perovskites import (
    generate_perovskite,
    generate_double_perovskite,
    generate_layered_perovskite,
)

# Clathrates
from .clathrates import (
    generate_clathrate,
    CLATHRATE_TYPES,
)

# Disorder and alloys
from .disorder import (
    generate_sqs,
    generate_random_alloy,
    calculate_short_range_order,
)

# Boundaries
from .boundaries import (
    generate_twin_boundary,
    generate_antiphase_boundary,
    generate_stacking_fault,
)

# Point defects
from .defects import (
    generate_vacancy,
    generate_interstitial,
    generate_antisite,
    generate_dopant_distribution,
)

# Zeolites
from .zeolites import (
    generate_zeolite,
    get_zeolite_info,
    list_zeolite_topologies,
    ZEOLITE_DATABASE,
)

# Cuprates
from .cuprates import (
    generate_cuprate,
    generate_custom_cuprate,
    CUPRATE_DATABASE,
)

# Strain and transformations
from .strain import (
    apply_strain,
    apply_strain_tensor,
    generate_bain_path,
    generate_epitaxial_strain,
)

# Special phases
from .special_phases import (
    generate_superhydride,
    generate_high_pressure_phase,
    generate_nte_material,
    generate_prussian_blue,
    HIGH_PRESSURE_PHASES,
)

# Dislocations and extended defects
from .dislocations import (
    generate_elastic_dipole,
    generate_dislocation_dipole,
    generate_frank_loop,
    generate_prismatic_loop,
    generate_sft,
    DISLOCATION_TYPES,
)

# Quantum condensed phases
from .quantum_phases import (
    generate_cdw_phase,
    generate_peierls_distortion,
    generate_jahn_teller,
    generate_breathing_mode,
    generate_spin_peierls,
    CDW_MATERIALS,
)

# Energetic materials
from .energetics import (
    generate_energetic_material,
    generate_energetic_cocrystal,
    generate_porous_energetic,
    generate_pbx_model,
    ENERGEIC_MATERIALS,
)


__all__ = [
    # Base
    "BulkStructure", "structure_to_dict", "atoms_to_dict",
    # Spacegroups
    "generate_crystal", "generate_all_polymorphs", "get_wyckoff_positions",
    # Supercells
    "make_supercell", "make_rational_supercell", "make_commensurate_cell", "SUPERCELL_PRESETS",
    # Polytypes
    "generate_polytype", "generate_custom_stacking", "POLYTYPE_SEQUENCES",
    # Magnetic
    "generate_magnetic_ordering", "get_magnetic_space_group", "MAGNETIC_ORDERINGS",
    # Perovskites
    "generate_perovskite", "generate_double_perovskite", "generate_layered_perovskite",
    # Clathrates
    "generate_clathrate", "CLATHRATE_TYPES",
    # Disorder
    "generate_sqs", "generate_random_alloy", "calculate_short_range_order",
    # Boundaries
    "generate_twin_boundary", "generate_antiphase_boundary", "generate_stacking_fault",
    # Defects
    "generate_vacancy", "generate_interstitial", "generate_antisite", "generate_dopant_distribution",
    # Zeolites
    "generate_zeolite", "get_zeolite_info", "list_zeolite_topologies", "ZEOLITE_DATABASE",
    # Cuprates
    "generate_cuprate", "generate_custom_cuprate", "CUPRATE_DATABASE",
    # Strain
    "apply_strain", "apply_strain_tensor", "generate_bain_path", "generate_epitaxial_strain",
    # Special phases
    "generate_superhydride", "generate_high_pressure_phase", "generate_nte_material",
    "generate_prussian_blue", "HIGH_PRESSURE_PHASES",
    # Dislocations
    "generate_elastic_dipole", "generate_dislocation_dipole", "generate_frank_loop",
    "generate_prismatic_loop", "generate_sft", "DISLOCATION_TYPES",
    # Quantum phases
    "generate_cdw_phase", "generate_peierls_distortion", "generate_jahn_teller",
    "generate_breathing_mode", "generate_spin_peierls", "CDW_MATERIALS",
    # Energetics
    "generate_energetic_material", "generate_energetic_cocrystal",
    "generate_porous_energetic", "generate_pbx_model", "ENERGEIC_MATERIALS",
]
