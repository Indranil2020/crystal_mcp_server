"""
surface/ - Surface and Slab Structure Generation

Comprehensive module for surface structures.
Covers Category 4 of structure_catalogue.md (8 items).

Submodules:
- slabs: Miller-index surfaces, polar slabs
- reconstructions: Si(100)-2×1, Au herringbone, etc.
- nanoparticles: Wulff shapes, core-shell
- adatoms: Adatom/advacancy superlattices
- stepped: Vicinal surfaces, kinks, terraces
- alloys: Surface alloying, de-alloying
- interfaces: Polar interfaces (LAO/STO), oxide superlattices
"""

# Slabs and polar surfaces
from .slabs import generate_slab, generate_polar_slab

# Surface reconstructions
from .reconstructions import generate_reconstruction, RECONSTRUCTIONS

# Nanoparticles
from .nanoparticles import generate_wulff_nanoparticle, generate_core_shell

# Adatom superlattices
from .adatoms import (
    generate_adatom_superlattice,
    generate_advacancy_superlattice,
    ADATOM_PATTERNS,
    ADATOM_SYSTEMS,
)

# Stepped and vicinal surfaces
from .stepped import (
    generate_vicinal_surface,
    generate_stepped_surface_with_kinks,
    generate_roughened_surface,
    VICINAL_SURFACES,
)

# Surface alloying
from .alloys import (
    generate_surface_alloy,
    generate_dealloyed_surface,
    generate_segregation_profile,
    SURFACE_ALLOY_DATABASE,
)

# Polar interfaces
from .interfaces import (
    generate_polar_interface,
    generate_oxide_superlattice,
    generate_interface_with_defects,
    POLAR_INTERFACES,
    PEROVSKITE_DATABASE,
)


__all__ = [
    # Slabs
    "generate_slab", "generate_polar_slab",
    # Reconstructions
    "generate_reconstruction", "RECONSTRUCTIONS",
    # Nanoparticles
    "generate_wulff_nanoparticle", "generate_core_shell",
    # Adatoms
    "generate_adatom_superlattice", "generate_advacancy_superlattice",
    "ADATOM_PATTERNS", "ADATOM_SYSTEMS",
    # Stepped
    "generate_vicinal_surface", "generate_stepped_surface_with_kinks",
    "generate_roughened_surface", "VICINAL_SURFACES",
    # Alloys
    "generate_surface_alloy", "generate_dealloyed_surface",
    "generate_segregation_profile", "SURFACE_ALLOY_DATABASE",
    # Interfaces
    "generate_polar_interface", "generate_oxide_superlattice",
    "generate_interface_with_defects", "POLAR_INTERFACES", "PEROVSKITE_DATABASE",
]
