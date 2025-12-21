"""
meta_structures/artificial_lattices.py - Artificial Lattice Generation

Provides generation of artificial 2D lattices with exotic band structures:
- Kagome lattice (flat bands)
- Lieb lattice (flat bands at Fermi level)
- Checkerboard lattice
- Dice (T3) lattice
- Honeycomb variants (alpha-T3, etc.)
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


def generate_kagome_lattice(
    element: str = "Fe",
    a: float = 5.0,
    size: Tuple[int, int] = (3, 3),
    vacuum: float = 15.0,
    magnetic_order: Optional[str] = None,
    buckling: float = 0.0
) -> Dict[str, Any]:
    """
    Generate Kagome lattice structure.

    The Kagome lattice features flat bands and is studied for:
    - Frustrated magnetism
    - Superconductivity (e.g., AV3Sb5)
    - Topological properties

    Args:
        element: Element for lattice sites
        a: Lattice parameter (Angstroms)
        size: Supercell size (nx, ny)
        vacuum: Vacuum thickness for 2D slab
        magnetic_order: None, "FM", "AFM_120", "spin_liquid"
        buckling: Out-of-plane buckling amplitude

    Returns:
        Kagome lattice structure
    """
    from pymatgen.core import Structure, Lattice

    # Kagome has 3 atoms per unit cell in triangular Bravais lattice
    # Positions: (0,0), (1/2,0), (1/4,1/2)

    # Triangular lattice vectors
    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a * np.sqrt(3)/2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # Kagome basis positions
    kagome_basis = [
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.25, 0.5, 0.5]
    ]

    # Add buckling
    if buckling > 0:
        kagome_basis[0][2] = 0.5 + buckling / vacuum
        kagome_basis[1][2] = 0.5 - buckling / vacuum
        kagome_basis[2][2] = 0.5

    struct = Structure(
        lattice=lattice,
        species=[element] * 3,
        coords=kagome_basis,
        coords_are_cartesian=False
    )

    # Apply supercell
    struct = struct * [size[0], size[1], 1]

    # Add magnetic moments if specified
    magnetic_info = None
    if magnetic_order:
        moments = []
        if magnetic_order == "FM":
            moments = [[0, 0, 2.0]] * len(struct)
        elif magnetic_order == "AFM_120":
            # 120-degree order
            for i, site in enumerate(struct.sites):
                angle = (i % 3) * 2 * np.pi / 3
                moments.append([2.0 * np.cos(angle), 2.0 * np.sin(angle), 0])
        elif magnetic_order == "spin_liquid":
            moments = [[0, 0, 0]] * len(struct)  # Disordered

        magnetic_info = {
            "order": magnetic_order,
            "moments": moments
        }

    return {
        "success": True,
        "lattice_type": "kagome",
        "a": a,
        "n_atoms": len(struct),
        "size": list(size),
        "buckling_A": buckling,
        "band_structure_features": [
            "Flat band at E = 2t",
            "Dirac cone at K point",
            "Topologically non-trivial"
        ],
        "magnetic_order": magnetic_order,
        "magnetic_info": magnetic_info,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_lieb_lattice(
    elements: List[str] = ["Cu", "O"],
    a: float = 3.8,
    size: Tuple[int, int] = (3, 3),
    vacuum: float = 15.0,
    dimerization: float = 0.0
) -> Dict[str, Any]:
    """
    Generate Lieb lattice structure.

    The Lieb lattice has a flat band at E=0 (Fermi level) and
    is relevant for cuprate physics.

    Args:
        elements: [corner, edge] elements (e.g., ["Cu", "O"])
        a: Lattice parameter
        size: Supercell size
        vacuum: Vacuum thickness
        dimerization: Peierls dimerization amplitude

    Returns:
        Lieb lattice structure
    """
    from pymatgen.core import Structure, Lattice

    # Lieb lattice: square with atoms at corners and edge centers
    # 3 atoms per unit cell

    lattice = Lattice.orthorhombic(a, a, vacuum)

    # Basis: corner (0,0), edge-x (1/2,0), edge-y (0,1/2)
    lieb_basis = [
        [0.0, 0.0, 0.5],      # Corner
        [0.5, 0.0, 0.5],      # Edge x
        [0.0, 0.5, 0.5]       # Edge y
    ]

    # Apply dimerization (breathing mode)
    if dimerization > 0:
        delta = dimerization / a
        lieb_basis[1][0] = 0.5 - delta
        lieb_basis[2][1] = 0.5 - delta

    species = [elements[0], elements[1] if len(elements) > 1 else elements[0],
               elements[1] if len(elements) > 1 else elements[0]]

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=lieb_basis,
        coords_are_cartesian=False
    )

    struct = struct * [size[0], size[1], 1]

    return {
        "success": True,
        "lattice_type": "lieb",
        "a": a,
        "n_atoms": len(struct),
        "size": list(size),
        "elements": elements,
        "dimerization": dimerization,
        "band_structure_features": [
            "Flat band at E = 0",
            "Dirac cones at M points",
            "Bipartite lattice with imbalanced sublattices"
        ],
        "applications": [
            "CuO2 planes in cuprates",
            "Photonic lattices",
            "Cold atom systems"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_checkerboard_lattice(
    elements: List[str] = ["A", "B"],
    a: float = 4.0,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    crossing_type: str = "standard"
) -> Dict[str, Any]:
    """
    Generate checkerboard lattice.

    The checkerboard lattice is a decorated square lattice with
    crossed plaquettes.

    Args:
        elements: Elements for the two sublattices
        a: Lattice parameter
        size: Supercell size
        vacuum: Vacuum thickness
        crossing_type: "standard", "pyrochlore_projection"

    Returns:
        Checkerboard lattice structure
    """
    from pymatgen.core import Structure, Lattice

    # Checkerboard: 2 atoms per primitive cell on square lattice
    # with diagonal connections

    lattice = Lattice.orthorhombic(a, a, vacuum)

    # Standard checkerboard
    coords = [
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5]
    ]

    if crossing_type == "pyrochlore_projection":
        # Add extra sites at plaquette centers for pyrochlore-like
        coords.extend([
            [0.25, 0.25, 0.5],
            [0.75, 0.75, 0.5]
        ])
        species = [elements[0], elements[0], elements[1] if len(elements) > 1 else elements[0]] * 2
        species = species[:len(coords)]
    else:
        species = [elements[0], elements[1] if len(elements) > 1 else elements[0]]

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    struct = struct * [size[0], size[1], 1]

    return {
        "success": True,
        "lattice_type": "checkerboard",
        "crossing_type": crossing_type,
        "a": a,
        "n_atoms": len(struct),
        "size": list(size),
        "band_structure_features": [
            "Flat band for pyrochlore projection",
            "Frustrated magnetism possible"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_dice_lattice(
    element: str = "C",
    a: float = 2.5,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    alpha: float = 1.0
) -> Dict[str, Any]:
    """
    Generate dice (T3) lattice.

    The dice lattice is honeycomb with a hub atom at each hexagon center.
    It has a flat band and pseudospin-1 Dirac fermions.

    Args:
        element: Element for all sites
        a: Lattice parameter
        size: Supercell size
        vacuum: Vacuum thickness
        alpha: Hub coupling parameter (1 = full dice, 0 = honeycomb)

    Returns:
        Dice lattice structure
    """
    from pymatgen.core import Structure, Lattice

    # Dice lattice: honeycomb + hub at hexagon centers
    # 3 atoms per unit cell in triangular Bravais

    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a * np.sqrt(3)/2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # Basis: A (0,0), B (1/3,1/3), Hub (2/3,2/3) - or equivalently hexagon center
    dice_basis = [
        [0.0, 0.0, 0.5],      # A sublattice
        [1/3, 1/3, 0.5],      # B sublattice
        [2/3, 2/3, 0.5]       # Hub (center of hexagon)
    ]

    struct = Structure(
        lattice=lattice,
        species=[element] * 3,
        coords=dice_basis,
        coords_are_cartesian=False
    )

    struct = struct * [size[0], size[1], 1]

    return {
        "success": True,
        "lattice_type": "dice_T3",
        "a": a,
        "alpha": alpha,
        "n_atoms": len(struct),
        "size": list(size),
        "band_structure_features": [
            "Flat band touching Dirac cone at K",
            "Pseudospin-1 Dirac fermions",
            "Berry phase π (not 2π like graphene)"
        ],
        "applications": [
            "Optical lattices",
            "Photonic crystals",
            "Decorated graphene"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_honeycomb_variants(
    variant: str = "alpha_T3",
    element: str = "C",
    a: float = 2.46,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    alpha: float = 0.5
) -> Dict[str, Any]:
    """
    Generate honeycomb lattice variants.

    Various modifications of the honeycomb lattice with different
    topological and electronic properties.

    Args:
        variant: "alpha_T3", "kekulé", "semenoff", "haldane", "decorated"
        element: Element symbol
        a: Lattice parameter
        size: Supercell size
        vacuum: Vacuum thickness
        alpha: Parameter for alpha-T3 interpolation

    Returns:
        Honeycomb variant structure
    """
    from pymatgen.core import Structure, Lattice

    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a * np.sqrt(3)/2, 0])
    a3 = np.array([0, 0, vacuum])

    if variant == "alpha_T3":
        # Interpolates between honeycomb (α=0) and dice (α=1)
        lattice = Lattice([a1, a2, a3])

        coords = [
            [0.0, 0.0, 0.5],      # A
            [1/3, 1/3, 0.5],      # B
            [2/3, 2/3, 0.5]       # Hub (weight controlled by α)
        ]
        species = [element] * 3

        description = f"α-T3 lattice with α={alpha} (0=honeycomb, 1=dice)"

    elif variant == "kekule":
        # Kekulé distortion: 3x3 supercell with bond alternation
        lattice = Lattice([3*a1, 3*a2, a3])

        # Honeycomb with Kekulé pattern
        coords = []
        species = []

        for i in range(3):
            for j in range(3):
                coords.append([(i + 0.0) / 3, (j + 0.0) / 3, 0.5])
                coords.append([(i + 1/3) / 3, (j + 1/3) / 3, 0.5])
                species.extend([element, element])

        description = "Kekulé distorted honeycomb (opens gap at K)"

    elif variant == "semenoff":
        # Semenoff mass: staggered on-site potential
        lattice = Lattice([a1, a2, a3])

        coords = [
            [0.0, 0.0, 0.5],
            [1/3, 1/3, 0.5]
        ]
        # Different "elements" represent different on-site energies
        species = ["A", "B"]  # Would be same element with different potential

        description = "Semenoff mass term (trivial insulator)"

    elif variant == "haldane":
        # Haldane model: honeycomb with complex NNN hopping
        lattice = Lattice([a1, a2, a3])

        coords = [
            [0.0, 0.0, 0.5],
            [1/3, 1/3, 0.5]
        ]
        species = [element] * 2

        description = "Haldane model geometry (Chern insulator with time-reversal breaking)"

    elif variant == "decorated":
        # Decorated honeycomb: extra atoms on bonds
        lattice = Lattice([a1, a2, a3])

        # A and B sites plus 3 edge sites
        coords = [
            [0.0, 0.0, 0.5],      # A
            [1/3, 1/3, 0.5],      # B
            [1/6, 1/6, 0.5],      # Edge AB
            [1/2, 1/6, 0.5],      # Edge
            [1/6, 1/2, 0.5]       # Edge
        ]
        species = [element] * 5

        description = "Decorated honeycomb (5 atoms/cell)"

    else:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_VARIANT",
                "message": f"Unknown honeycomb variant: {variant}",
                "available": ["alpha_T3", "kekule", "semenoff", "haldane", "decorated"]
            }
        }

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    struct = struct * [size[0], size[1], 1]

    return {
        "success": True,
        "lattice_type": f"honeycomb_{variant}",
        "variant": variant,
        "description": description,
        "a": a,
        "alpha": alpha if variant == "alpha_T3" else None,
        "n_atoms": len(struct),
        "size": list(size),
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_ruby_lattice(
    element: str = "Fe",
    a: float = 5.0,
    size: Tuple[int, int] = (3, 3),
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate Ruby lattice (edge-centered Kagome).

    Ruby lattice is a Kagome lattice with additional edge sites,
    featuring multiple flat bands.

    Args:
        element: Element for sites
        a: Lattice parameter
        size: Supercell size
        vacuum: Vacuum thickness

    Returns:
        Ruby lattice structure
    """
    from pymatgen.core import Structure, Lattice

    # Ruby = Kagome + edge decoration
    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a * np.sqrt(3)/2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # Kagome basis + edge sites
    coords = [
        # Kagome sites
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.25, 0.5, 0.5],
        # Edge decorations
        [0.25, 0.0, 0.5],
        [0.625, 0.25, 0.5],
        [0.125, 0.25, 0.5]
    ]

    struct = Structure(
        lattice=lattice,
        species=[element] * 6,
        coords=coords,
        coords_are_cartesian=False
    )

    struct = struct * [size[0], size[1], 1]

    return {
        "success": True,
        "lattice_type": "ruby",
        "a": a,
        "n_atoms": len(struct),
        "size": list(size),
        "band_structure_features": [
            "Multiple flat bands",
            "Related to Kagome by decoration"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }
