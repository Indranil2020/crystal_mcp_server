"""
two_d/hea_2d.py - 2D High-Entropy Materials Generation

Provides generation of 2D high-entropy materials:
- High-entropy oxides (HEO)
- High-entropy carbides (HEC)
- High-entropy nitrides (HEN)
- High-entropy TMDs
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


def generate_2d_hea(
    elements: List[str],
    anion: str = "O",
    structure_type: str = "oxide",
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy material.

    Creates high-entropy alloys in 2D form for various applications
    including catalysis, energy storage, and electronics.

    Args:
        elements: List of cation elements (e.g., ["Ti", "V", "Cr", "Nb", "Ta"])
        anion: Anion element ("O", "C", "N", "S", "Se")
        structure_type: "oxide", "carbide", "nitride", "tmd", "mxene"
        composition: Molar fractions (equimolar if None)
        size: Supercell size (nx, ny)
        vacuum: Vacuum thickness in Angstroms
        ordering: "random", "sqs", "ordered"
        seed: Random seed for reproducibility

    Returns:
        2D high-entropy material structure
    """
    if structure_type == "oxide":
        return generate_2d_heo(elements, composition, size, vacuum, ordering, seed)
    elif structure_type == "carbide":
        return generate_2d_hec(elements, composition, size, vacuum, ordering, seed)
    elif structure_type == "nitride":
        return generate_2d_hen(elements, composition, size, vacuum, ordering, seed)
    elif structure_type == "tmd":
        return generate_2d_he_tmd(elements, anion, composition, size, vacuum, ordering, seed)
    elif structure_type == "mxene":
        return generate_2d_he_mxene(elements, composition, size, vacuum, ordering, seed)
    else:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_TYPE",
                "message": f"Unknown structure type: {structure_type}",
                "available": ["oxide", "carbide", "nitride", "tmd", "mxene"]
            }
        }


def generate_2d_heo(
    elements: List[str] = ["Ti", "V", "Cr", "Mn", "Fe"],
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy oxide.

    Creates monolayer or few-layer high-entropy oxide structures.

    Args:
        elements: Cation elements (typically 5+)
        composition: Molar fractions
        size: Supercell size
        vacuum: Vacuum thickness
        ordering: "random" or "sqs"
        seed: Random seed

    Returns:
        2D HEO structure
    """
    from pymatgen.core import Structure, Lattice

    if seed is not None:
        np.random.seed(seed)

    n_elements = len(elements)
    if n_elements < 4:
        return {
            "success": False,
            "error": {
                "code": "INSUFFICIENT_ELEMENTS",
                "message": "High-entropy materials require at least 4 elements"
            }
        }

    # Default to equimolar
    if composition is None:
        composition = [1.0 / n_elements] * n_elements

    # Rocksalt (111) monolayer as base
    a = 4.2  # Approximate average rocksalt parameter

    # Hexagonal surface unit cell from (111) of rocksalt
    a_hex = a / np.sqrt(2)

    a1 = np.array([a_hex, 0, 0])
    a2 = np.array([a_hex / 2, a_hex * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # Single (111) layer: cation-oxygen-cation trilayer
    base_coords = [
        [0.0, 0.0, 0.45],     # Cation bottom
        [1/3, 1/3, 0.50],     # Oxygen
        [2/3, 2/3, 0.55]      # Cation top
    ]

    # Create supercell manually
    all_coords = []
    all_species = []

    total_cation_sites = size[0] * size[1] * 2  # Two cation sites per cell

    # Assign cations based on composition
    cation_assignments = []
    for i, (elem, frac) in enumerate(zip(elements, composition)):
        n_atoms = int(round(frac * total_cation_sites))
        cation_assignments.extend([elem] * n_atoms)

    # Pad or trim to exact count
    while len(cation_assignments) < total_cation_sites:
        cation_assignments.append(elements[0])
    cation_assignments = cation_assignments[:total_cation_sites]

    if ordering == "random":
        np.random.shuffle(cation_assignments)

    cation_idx = 0
    for i in range(size[0]):
        for j in range(size[1]):
            # Bottom cation
            all_coords.append([
                (base_coords[0][0] + i) / size[0],
                (base_coords[0][1] + j) / size[1],
                base_coords[0][2]
            ])
            all_species.append(cation_assignments[cation_idx])
            cation_idx += 1

            # Oxygen
            all_coords.append([
                (base_coords[1][0] + i) / size[0],
                (base_coords[1][1] + j) / size[1],
                base_coords[1][2]
            ])
            all_species.append("O")

            # Top cation
            all_coords.append([
                (base_coords[2][0] + i) / size[0],
                (base_coords[2][1] + j) / size[1],
                base_coords[2][2]
            ])
            all_species.append(cation_assignments[min(cation_idx, len(cation_assignments)-1)])
            cation_idx = min(cation_idx + 1, len(cation_assignments) - 1)

    # Scale lattice for supercell
    super_lattice = Lattice([
        size[0] * a1,
        size[1] * a2,
        a3
    ])

    # Adjust coordinates for supercell
    super_coords = [[c[0], c[1], c[2]] for c in all_coords]

    struct = Structure(
        lattice=super_lattice,
        species=all_species,
        coords=super_coords,
        coords_are_cartesian=False
    )

    # Calculate configurational entropy
    from collections import Counter
    cation_counts = Counter([s for s in all_species if s != "O"])
    total_cations = sum(cation_counts.values())
    entropy = 0
    for count in cation_counts.values():
        p = count / total_cations
        if p > 0:
            entropy -= p * np.log(p)
    entropy_R = entropy  # In units of R

    return {
        "success": True,
        "material_type": "2D_HEO",
        "elements": elements,
        "composition": composition,
        "n_elements": n_elements,
        "ordering": ordering,
        "n_atoms": len(struct),
        "cation_distribution": dict(cation_counts),
        "configurational_entropy_R": entropy_R,
        "is_high_entropy": entropy_R >= 1.5,  # Often defined as S >= 1.5R
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_2d_hec(
    elements: List[str] = ["Ti", "Zr", "Hf", "Nb", "Ta"],
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy carbide (MXene-like).

    Args:
        elements: Metal elements
        composition: Molar fractions
        size: Supercell size
        vacuum: Vacuum thickness
        ordering: "random" or "sqs"
        seed: Random seed

    Returns:
        2D HEC structure
    """
    from pymatgen.core import Structure, Lattice

    if seed is not None:
        np.random.seed(seed)

    n_elements = len(elements)
    if composition is None:
        composition = [1.0 / n_elements] * n_elements

    # M2C MXene-like structure
    a = 3.1  # Approximate

    a1 = np.array([a, 0, 0])
    a2 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # M2C layer: M-C-M sandwich
    base_coords = [
        [0.0, 0.0, 0.45],     # Metal bottom
        [1/3, 1/3, 0.50],     # Carbon
        [2/3, 2/3, 0.55]      # Metal top
    ]

    all_coords = []
    all_species = []

    total_metal_sites = size[0] * size[1] * 2

    metal_assignments = []
    for elem, frac in zip(elements, composition):
        n_atoms = int(round(frac * total_metal_sites))
        metal_assignments.extend([elem] * n_atoms)

    while len(metal_assignments) < total_metal_sites:
        metal_assignments.append(elements[0])
    metal_assignments = metal_assignments[:total_metal_sites]

    if ordering == "random":
        np.random.shuffle(metal_assignments)

    metal_idx = 0
    for i in range(size[0]):
        for j in range(size[1]):
            all_coords.append([
                (base_coords[0][0] + i) / size[0],
                (base_coords[0][1] + j) / size[1],
                base_coords[0][2]
            ])
            all_species.append(metal_assignments[metal_idx])
            metal_idx += 1

            all_coords.append([
                (base_coords[1][0] + i) / size[0],
                (base_coords[1][1] + j) / size[1],
                base_coords[1][2]
            ])
            all_species.append("C")

            all_coords.append([
                (base_coords[2][0] + i) / size[0],
                (base_coords[2][1] + j) / size[1],
                base_coords[2][2]
            ])
            all_species.append(metal_assignments[min(metal_idx, len(metal_assignments)-1)])
            metal_idx = min(metal_idx + 1, len(metal_assignments) - 1)

    super_lattice = Lattice([size[0] * a1, size[1] * a2, a3])

    struct = Structure(
        lattice=super_lattice,
        species=all_species,
        coords=all_coords,
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "material_type": "2D_HEC",
        "elements": elements,
        "composition": composition,
        "n_elements": n_elements,
        "structure_type": "M2C_MXene",
        "n_atoms": len(struct),
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_2d_hen(
    elements: List[str] = ["Ti", "V", "Cr", "Nb", "Ta"],
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy nitride.

    Args:
        elements: Metal elements
        composition: Molar fractions
        size: Supercell size
        vacuum: Vacuum thickness
        ordering: "random" or "sqs"
        seed: Random seed

    Returns:
        2D HEN structure
    """
    # Similar to HEC but with N instead of C
    result = generate_2d_hec(elements, composition, size, vacuum, ordering, seed)

    if result["success"]:
        # Replace C with N in the structure
        struct_dict = result["pymatgen_structure"]
        for site in struct_dict["sites"]:
            if site["species"][0]["element"] == "C":
                site["species"][0]["element"] = "N"

        result["material_type"] = "2D_HEN"
        result["structure_type"] = "M2N_nitride"

        # Update atoms list
        for atom in result["structure"]["atoms"]:
            if atom["element"] == "C":
                atom["element"] = "N"

    return result


def generate_2d_he_tmd(
    metals: List[str] = ["Mo", "W", "Nb", "Ta", "V"],
    chalcogen: str = "S",
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy TMD.

    Args:
        metals: Transition metal elements
        chalcogen: "S", "Se", or "Te"
        composition: Molar fractions
        size: Supercell size
        vacuum: Vacuum thickness
        ordering: "random" or "sqs"
        seed: Random seed

    Returns:
        2D HE-TMD structure
    """
    from pymatgen.core import Structure, Lattice

    if seed is not None:
        np.random.seed(seed)

    n_elements = len(metals)
    if composition is None:
        composition = [1.0 / n_elements] * n_elements

    # 1H-MX2 structure
    a = 3.2

    a1 = np.array([a, 0, 0])
    a2 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # 1H trilayer: X-M-X
    base_coords = [
        [1/3, 2/3, 0.45],     # Bottom chalcogen
        [0.0, 0.0, 0.50],     # Metal
        [1/3, 2/3, 0.55]      # Top chalcogen
    ]

    all_coords = []
    all_species = []

    total_metal_sites = size[0] * size[1]

    metal_assignments = []
    for elem, frac in zip(metals, composition):
        n_atoms = int(round(frac * total_metal_sites))
        metal_assignments.extend([elem] * n_atoms)

    while len(metal_assignments) < total_metal_sites:
        metal_assignments.append(metals[0])
    metal_assignments = metal_assignments[:total_metal_sites]

    if ordering == "random":
        np.random.shuffle(metal_assignments)

    metal_idx = 0
    for i in range(size[0]):
        for j in range(size[1]):
            # Bottom chalcogen
            all_coords.append([
                (base_coords[0][0] + i) / size[0],
                (base_coords[0][1] + j) / size[1],
                base_coords[0][2]
            ])
            all_species.append(chalcogen)

            # Metal
            all_coords.append([
                (base_coords[1][0] + i) / size[0],
                (base_coords[1][1] + j) / size[1],
                base_coords[1][2]
            ])
            all_species.append(metal_assignments[metal_idx])
            metal_idx += 1

            # Top chalcogen
            all_coords.append([
                (base_coords[2][0] + i) / size[0],
                (base_coords[2][1] + j) / size[1],
                base_coords[2][2]
            ])
            all_species.append(chalcogen)

    super_lattice = Lattice([size[0] * a1, size[1] * a2, a3])

    struct = Structure(
        lattice=super_lattice,
        species=all_species,
        coords=all_coords,
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "material_type": "2D_HE_TMD",
        "metals": metals,
        "chalcogen": chalcogen,
        "composition": composition,
        "phase": "1H",
        "n_atoms": len(struct),
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_2d_he_mxene(
    metals: List[str] = ["Ti", "V", "Nb", "Ta", "Mo"],
    composition: Optional[List[float]] = None,
    size: Tuple[int, int] = (4, 4),
    vacuum: float = 15.0,
    termination: str = "O",
    ordering: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate 2D high-entropy MXene.

    Args:
        metals: Transition metal elements
        composition: Molar fractions
        size: Supercell size
        vacuum: Vacuum thickness
        termination: Surface termination ("O", "OH", "F", "bare")
        ordering: "random" or "sqs"
        seed: Random seed

    Returns:
        2D HE-MXene structure
    """
    from pymatgen.core import Structure, Lattice

    if seed is not None:
        np.random.seed(seed)

    n_elements = len(metals)
    if composition is None:
        composition = [1.0 / n_elements] * n_elements

    # M2CTx structure
    a = 3.1

    a1 = np.array([a, 0, 0])
    a2 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, vacuum])

    lattice = Lattice([a1, a2, a3])

    # Full M2CTx layer with terminations
    if termination == "bare":
        base_coords = [
            [0.0, 0.0, 0.45],     # Metal bottom
            [1/3, 1/3, 0.50],     # Carbon
            [2/3, 2/3, 0.55]      # Metal top
        ]
        base_species = ["M", "C", "M"]
    else:
        base_coords = [
            [1/3, 2/3, 0.40],     # Bottom termination
            [0.0, 0.0, 0.45],     # Metal bottom
            [1/3, 1/3, 0.50],     # Carbon
            [2/3, 2/3, 0.55],     # Metal top
            [1/3, 2/3, 0.60]      # Top termination
        ]
        base_species = [termination, "M", "C", "M", termination]

    all_coords = []
    all_species = []

    metal_count = sum(1 for s in base_species if s == "M")
    total_metal_sites = size[0] * size[1] * metal_count

    metal_assignments = []
    for elem, frac in zip(metals, composition):
        n_atoms = int(round(frac * total_metal_sites))
        metal_assignments.extend([elem] * n_atoms)

    while len(metal_assignments) < total_metal_sites:
        metal_assignments.append(metals[0])
    metal_assignments = metal_assignments[:total_metal_sites]

    if ordering == "random":
        np.random.shuffle(metal_assignments)

    metal_idx = 0
    for i in range(size[0]):
        for j in range(size[1]):
            for k, (coord, sp) in enumerate(zip(base_coords, base_species)):
                all_coords.append([
                    (coord[0] + i) / size[0],
                    (coord[1] + j) / size[1],
                    coord[2]
                ])
                if sp == "M":
                    all_species.append(metal_assignments[metal_idx])
                    metal_idx += 1
                else:
                    all_species.append(sp)

    super_lattice = Lattice([size[0] * a1, size[1] * a2, a3])

    struct = Structure(
        lattice=super_lattice,
        species=all_species,
        coords=all_coords,
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "material_type": "2D_HE_MXene",
        "metals": metals,
        "composition": composition,
        "structure_type": "M2CTx",
        "termination": termination,
        "n_atoms": len(struct),
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }
