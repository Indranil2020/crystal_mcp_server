"""
battery/anodes.py - Battery Anode Structure Generation

Provides anode structure generation for:
- Conversion anodes (Si, Ge, Sn, P)
- Alloy anodes and lithiation stages
- Metal anodes with SEI and dendrites
- Sulfur cathodes (S8 in hosts)
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


def generate_silicon_anode(
    phase: str = "crystalline",
    lithiation: float = 0.0,
    supercell: List[int] = [2, 2, 2],
    include_surface: bool = False,
    surface_miller: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate silicon anode structure at various lithiation stages.

    Silicon undergoes massive volume expansion during lithiation,
    forming various LixSi phases.

    Args:
        phase: "crystalline", "amorphous", "nanoparticle"
        lithiation: Li/Si ratio (0 = pure Si, 3.75 = Li15Si4, 4.4 = Li22Si5)
        supercell: Supercell dimensions
        include_surface: Generate surface slab
        surface_miller: Surface orientation for slab

    Returns:
        Silicon anode structure
    """
    from pymatgen.core import Structure, Lattice

    # Silicon phases during lithiation
    si_phases = {
        0.0: {  # Pure Si
            "name": "diamond_Si",
            "lattice": Lattice.cubic(5.431),
            "species": ["Si"] * 8,
            "coords": [
                [0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
                [0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
                [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
            ],
            "spacegroup": "Fd-3m",
            "volume_expansion": 1.0
        },
        1.0: {  # LiSi
            "name": "LiSi_I41/a",
            "a": 9.36, "c": 5.75,
            "species": ["Li"] * 8 + ["Si"] * 8,
            "volume_expansion": 2.0
        },
        2.33: {  # Li7Si3
            "name": "Li7Si3_R-3m",
            "volume_expansion": 2.5
        },
        3.25: {  # Li13Si4
            "name": "Li13Si4_Pbam",
            "volume_expansion": 3.2
        },
        3.75: {  # Li15Si4
            "name": "Li15Si4_I-43d",
            "a": 10.63,
            "species": ["Li"] * 60 + ["Si"] * 16,
            "volume_expansion": 3.5
        },
        4.4: {  # Li22Si5 (max lithiation)
            "name": "Li22Si5_F-43m",
            "a": 18.75,
            "volume_expansion": 4.0
        }
    }

    # Find closest phase
    closest_x = min(si_phases.keys(), key=lambda x: abs(x - lithiation))
    params = si_phases[closest_x]

    if closest_x == 0.0:
        # Pure silicon
        struct = Structure(
            lattice=params["lattice"],
            species=params["species"],
            coords=params["coords"],
            coords_are_cartesian=False
        )
    elif closest_x == 3.75:
        # Li15Si4 - important fully lithiated phase
        lattice = Lattice.cubic(params["a"])
        # Simplified structure (actual has 76 atoms per cell)
        struct = Structure(
            lattice=lattice,
            species=["Li"] * 15 + ["Si"] * 4,
            coords=[
                # Li positions (simplified)
                *[[i*0.1, j*0.1, k*0.1] for i in range(3) for j in range(3) for k in range(2)][:15],
                # Si positions
                [0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
                [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
            ],
            coords_are_cartesian=False
        )
    else:
        # Interpolated phase - simplified model
        # Use linear interpolation of Si structure with Li added
        base_struct = Structure(
            lattice=Lattice.cubic(5.431 * (params["volume_expansion"] ** (1/3))),
            species=["Si"] * 8,
            coords=[
                [0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
                [0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
                [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
            ],
            coords_are_cartesian=False
        )

        # Add Li at interstitial sites
        n_li = int(lithiation * 8)
        li_sites = [
            [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5],
            [0.5, 0.5, 0.5], [0.125, 0.125, 0.125], [0.375, 0.375, 0.375],
            [0.625, 0.625, 0.625], [0.875, 0.875, 0.875]
        ]

        species = ["Si"] * 8 + ["Li"] * min(n_li, len(li_sites))
        coords = list(base_struct.frac_coords) + li_sites[:min(n_li, len(li_sites))]

        struct = Structure(
            lattice=base_struct.lattice,
            species=species,
            coords=coords,
            coords_are_cartesian=False
        )

    # Apply supercell
    struct = struct * supercell

    # Handle amorphous or nanoparticle
    if phase == "amorphous":
        # Add random displacements
        new_coords = []
        for site in struct.sites:
            displacement = np.random.normal(0, 0.1, 3)
            new_coords.append(site.coords + displacement)
        struct = Structure(
            lattice=struct.lattice,
            species=[str(s.specie) for s in struct.sites],
            coords=new_coords,
            coords_are_cartesian=True
        )

    return {
        "success": True,
        "material": "Si",
        "phase": phase,
        "lithiation_x": lithiation,
        "closest_phase": params["name"],
        "volume_expansion": params["volume_expansion"],
        "Li_Si_ratio": lithiation,
        "n_atoms": len(struct),
        "composition": struct.composition.reduced_formula,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_germanium_anode(
    lithiation: float = 0.0,
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate germanium anode structure.

    Ge has better rate capability than Si due to faster Li diffusion.

    Args:
        lithiation: Li/Ge ratio (0 to 4.4)
        supercell: Supercell dimensions

    Returns:
        Ge anode structure
    """
    from pymatgen.core import Structure, Lattice

    # Ge has diamond structure like Si
    ge_a = 5.658  # Lattice parameter

    # Volume expansion with lithiation
    volume_expansion = 1 + lithiation * 0.75  # ~3.7x at full lithiation

    lattice = Lattice.cubic(ge_a * (volume_expansion ** (1/3)))

    # Diamond structure base
    ge_coords = [
        [0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
        [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
    ]

    # Li interstitial sites
    li_sites = [
        [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5], [0.125, 0.125, 0.125], [0.375, 0.375, 0.375]
    ]

    n_li = int(lithiation * 8)
    species = ["Ge"] * 8 + ["Li"] * min(n_li, len(li_sites) * 4)

    # Add Li sites (repeat if needed)
    all_li = li_sites * 4
    coords = ge_coords + all_li[:min(n_li, len(all_li))]

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    struct = struct * supercell

    return {
        "success": True,
        "material": "Ge",
        "lithiation_x": lithiation,
        "volume_expansion": volume_expansion,
        "n_atoms": len(struct),
        "composition": struct.composition.reduced_formula,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_tin_anode(
    phase: str = "beta",
    lithiation: float = 0.0,
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate tin anode structure.

    Sn forms various alloys with Li up to Li22Sn5.

    Args:
        phase: "alpha" (gray), "beta" (white/metallic)
        lithiation: Li/Sn ratio (0 to 4.4)
        supercell: Supercell dimensions

    Returns:
        Sn anode structure
    """
    from pymatgen.core import Structure, Lattice

    if phase == "beta":
        # Beta-Sn (BCT, white tin)
        a = 5.832
        c = 3.182
        lattice = Lattice.tetragonal(a, c)
        coords = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ]
        species = ["Sn"] * 2
    else:
        # Alpha-Sn (diamond, gray tin)
        a = 6.489
        lattice = Lattice.cubic(a)
        coords = [
            [0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
            [0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
        ]
        species = ["Sn"] * 8

    # Volume expansion
    volume_expansion = 1 + lithiation * 0.68

    struct = Structure(
        lattice=Lattice(lattice.matrix * (volume_expansion ** (1/3))),
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    # Add Li if lithiated
    if lithiation > 0:
        n_li = int(lithiation * len([s for s in species if s == "Sn"]))
        li_sites = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]] * (n_li // 2 + 1)

        new_species = list(species) + ["Li"] * n_li
        new_coords = list(coords) + li_sites[:n_li]

        struct = Structure(
            lattice=struct.lattice,
            species=new_species,
            coords=new_coords,
            coords_are_cartesian=False
        )

    struct = struct * supercell

    return {
        "success": True,
        "material": "Sn",
        "phase": phase,
        "lithiation_x": lithiation,
        "volume_expansion": volume_expansion,
        "n_atoms": len(struct),
        "composition": struct.composition.reduced_formula,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_sulfur_cathode(
    host: str = "carbon",
    sulfur_loading: str = "S8",
    host_structure: str = "microporous",
    supercell: List[int] = [2, 2, 1]
) -> Dict[str, Any]:
    """
    Generate sulfur cathode in host material.

    For Li-S batteries, sulfur is encapsulated in conductive hosts.

    Args:
        host: "carbon", "MXene", "MOF", "graphene"
        sulfur_loading: "S8" (orthorhombic), "Li2S", "polysulfides"
        host_structure: "microporous", "mesoporous", "layered"
        supercell: Supercell dimensions

    Returns:
        Sulfur cathode structure
    """
    from pymatgen.core import Structure, Lattice, Molecule

    # S8 crown molecule coordinates
    s8_coords = []
    r = 2.05  # S-S bond length
    theta = 108  # Dihedral angle
    for i in range(8):
        angle = i * 2 * np.pi / 8
        z = 0.5 * np.sin(i * np.pi / 4)  # Crown shape
        s8_coords.append([
            r * np.cos(angle),
            r * np.sin(angle),
            z
        ])

    if host == "carbon":
        if host_structure == "microporous":
            # Microporous carbon cage
            a = 12.0
            lattice = Lattice.cubic(a)

            # Carbon cage atoms (simplified)
            c_coords = []
            for i in range(4):
                for j in range(4):
                    for k in range(4):
                        if (i + j + k) % 2 == 0:
                            c_coords.append([i/4, j/4, k/4])

            # Place S8 in center
            s8_centered = [[c[0]/a + 0.5, c[1]/a + 0.5, c[2]/a + 0.5] for c in s8_coords]

            species = ["C"] * len(c_coords) + ["S"] * 8
            coords = c_coords + s8_centered

        elif host_structure == "layered":
            # Graphene layers with S8 between
            a = 15.0
            c = 10.0
            lattice = Lattice.tetragonal(a, c)

            # Two graphene layers
            c_coords = []
            for i in range(5):
                for j in range(5):
                    c_coords.append([i/5, j/5, 0.1])  # Bottom layer
                    c_coords.append([i/5 + 0.1, j/5 + 0.1, 0.9])  # Top layer

            # S8 between layers
            s8_centered = [[c[0]/a + 0.5, c[1]/a + 0.5, c[2]/c + 0.5] for c in s8_coords]

            species = ["C"] * len(c_coords) + ["S"] * 8
            coords = c_coords + s8_centered

    elif host == "MXene":
        # Ti3C2 MXene layers
        a = 12.0
        c = 15.0
        lattice = Lattice.tetragonal(a, c)

        # Simplified MXene structure
        ti_coords = [[0.0, 0.0, 0.1], [0.5, 0.0, 0.1], [0.0, 0.5, 0.1], [0.5, 0.5, 0.1]]
        c_coords = [[0.25, 0.25, 0.15], [0.75, 0.75, 0.15]]

        # S8 in gallery
        s8_centered = [[c[0]/a + 0.5, c[1]/a + 0.5, 0.5] for c in s8_coords]

        species = ["Ti"] * 4 + ["C"] * 2 + ["S"] * 8
        coords = ti_coords + c_coords + s8_centered

    else:
        # Default graphene host
        a = 15.0
        c = 10.0
        lattice = Lattice.tetragonal(a, c)
        s8_centered = [[c[0]/a + 0.5, c[1]/a + 0.5, 0.5] for c in s8_coords]
        species = ["S"] * 8
        coords = s8_centered

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    struct = struct * supercell

    return {
        "success": True,
        "material": "S cathode",
        "host": host,
        "host_structure": host_structure,
        "sulfur_loading": sulfur_loading,
        "n_atoms": len(struct),
        "composition": struct.composition.reduced_formula,
        "sulfur_content_wt_percent": 8 * 32.06 / (sum(s.specie.atomic_mass for s in struct.sites)) * 100,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_metal_anode(
    metal: str = "Li",
    surface_miller: List[int] = [1, 0, 0],
    n_layers: int = 6,
    include_sei: bool = False,
    sei_composition: str = "LiF",
    include_dendrite: bool = False,
    dendrite_height_A: float = 10.0,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate metal anode with optional SEI and dendrite.

    Creates Li/Na/K metal surfaces with solid electrolyte interphase
    and dendrite morphology for safety studies.

    Args:
        metal: "Li", "Na", "K"
        surface_miller: Surface orientation
        n_layers: Number of metal layers
        include_sei: Add SEI layer
        sei_composition: SEI material ("LiF", "Li2CO3", "Li2O")
        include_dendrite: Add dendrite on surface
        dendrite_height_A: Height of dendrite
        vacuum: Vacuum thickness

    Returns:
        Metal anode structure
    """
    from pymatgen.core import Structure, Lattice
    from pymatgen.core.surface import SlabGenerator

    # Metal lattice parameters (BCC)
    metal_params = {
        "Li": {"a": 3.51, "structure": "bcc"},
        "Na": {"a": 4.29, "structure": "bcc"},
        "K": {"a": 5.33, "structure": "bcc"}
    }

    if metal not in metal_params:
        return {
            "success": False,
            "error": {
                "code": "INVALID_METAL",
                "message": f"Unknown metal '{metal}'",
                "available": list(metal_params.keys())
            }
        }

    params = metal_params[metal]

    # Create bulk metal
    lattice = Lattice.cubic(params["a"])
    bulk = Structure(
        lattice=lattice,
        species=[metal, metal],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
        coords_are_cartesian=False
    )

    # Generate surface slab
    slabgen = SlabGenerator(
        bulk,
        miller_index=tuple(surface_miller),
        min_slab_size=n_layers * params["a"] / 2,
        min_vacuum_size=vacuum,
        center_slab=True
    )

    slabs = slabgen.get_slabs()
    if not slabs:
        return {
            "success": False,
            "error": {"message": "Could not generate slab"}
        }

    slab = slabs[0]

    # Add SEI layer if requested
    sei_info = None
    if include_sei:
        # SEI layer on top of metal
        sei_thickness_A = 3.0

        sei_materials = {
            "LiF": {"species": ["Li", "F"], "coords": [[0, 0], [0.5, 0.5]], "a": 4.03},
            "Li2CO3": {"species": ["Li", "Li", "C", "O", "O", "O"], "a": 5.0},
            "Li2O": {"species": ["Li", "Li", "O"], "coords": [[0.25, 0.25], [0.75, 0.75], [0, 0]], "a": 4.62}
        }

        sei_params = sei_materials.get(sei_composition, sei_materials["LiF"])

        # Add SEI atoms on top of slab
        max_z = max(s.coords[2] for s in slab.sites)

        new_species = [str(s.specie) for s in slab.sites]
        new_coords = [s.coords.tolist() for s in slab.sites]

        # Add SEI layer
        for i in range(2):
            for j in range(2):
                for k, sp in enumerate(sei_params["species"]):
                    z = max_z + sei_thickness_A * (k + 1) / len(sei_params["species"])
                    new_species.append(sp)
                    new_coords.append([
                        i * slab.lattice.a / 2 + 1,
                        j * slab.lattice.b / 2 + 1,
                        z
                    ])

        slab = Structure(
            lattice=slab.lattice,
            species=new_species,
            coords=new_coords,
            coords_are_cartesian=True
        )

        sei_info = {
            "composition": sei_composition,
            "thickness_A": sei_thickness_A
        }

    # Add dendrite if requested
    dendrite_info = None
    if include_dendrite:
        max_z = max(s.coords[2] for s in slab.sites)

        new_species = [str(s.specie) for s in slab.sites]
        new_coords = [s.coords.tolist() for s in slab.sites]

        # Simple needle dendrite
        n_dendrite_atoms = int(dendrite_height_A / (params["a"] / 2))
        dendrite_x = slab.lattice.a / 2
        dendrite_y = slab.lattice.b / 2

        for i in range(n_dendrite_atoms):
            z = max_z + (i + 1) * params["a"] / 2
            new_species.append(metal)
            new_coords.append([dendrite_x, dendrite_y, z])

        slab = Structure(
            lattice=slab.lattice,
            species=new_species,
            coords=new_coords,
            coords_are_cartesian=True
        )

        dendrite_info = {
            "height_A": dendrite_height_A,
            "n_atoms": n_dendrite_atoms,
            "morphology": "needle"
        }

    return {
        "success": True,
        "metal": metal,
        "surface_miller": list(surface_miller),
        "n_layers": n_layers,
        "n_atoms": len(slab),
        "has_sei": include_sei,
        "sei_info": sei_info,
        "has_dendrite": include_dendrite,
        "dendrite_info": dendrite_info,
        "structure": {
            "lattice": {"matrix": slab.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in slab.sites
            ]
        },
        "pymatgen_structure": slab.as_dict()
    }


def generate_phosphorus_anode(
    allotrope: str = "black",
    lithiation: float = 0.0,
    supercell: List[int] = [2, 2, 1]
) -> Dict[str, Any]:
    """
    Generate phosphorus anode structure.

    Black phosphorus has layered structure suitable for Li intercalation.

    Args:
        allotrope: "black" (layered), "red" (amorphous)
        lithiation: Li/P ratio (0 to 3 for Li3P)
        supercell: Supercell dimensions

    Returns:
        Phosphorus anode structure
    """
    from pymatgen.core import Structure, Lattice

    if allotrope == "black":
        # Black phosphorus (orthorhombic, layered)
        a = 3.31
        b = 10.48
        c = 4.38

        lattice = Lattice.orthorhombic(a, b, c)

        # Phosphorus positions in puckered layer
        p_coords = [
            [0.0, 0.1, 0.08],
            [0.0, 0.4, 0.92],
            [0.5, 0.6, 0.08],
            [0.5, 0.9, 0.92]
        ]

        species = ["P"] * 4

        if lithiation > 0:
            # Add Li between layers
            n_li = int(lithiation * 4)
            li_sites = [[0, 0.25, 0.5], [0.5, 0.75, 0.5]] * (n_li // 2 + 1)
            species += ["Li"] * n_li
            p_coords += li_sites[:n_li]

        struct = Structure(
            lattice=lattice,
            species=species,
            coords=p_coords,
            coords_are_cartesian=False
        )

    else:
        # Red phosphorus (simplified as amorphous P4 clusters)
        a = 10.0
        lattice = Lattice.cubic(a)

        # P4 tetrahedron
        p_coords = [
            [0.5, 0.5, 0.6],
            [0.6, 0.4, 0.4],
            [0.4, 0.6, 0.4],
            [0.4, 0.4, 0.6]
        ]

        species = ["P"] * 4

        struct = Structure(
            lattice=lattice,
            species=species,
            coords=p_coords,
            coords_are_cartesian=False
        )

    struct = struct * supercell

    return {
        "success": True,
        "material": "P",
        "allotrope": allotrope,
        "lithiation_x": lithiation,
        "n_atoms": len(struct),
        "composition": struct.composition.reduced_formula,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }
