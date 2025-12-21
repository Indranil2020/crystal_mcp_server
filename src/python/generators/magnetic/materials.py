"""
magnetic/materials.py - Magnetic Materials

Comprehensive magnetic material generation:
- Ferromagnets, antiferromagnets, ferrimagnets
- Magnetic orderings (FM, AFM, FiM, helical, skyrmion)
- Magnetic anisotropy
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Magnetic material database
MAGNETIC_MATERIAL_DATABASE = {
    # Elemental ferromagnets
    "Fe_bcc": {"a": 2.87, "ordering": "FM", "Tc_K": 1043, "moment_muB": 2.22,
               "anisotropy": "cubic", "K1_MJ_m3": 0.048},
    "Co_hcp": {"a": 2.51, "c": 4.07, "ordering": "FM", "Tc_K": 1388, "moment_muB": 1.72,
               "anisotropy": "uniaxial", "K1_MJ_m3": 0.53},
    "Ni_fcc": {"a": 3.52, "ordering": "FM", "Tc_K": 631, "moment_muB": 0.62,
               "anisotropy": "cubic", "K1_MJ_m3": -0.005},
    "Gd_hcp": {"a": 3.63, "c": 5.78, "ordering": "FM", "Tc_K": 293, "moment_muB": 7.63,
               "rare_earth": True},
    
    # Antiferromagnets
    "Cr": {"a": 2.88, "ordering": "AFM", "TN_K": 311, "type": "SDW"},
    "Mn_alpha": {"ordering": "AFM", "TN_K": 95, "complex": True},
    "FeO": {"a": 4.31, "ordering": "AFM", "TN_K": 198, "structure": "rocksalt"},
    "NiO": {"a": 4.17, "ordering": "AFM", "TN_K": 523, "structure": "rocksalt",
            "exchange_bias": True},
    "MnO": {"a": 4.44, "ordering": "AFM", "TN_K": 118, "structure": "rocksalt"},
    "CoO": {"a": 4.26, "ordering": "AFM", "TN_K": 291, "structure": "rocksalt"},
    
    # Ferrimagnets
    "Fe3O4": {"a": 8.40, "ordering": "FiM", "Tc_K": 858, "moment_muB": 4.0,
              "structure": "spinel", "half_metal": True},
    "Y3Fe5O12": {"a": 12.38, "ordering": "FiM", "Tc_K": 560, "structure": "garnet",
                 "name": "YIG", "damping_low": True},
    "CoFe2O4": {"a": 8.39, "ordering": "FiM", "Tc_K": 793, "structure": "spinel",
                "high_coercivity": True},
    
    # Hard magnets
    "Nd2Fe14B": {"a": 8.80, "c": 12.20, "ordering": "FM", "Tc_K": 585,
                 "moment_muB": 31.8, "BHmax_MGOe": 56, "hard_magnet": True},
    "SmCo5": {"a": 5.00, "c": 3.97, "ordering": "FM", "Tc_K": 1000,
              "moment_muB": 9.0, "BHmax_MGOe": 28, "hard_magnet": True},
    "FePt_L10": {"a": 3.86, "c": 3.71, "ordering": "FM", "Tc_K": 750,
                 "moment_muB": 3.0, "K1_MJ_m3": 6.6, "hard_magnet": True},
    
    # Spin glass and frustration
    "CuMn": {"ordering": "spin_glass", "Tf_K": 10, "dilute": True},
    "AuFe": {"ordering": "spin_glass", "Tf_K": 15, "dilute": True},
    
    # 2D magnets
    "CrI3": {"a": 6.87, "ordering": "FM", "Tc_K": 61, "2D": True,
             "monolayer_FM": True, "Ising": True},
    "CrGeTe3": {"a": 6.83, "ordering": "FM", "Tc_K": 68, "2D": True},
    "Fe3GeTe2": {"a": 3.99, "c": 16.33, "ordering": "FM", "Tc_K": 230, "2D": True,
                 "itinerant": True},
    
    # Antiferromagnetic for spintronics
    "Mn3Sn": {"a": 5.67, "c": 4.53, "ordering": "AFM", "TN_K": 420,
              "Weyl_AFM": True, "AHE": True},
    "Mn3Pt": {"ordering": "AFM", "TN_K": 475, "kagome": True},
}


# Magnetic ordering types
MAGNETIC_ORDERINGS = {
    "FM": {"type": "ferromagnetic", "parallel": True, "net_moment": True},
    "AFM_A": {"type": "antiferromagnetic", "layered": True, "planes_parallel": True},
    "AFM_C": {"type": "antiferromagnetic", "checkerboard": True},
    "AFM_G": {"type": "antiferromagnetic", "all_antiparallel": True},
    "FiM": {"type": "ferrimagnetic", "sublattices": 2, "net_moment": True},
    "canted": {"type": "canted_AFM", "weak_FM": True},
    "helical": {"type": "helical", "spiral": True, "wavevector": True},
    "conical": {"type": "conical", "helical_component": True, "FM_component": True},
    "skyrmion_lattice": {"type": "skyrmion", "topological": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_magnetic_structure(
    material: str = None,
    ordering: str = None,
    supercell: List[int] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate magnetic material with specified ordering.

    Args:
        material: Material from database
        ordering: Magnetic ordering type
        supercell: Supercell dimensions
        **kwargs: Accepts aliases (magnetic_ordering, base_structure_sg, elements, composition, a)

    Returns:
        Magnetic structure with spin information
    """
    # Handle parameter aliases
    if ordering is None:
        ordering = kwargs.get('magnetic_ordering', 'FM')
    if supercell is None:
        supercell = kwargs.get('supercell', [2, 2, 2])

    # Import at function level to avoid UnboundLocalError
    from pymatgen.core import Structure, Lattice

    # If user provides elements/composition/a, generate custom structure
    if kwargs.get('elements') or kwargs.get('base_structure_sg'):
        # Custom magnetic structure from scratch

        elements = kwargs.get('elements', ['Fe'])
        composition = kwargs.get('composition', [4])
        a = kwargs.get('a', 2.87)
        sg = kwargs.get('base_structure_sg', 229)

        # Create BCC Fe-like structure
        if sg in [229, 225]:  # BCC or FCC cubic
            lattice = Lattice.cubic(a)
            if sg == 229:  # Im-3m BCC
                species = []
                coords = []
                for i, (elem, count) in enumerate(zip(elements, composition)):
                    if count >= 2:
                        species.extend([elem] * 2)
                        coords.append([0, 0, 0])
                        coords.append([0.5, 0.5, 0.5])
                    else:
                        species.append(elem)
                        coords.append([0, 0, 0])
            else:  # Fm-3m FCC
                species = []
                coords = []
                for elem, count in zip(elements, composition):
                    for _ in range(min(count, 4)):
                        species.append(elem)
                    coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]][:len(species)]

            structure = Structure(lattice, species, coords)

            # Assign spins
            moment = 2.2
            spins = []
            for i in range(len(structure)):
                if ordering.lower() in ['fm', 'ferromagnetic']:
                    spins.append([0, 0, moment])
                elif ordering.lower() in ['afm', 'antiferromagnetic']:
                    spins.append([0, 0, moment if i % 2 == 0 else -moment])
                else:
                    spins.append([0, 0, moment])

            return {
                "success": True,
                "material": f"custom_{elements[0]}",
                "ordering": ordering,
                "space_group": sg,
                "structure": structure_to_dict(structure),
                "magnetic_moments": spins,
                "n_atoms": len(structure)
            }

    # Default material
    if material is None:
        material = "Fe_bcc"

    if material not in MAGNETIC_MATERIAL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material",
                      "available": list(MAGNETIC_MATERIAL_DATABASE.keys())}
        }
    
    info = MAGNETIC_MATERIAL_DATABASE[material]
    a = info.get("a", 4.0)
    c = info.get("c", a)
    
    # Create structure
    if "hcp" in material or info.get("structure") == "hcp":
        lattice = Lattice.hexagonal(a, c)
        base_coords = [[0, 0, 0], [1/3, 2/3, 0.5]]
    elif info.get("structure") == "rocksalt":
        lattice = Lattice.cubic(a)
        base_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    elif info.get("structure") == "spinel":
        lattice = Lattice.cubic(a)
        base_coords = [[0.125, 0.125, 0.125], [0.5, 0.5, 0.5], 
                       [0.25, 0.25, 0.25]]
    else:  # BCC or FCC
        lattice = Lattice.cubic(a)
        if "bcc" in material:
            base_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        else:
            base_coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    # Extract element
    elem = material.split("_")[0].replace("3", "").replace("2", "")
    if len(elem) > 2:
        elem = elem[:2]
    
    species = []
    coords = []
    spins = []
    
    nx, ny, nz = supercell
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for idx, base in enumerate(base_coords):
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz
                    
                    species.append(elem)
                    coords.append([x % 1, y % 1, z % 1])
                    
                    # Assign spin based on ordering
                    moment = info.get("moment_muB", 2.0)
                    
                    if ordering == "FM":
                        spin = [0, 0, moment]
                    elif ordering == "AFM_A":
                        spin = [0, 0, moment if k % 2 == 0 else -moment]
                    elif ordering == "AFM_C":
                        spin = [0, 0, moment if (i + j) % 2 == 0 else -moment]
                    elif ordering == "AFM_G":
                        spin = [0, 0, moment if (i + j + k) % 2 == 0 else -moment]
                    elif ordering == "FiM":
                        spin = [0, 0, moment if idx == 0 else -moment * 0.5]
                    elif ordering == "helical":
                        angle = 2 * np.pi * k / 8
                        spin = [moment * np.cos(angle), moment * np.sin(angle), 0]
                    else:
                        spin = [0, 0, moment]
                    
                    spins.append(spin)
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [2, 2, 2]:
        pass  # Already built with supercell
    
    return {
        "success": True,
        "material": material,
        "ordering": ordering,
        "Tc_or_TN_K": info.get("Tc_K", info.get("TN_K", 0)),
        "moment_muB_per_atom": info.get("moment_muB", 0),
        "is_hard_magnet": info.get("hard_magnet", False),
        "is_2D": info.get("2D", False),
        "n_atoms": len(structure),
        "spins": spins,
        "structure": structure_to_dict(structure)
    }


def generate_antiferromagnet(
    material: str = "NiO",
    ordering: str = "AFM_G",
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """Generate antiferromagnetic structure."""
    result = generate_magnetic_structure(material, ordering, supercell)
    if result["success"]:
        result["is_antiferromagnet"] = True
        result["has_exchange_bias"] = MAGNETIC_MATERIAL_DATABASE.get(material, {}).get("exchange_bias", False)
    return result


def get_magnetic_database() -> Dict[str, Any]:
    """Get database organized by type."""
    return {
        "success": True,
        "by_ordering": {
            "ferromagnets": [k for k, v in MAGNETIC_MATERIAL_DATABASE.items() if v["ordering"] == "FM"],
            "antiferromagnets": [k for k, v in MAGNETIC_MATERIAL_DATABASE.items() if v["ordering"] == "AFM"],
            "ferrimagnets": [k for k, v in MAGNETIC_MATERIAL_DATABASE.items() if v["ordering"] == "FiM"],
            "hard_magnets": [k for k, v in MAGNETIC_MATERIAL_DATABASE.items() if v.get("hard_magnet")],
            "2D_magnets": [k for k, v in MAGNETIC_MATERIAL_DATABASE.items() if v.get("2D")],
        }
    }


# Heusler alloy database
HEUSLER_DATABASE = {
    # Full Heusler (X2YZ)
    "Cu2MnAl": {"X": "Cu", "Y": "Mn", "Z": "Al", "a": 5.95, "type": "full", "moment_muB": 4.0},
    "Ni2MnGa": {"X": "Ni", "Y": "Mn", "Z": "Ga", "a": 5.82, "type": "full", "shape_memory": True},
    "Co2MnSi": {"X": "Co", "Y": "Mn", "Z": "Si", "a": 5.65, "type": "full", "half_metal": True},
    "Co2FeSi": {"X": "Co", "Y": "Fe", "Z": "Si", "a": 5.64, "type": "full", "half_metal": True},
    "Co2MnGe": {"X": "Co", "Y": "Mn", "Z": "Ge", "a": 5.74, "type": "full", "half_metal": True},
    "Fe2VAl": {"X": "Fe", "Y": "V", "Z": "Al", "a": 5.76, "type": "full", "thermoelectric": True},
    # Half Heusler (XYZ)
    "NiMnSb": {"X": "Ni", "Y": "Mn", "Z": "Sb", "a": 5.92, "type": "half", "half_metal": True},
    "CoMnSb": {"X": "Co", "Y": "Mn", "Z": "Sb", "a": 5.88, "type": "half"},
    "PtMnSb": {"X": "Pt", "Y": "Mn", "Z": "Sb", "a": 6.20, "type": "half"},
    "GdPtBi": {"X": "Gd", "Y": "Pt", "Z": "Bi", "a": 6.68, "type": "half", "topological": True},
    "LaPtBi": {"X": "La", "Y": "Pt", "Z": "Bi", "a": 6.85, "type": "half", "topological": True},
}


def generate_heusler(
    compound: str = "Co2MnSi",
    supercell: List[int] = None
) -> Dict[str, Any]:
    """
    Generate Heusler alloy structure.

    Heusler alloys are intermetallics with promising spintronic properties,
    including half-metallicity, high spin polarization, and shape memory effects.

    Args:
        compound: Heusler compound name (Cu2MnAl, NiMnSb, etc.)
        supercell: Supercell dimensions [nx, ny, nz]

    Returns:
        Heusler structure with magnetic properties

    Examples:
        >>> result = generate_heusler('Co2MnSi')
        >>> result["is_half_metal"]
        True
    """
    if supercell is None:
        supercell = [1, 1, 1]

    if compound not in HEUSLER_DATABASE:
        return {
            "success": False,
            "error": {
                "code": "INVALID_HEUSLER",
                "message": f"Unknown Heusler compound '{compound}'",
                "available": list(HEUSLER_DATABASE.keys())
            }
        }

    params = HEUSLER_DATABASE[compound]
    X = params["X"]
    Y = params["Y"]
    Z = params["Z"]
    a = params["a"]
    heusler_type = params["type"]

    lattice = Lattice.cubic(a)

    if heusler_type == "full":
        # Full Heusler: X2YZ in L21 structure (Fm-3m, #225)
        # X at (1/4, 1/4, 1/4) and (3/4, 3/4, 3/4)
        # Y at (0, 0, 0) and (1/2, 1/2, 1/2)
        # Z at (1/2, 0, 0)... but simpler: Y at origin, Z at face centers
        species = [Y, Z, X, X]
        coords = [
            [0, 0, 0],       # Y at origin
            [0.5, 0.5, 0.5], # Z at body center
            [0.25, 0.25, 0.25],   # X
            [0.75, 0.75, 0.75]    # X
        ]
    else:
        # Half Heusler: XYZ in C1b structure (F-43m, #216)
        # X at (1/4, 1/4, 1/4), Y at (0, 0, 0), Z at (1/2, 1/2, 1/2)
        species = [Y, Z, X]
        coords = [
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.25]
        ]

    structure = Structure(lattice, species, coords)

    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)

    return {
        "success": True,
        "compound": compound,
        "heusler_type": heusler_type,
        "X_element": X,
        "Y_element": Y,
        "Z_element": Z,
        "a_angstrom": a,
        "is_half_metal": params.get("half_metal", False),
        "is_shape_memory": params.get("shape_memory", False),
        "is_topological": params.get("topological", False),
        "moment_muB": params.get("moment_muB", 0),
        "space_group": 225 if heusler_type == "full" else 216,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
