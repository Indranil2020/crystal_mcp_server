"""
surface/slabs.py - Surface Slab Generation

Comprehensive slab generation per structure_catalogue.md Category 4:
(i) Miller-index surfaces (hkl) up to ≥ 15 Å thick
(vii) Roughened / stepped surfaces (kinks, ad-terraces, vicinal surfaces)
"""

from typing import Dict, Any, List, Optional, Union
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator


# Surface energy database (J/m^2)
SURFACE_ENERGY_DATABASE = {
    "Au": {"(111)": 1.28, "(100)": 1.63, "(110)": 1.70},
    "Ag": {"(111)": 1.17, "(100)": 1.20, "(110)": 1.24},
    "Cu": {"(111)": 1.79, "(100)": 2.17, "(110)": 2.24},
    "Pt": {"(111)": 2.30, "(100)": 2.73, "(110)": 2.82},
    "Pd": {"(111)": 2.00, "(100)": 2.33, "(110)": 2.45},
    "Ni": {"(111)": 2.38, "(100)": 2.43, "(110)": 2.37},
    "Fe_bcc": {"(110)": 2.41, "(100)": 2.22, "(111)": 2.69},
    "Al": {"(111)": 1.14, "(100)": 1.35, "(110)": 1.35},
    "Si": {"(111)": 1.23, "(100)": 2.13, "(110)": 1.51},
    "Ge": {"(111)": 1.01, "(100)": 1.71, "(110)": 1.23},
    "TiO2_rutile": {"(110)": 0.31, "(100)": 0.53, "(001)": 0.90},
    "TiO2_anatase": {"(101)": 0.44, "(001)": 0.90, "(100)": 0.53},
    "MgO": {"(100)": 1.04, "(110)": 1.35, "(111)": 1.87},
    "Al2O3": {"(0001)": 1.59, "(10-12)": 2.03, "(11-20)": 2.17},
}


# Common slab configurations
SLAB_CONFIGURATIONS = {
    "Au_111": {"element": "Au", "structure": "fcc", "a": 4.08, "miller": [1, 1, 1]},
    "Pt_111": {"element": "Pt", "structure": "fcc", "a": 3.92, "miller": [1, 1, 1]},
    "Cu_111": {"element": "Cu", "structure": "fcc", "a": 3.61, "miller": [1, 1, 1]},
    "Ag_111": {"element": "Ag", "structure": "fcc", "a": 4.09, "miller": [1, 1, 1]},
    "Pd_111": {"element": "Pd", "structure": "fcc", "a": 3.89, "miller": [1, 1, 1]},
    "Ni_111": {"element": "Ni", "structure": "fcc", "a": 3.52, "miller": [1, 1, 1]},
    "Fe_110": {"element": "Fe", "structure": "bcc", "a": 2.87, "miller": [1, 1, 0]},
    "Si_111": {"element": "Si", "structure": "diamond", "a": 5.43, "miller": [1, 1, 1]},
    "Ge_111": {"element": "Ge", "structure": "diamond", "a": 5.66, "miller": [1, 1, 1]},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_slab(
    bulk_structure: Union[Structure, str] = None,
    miller_index: List[int] = None,
    min_slab_thickness: float = 10.0,
    min_vacuum: float = 15.0,
    in_unit_planes: bool = False,
    center_slab: bool = True,
    symmetric: bool = True,
    primitive: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate slab from bulk structure.

    Args:
        bulk_structure: Bulk structure or preset name
        miller_index: Miller indices [h, k, l]
        min_slab_thickness: Minimum slab thickness in Å
        min_vacuum: Vacuum thickness in Å
        in_unit_planes: Interpret thickness in unit planes
        center_slab: Center slab in cell
        symmetric: Make slab symmetric
        primitive: Use primitive cell
        **kwargs: Accepts aliases (structure, miller_indices, min_vacuum_thickness)

    Returns:
        Slab structure
    """
    # Handle parameter aliases
    if bulk_structure is None:
        bulk_structure = kwargs.get('structure')
    if miller_index is None:
        miller_index = kwargs.get('miller_indices', [1, 1, 1])
    if 'min_vacuum_thickness' in kwargs:
        min_vacuum = kwargs['min_vacuum_thickness']

    if bulk_structure is None:
        return {
            "success": False,
            "error": {"code": "MISSING_STRUCTURE",
                      "message": "Must provide 'bulk_structure' or 'structure'"}
        }

    # Handle preset configurations
    if isinstance(bulk_structure, str):
        if bulk_structure in SLAB_CONFIGURATIONS:
            config = SLAB_CONFIGURATIONS[bulk_structure]
            a = config["a"]
            element = config["element"]
            structure_type = config["structure"]
            miller_index = config["miller"]
            
            if structure_type == "fcc":
                lattice = Lattice.cubic(a)
                species = [element] * 4
                coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
            elif structure_type == "bcc":
                lattice = Lattice.cubic(a)
                species = [element] * 2
                coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
            elif structure_type == "diamond":
                lattice = Lattice.cubic(a)
                species = [element] * 8
                coords = [
                    [0, 0, 0], [0.25, 0.25, 0.25],
                    [0.5, 0.5, 0], [0.75, 0.75, 0.25],
                    [0.5, 0, 0.5], [0.75, 0.25, 0.75],
                    [0, 0.5, 0.5], [0.25, 0.75, 0.75]
                ]
            else:
                lattice = Lattice.cubic(a)
                species = [element]
                coords = [[0, 0, 0]]
            
            bulk_structure = Structure(lattice, species, coords)
        else:
            return {
                "success": False,
                "error": {"code": "INVALID_PRESET", "message": f"Unknown preset '{bulk_structure}'",
                          "available": list(SLAB_CONFIGURATIONS.keys())}
            }
    
    # Generate slab using SlabGenerator
    slabgen = SlabGenerator(
        bulk_structure,
        miller_index,
        min_slab_size=min_slab_thickness,
        min_vacuum_size=min_vacuum,
        in_unit_planes=in_unit_planes,
        center_slab=center_slab,
        primitive=primitive
    )
    
    slabs = slabgen.get_slabs(symmetrize=symmetric)
    
    if not slabs:
        return {
            "success": False,
            "error": {"code": "NO_SLABS", "message": "Could not generate slab for given parameters"}
        }
    
    # Take first slab
    slab = slabs[0]
    
    # Calculate slab properties
    slab_thickness = max(s.coords[2] for s in slab) - min(s.coords[2] for s in slab)
    
    return {
        "success": True,
        "miller_index": miller_index,
        "slab_thickness_A": round(slab_thickness, 2),
        "vacuum_thickness_A": min_vacuum,
        "is_symmetric": slab.is_symmetric(),
        "is_polar": slab.is_polar(),
        "n_atoms": len(slab),
        "surface_area_A2": round(slab.surface_area, 2),
        "structure": structure_to_dict(slab)
    }


def generate_vicinal_surface(
    base_surface: str = "Au_111",
    miscut_angle_deg: float = 5.0,
    miscut_direction: str = "[11-2]",
    n_terraces: int = 3
) -> Dict[str, Any]:
    """
    Generate vicinal (stepped) surface.
    
    Args:
        base_surface: Base low-index surface
        miscut_angle_deg: Miscut angle in degrees
        miscut_direction: Direction of miscut
        n_terraces: Number of terraces
    
    Returns:
        Vicinal surface structure
    """
    if base_surface not in SLAB_CONFIGURATIONS:
        return {
            "success": False,
            "error": {"code": "INVALID_SURFACE", "message": f"Unknown surface",
                      "available": list(SLAB_CONFIGURATIONS.keys())}
        }
    
    config = SLAB_CONFIGURATIONS[base_surface]
    a = config["a"]
    element = config["element"]
    base_miller = config["miller"]
    
    # Calculate terrace width from miscut angle
    d_spacing = a / np.sqrt(sum(m**2 for m in base_miller))
    terrace_width = d_spacing / np.tan(np.radians(miscut_angle_deg))
    
    # Build stepped surface
    terrace_atoms = 10  # Atoms per terrace row
    step_height = d_spacing
    
    species = []
    coords = []
    
    box_a = terrace_width * n_terraces + 5
    box_b = a * 3
    box_c = step_height * n_terraces + 30  # With vacuum
    
    lattice = Lattice.orthorhombic(box_a, box_b, box_c)
    
    for terrace in range(n_terraces):
        z_base = terrace * step_height / box_c + 0.3
        x_offset = terrace * terrace_width / box_a
        
        for layer in range(3):  # 3 atomic layers per terrace
            z = z_base + layer * (d_spacing / 3) / box_c
            
            for i in range(terrace_atoms):
                x = x_offset + (i + 0.5 * (layer % 2)) * (terrace_width / terrace_atoms) / box_a
                for j in range(3):
                    y = (j + 0.5 * (layer % 2)) / 3
                    
                    if x < 1 and y < 1:
                        species.append(element)
                        coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    # Estimate step density
    step_density = n_terraces / (terrace_width * n_terraces) * 1e8  # steps/cm
    
    return {
        "success": True,
        "base_surface": base_surface,
        "miscut_angle_deg": miscut_angle_deg,
        "miscut_direction": miscut_direction,
        "terrace_width_A": round(terrace_width, 2),
        "step_height_A": round(step_height, 2),
        "step_density_per_cm": round(step_density, 0),
        "n_terraces": n_terraces,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_kinked_surface(
    base_surface: str = "Cu_111",
    kink_density: float = 0.2,
    kink_type: str = "monatomic"
) -> Dict[str, Any]:
    """
    Generate surface with kink sites.
    
    Args:
        base_surface: Base surface
        kink_density: Fraction of step edge with kinks
        kink_type: 'monatomic' or 'diatomic'
    
    Returns:
        Kinked surface structure
    """
    # First generate stepped surface
    stepped = generate_vicinal_surface(base_surface, 8.0, "[11-2]", 4)
    
    if not stepped["success"]:
        return stepped
    
    atoms = stepped["structure"]["atoms"]
    lattice = stepped["structure"]["lattice"]
    
    # Find step edge atoms and add kinks
    z_coords = [a["coords"][2] for a in atoms]
    unique_z = sorted(set(z_coords))
    
    n_kinks = 0
    if len(unique_z) > 2:
        # Find atoms at step edges (highest z per terrace region)
        for i, atom in enumerate(atoms):
            z = atom["coords"][2]
            x = atom["coords"][0]
            
            # Randomly add kinks at step edges
            np.random.seed(42 + i)
            if np.random.random() < kink_density:
                # Shift atom to create kink
                if kink_type == "monatomic":
                    atom["coords"][2] += 0.02
                else:
                    atom["coords"][2] += 0.04
                atom["is_kink"] = True
                n_kinks += 1
    
    return {
        "success": True,
        "base_surface": base_surface,
        "kink_density": kink_density,
        "kink_type": kink_type,
        "n_kink_sites": n_kinks,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": lattice}
    }


def generate_surface_with_adterrace(
    base_surface: str = "Pt_111",
    adterrace_size: int = 10,
    adterrace_height: int = 1
) -> Dict[str, Any]:
    """
    Generate surface with ad-terrace (island).
    
    Args:
        base_surface: Base surface
        adterrace_size: Number of atoms in ad-terrace
        adterrace_height: Height of ad-terrace in layers
    
    Returns:
        Surface with ad-terrace
    """
    # Generate base slab
    slab = generate_slab(base_surface, min_slab_thickness=12.0, min_vacuum=20.0)
    
    if not slab["success"]:
        return slab
    
    atoms = slab["structure"]["atoms"]
    lattice_info = slab["structure"]["lattice"]
    
    # Find top surface atoms
    z_coords = [a["coords"][2] for a in atoms]
    z_max = max(z_coords)
    
    # Add ad-terrace atoms above surface
    config = SLAB_CONFIGURATIONS[base_surface]
    element = config["element"]
    a = config["a"]
    
    d_interlayer = a / np.sqrt(3) / lattice_info["c"]  # Approximate interlayer spacing
    
    # Create hexagonal island
    n_ring = int(np.sqrt(adterrace_size))
    for i in range(n_ring):
        for j in range(n_ring):
            if (i - n_ring/2)**2 + (j - n_ring/2)**2 <= (n_ring/2)**2:
                for h in range(adterrace_height):
                    x = 0.5 + (i - n_ring/2) * 0.05
                    y = 0.5 + (j - n_ring/2) * 0.05
                    z = z_max + d_interlayer * (h + 1)
                    
                    atoms.append({
                        "element": element,
                        "coords": [x, y, z],
                        "ad_terrace": True
                    })
    
    adterrace_atoms = len(atoms) - slab["n_atoms"]
    
    return {
        "success": True,
        "base_surface": base_surface,
        "adterrace_size": adterrace_atoms,
        "adterrace_height": adterrace_height,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": lattice_info}
    }
