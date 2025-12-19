"""
surface/nanoparticles.py - Nanoparticle Generation

Comprehensive nanoparticle generation per structure_catalogue.md Category 4:
(iv) Nanoparticles – Wulff shapes, cuboctahedra, truncated octahedra, icosahedra
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Nanoparticle shape database
NANOPARTICLE_SHAPES = {
    "sphere": {"n_facets": "infinite", "description": "Spherical nanoparticle"},
    "cube": {"n_facets": 6, "faces": "(100)", "description": "Cubic nanoparticle"},
    "octahedron": {"n_facets": 8, "faces": "(111)", "description": "Octahedral"},
    "cuboctahedron": {"n_facets": 14, "faces": "(100)+(111)", "description": "Truncated octahedron"},
    "truncated_cube": {"n_facets": 14, "description": "Truncated cube"},
    "truncated_octahedron": {"n_facets": 14, "faces": "(100)+(111)", "magic": True},
    "icosahedron": {"n_facets": 20, "description": "Icosahedral", "twinned": True},
    "decahedron": {"n_facets": 10, "description": "Five-fold twinned", "twinned": True},
    "rhombic_dodecahedron": {"n_facets": 12, "faces": "(110)", "description": "Rhombic dodecahedral"},
    "tetrahedron": {"n_facets": 4, "faces": "(111)", "description": "Tetrahedral"},
    "hexagonal_prism": {"n_facets": 8, "description": "Hexagonal prism (HCP)"},
}


# Surface energy ratios for Wulff construction
SURFACE_ENERGY_RATIOS = {
    "Au": {"(111)": 1.0, "(100)": 1.27, "(110)": 1.33},
    "Ag": {"(111)": 1.0, "(100)": 1.03, "(110)": 1.06},
    "Cu": {"(111)": 1.0, "(100)": 1.21, "(110)": 1.25},
    "Pt": {"(111)": 1.0, "(100)": 1.19, "(110)": 1.23},
    "Pd": {"(111)": 1.0, "(100)": 1.17, "(110)": 1.23},
    "Ni": {"(111)": 1.0, "(100)": 1.02, "(110)": 1.00},
    "Fe_bcc": {"(110)": 1.0, "(100)": 0.92, "(111)": 1.12},
}


# Magic number clusters
MAGIC_NUMBERS = {
    "icosahedron": [13, 55, 147, 309, 561, 923, 1415],
    "cuboctahedron": [13, 55, 147, 309, 561, 923, 1415],
    "truncated_octahedron": [38, 116, 260, 490, 826, 1288],
    "octahedron": [6, 19, 44, 85, 146, 231, 344],
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_nanoparticle(
    element: str = "Au",
    shape: str = "cuboctahedron",
    size_nm: float = 3.0,
    use_magic_number: bool = False
) -> Dict[str, Any]:
    """
    Generate nanoparticle structure.
    
    Args:
        element: Metal element
        shape: Nanoparticle shape
        size_nm: Approximate diameter in nm
        use_magic_number: Use closest magic number
    
    Returns:
        Nanoparticle structure
    """
    if shape not in NANOPARTICLE_SHAPES:
        return {
            "success": False,
            "error": {"code": "INVALID_SHAPE", "message": f"Unknown shape",
                      "available": list(NANOPARTICLE_SHAPES.keys())}
        }
    
    shape_info = NANOPARTICLE_SHAPES[shape]
    
    # Lattice constant
    lattice_params = {"Au": 4.08, "Ag": 4.09, "Cu": 3.61, "Pt": 3.92, 
                      "Pd": 3.89, "Ni": 3.52, "Fe": 2.87, "Al": 4.05}
    a = lattice_params.get(element, 4.0)
    
    diameter = size_nm * 10  # Convert to Angstrom
    radius = diameter / 2
    
    # Estimate number of atoms
    volume = (4/3) * np.pi * radius**3
    atoms_per_cell = 4 if element not in ["Fe"] else 2
    n_atoms_estimate = int(volume / (a**3) * atoms_per_cell)
    
    # Use magic number if requested
    if use_magic_number and shape in MAGIC_NUMBERS:
        magic = MAGIC_NUMBERS[shape]
        n_atoms_estimate = min(magic, key=lambda x: abs(x - n_atoms_estimate))
    
    # Box size
    box_size = diameter + 20
    lattice = Lattice.cubic(box_size)
    
    species = []
    coords = []
    
    # Generate atoms based on shape
    if shape == "sphere":
        n_shell = int(radius / (a / np.sqrt(2)))
        for i in range(-n_shell, n_shell + 1):
            for j in range(-n_shell, n_shell + 1):
                for k in range(-n_shell, n_shell + 1):
                    # FCC basis
                    for basis in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
                        x = (i + basis[0]) * a / np.sqrt(2)
                        y = (j + basis[1]) * a / np.sqrt(2)
                        z = (k + basis[2]) * a / np.sqrt(2)
                        
                        r = np.sqrt(x**2 + y**2 + z**2)
                        if r <= radius:
                            fx = (x + box_size/2) / box_size
                            fy = (y + box_size/2) / box_size
                            fz = (z + box_size/2) / box_size
                            species.append(element)
                            coords.append([fx, fy, fz])
    
    elif shape == "cube":
        n_layers = int(diameter / a)
        for i in range(n_layers):
            for j in range(n_layers):
                for k in range(n_layers):
                    x = (i - n_layers/2) * a
                    y = (j - n_layers/2) * a
                    z = (k - n_layers/2) * a
                    
                    fx = (x + box_size/2) / box_size
                    fy = (y + box_size/2) / box_size
                    fz = (z + box_size/2) / box_size
                    species.append(element)
                    coords.append([fx, fy, fz])
    
    elif shape in ["octahedron", "cuboctahedron", "truncated_octahedron"]:
        n_shell = int(radius / (a / np.sqrt(2)))
        
        for i in range(-n_shell, n_shell + 1):
            for j in range(-n_shell, n_shell + 1):
                for k in range(-n_shell, n_shell + 1):
                    for basis in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
                        x = (i + basis[0]) * a / np.sqrt(2)
                        y = (j + basis[1]) * a / np.sqrt(2)
                        z = (k + basis[2]) * a / np.sqrt(2)
                        
                        # Apply shape constraint
                        if shape == "octahedron":
                            constraint = abs(x) + abs(y) + abs(z) <= radius
                        elif shape == "cuboctahedron":
                            # Combined (111) and (100) constraint
                            oct_dist = abs(x) + abs(y) + abs(z)
                            cube_dist = max(abs(x), abs(y), abs(z))
                            constraint = oct_dist <= radius * 1.5 and cube_dist <= radius * 0.8
                        else:  # truncated_octahedron
                            oct_dist = abs(x) + abs(y) + abs(z)
                            cube_dist = max(abs(x), abs(y), abs(z))
                            constraint = oct_dist <= radius * 1.3 and cube_dist <= radius * 0.9
                        
                        if constraint:
                            fx = (x + box_size/2) / box_size
                            fy = (y + box_size/2) / box_size
                            fz = (z + box_size/2) / box_size
                            species.append(element)
                            coords.append([fx, fy, fz])
    
    elif shape == "icosahedron":
        # Icosahedral structure with 5-fold twins
        phi = (1 + np.sqrt(5)) / 2
        
        # Generate concentric icosahedral shells
        n_shells = max(1, int(radius / (a * 0.8)))
        
        for shell in range(1, n_shells + 1):
            r_shell = shell * a * 0.8
            n_verts = 12 * shell**2 if shell > 0 else 1
            
            if shell == 0:
                species.append(element)
                coords.append([0.5, 0.5, 0.5])
            else:
                for i in range(min(n_verts, 100)):
                    theta = np.arccos(1 - 2 * (i + 0.5) / n_verts)
                    phi_angle = np.pi * (1 + np.sqrt(5)) * i
                    
                    x = r_shell * np.sin(theta) * np.cos(phi_angle)
                    y = r_shell * np.sin(theta) * np.sin(phi_angle)
                    z = r_shell * np.cos(theta)
                    
                    fx = (x + box_size/2) / box_size
                    fy = (y + box_size/2) / box_size
                    fz = (z + box_size/2) / box_size
                    species.append(element)
                    coords.append([fx, fy, fz])
    
    else:
        # Generic spherical for other shapes
        for i in range(n_atoms_estimate):
            theta = np.arccos(1 - 2 * (i + 0.5) / n_atoms_estimate)
            phi = np.pi * (1 + np.sqrt(5)) * i
            r = radius * (i / n_atoms_estimate) ** (1/3)
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            fx = (x + box_size/2) / box_size
            fy = (y + box_size/2) / box_size
            fz = (z + box_size/2) / box_size
            species.append(element)
            coords.append([fx, fy, fz])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate surface atoms
    n_total = len(structure)
    n_surface = int(4 * np.pi * radius**2 / (a**2 / 2))  # Approximate
    
    return {
        "success": True,
        "element": element,
        "shape": shape,
        "diameter_nm": size_nm,
        "n_atoms": n_total,
        "n_surface_atoms": min(n_surface, n_total),
        "surface_fraction": round(min(n_surface, n_total) / n_total, 2),
        "is_magic_number": n_total in MAGIC_NUMBERS.get(shape, []),
        "is_twinned": shape_info.get("twinned", False),
        "exposed_facets": shape_info.get("faces", ""),
        "structure": structure_to_dict(structure)
    }


def generate_wulff_nanoparticle(
    element: str = "Au",
    size_nm: float = 5.0,
    custom_ratios: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    """
    Generate Wulff-shape nanoparticle based on surface energies.
    
    Args:
        element: Metal element
        size_nm: Target size
        custom_ratios: Custom surface energy ratios
    
    Returns:
        Wulff nanoparticle
    """
    if element not in SURFACE_ENERGY_RATIOS and custom_ratios is None:
        return {
            "success": False,
            "error": {"code": "INVALID_ELEMENT", "message": f"No surface energy data",
                      "available": list(SURFACE_ENERGY_RATIOS.keys())}
        }
    
    ratios = custom_ratios or SURFACE_ENERGY_RATIOS[element]
    
    # Determine dominant shape based on ratios
    gamma_111 = ratios.get("(111)", 1.0)
    gamma_100 = ratios.get("(100)", 1.2)
    
    ratio = gamma_100 / gamma_111
    
    if ratio > 1.4:
        shape = "octahedron"
    elif ratio < 0.8:
        shape = "cube"
    else:
        shape = "truncated_octahedron"
    
    result = generate_nanoparticle(element, shape, size_nm)
    
    if result["success"]:
        result["is_wulff"] = True
        result["surface_energy_ratios"] = ratios
        result["predicted_shape"] = shape
    
    return result
