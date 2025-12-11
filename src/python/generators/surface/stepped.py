"""
surface/stepped.py - Stepped and Vicinal Surfaces

Generates stepped surfaces with terraces and kinks:
- Vicinal surfaces with controlled step density
- Stepped surfaces with kinks
- Ad-terraces and ad-islands
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Common vicinal surface database
VICINAL_SURFACES = {
    # Pt vicinal surfaces
    "Pt(211)": {"base": "(111)", "step": "(100)", "terrace_width": 3, "step_height": 1, "a": 3.92},
    "Pt(533)": {"base": "(111)", "step": "(100)", "terrace_width": 4, "step_height": 1, "a": 3.92},
    "Pt(755)": {"base": "(111)", "step": "(100)", "terrace_width": 6, "step_height": 1, "a": 3.92},
    "Pt(332)": {"base": "(111)", "step": "(111)", "terrace_width": 5, "step_height": 1, "a": 3.92},
    "Pt(997)": {"base": "(111)", "step": "(111)", "terrace_width": 8, "step_height": 1, "a": 3.92},
    
    # Cu vicinal surfaces
    "Cu(211)": {"base": "(111)", "step": "(100)", "terrace_width": 3, "step_height": 1, "a": 3.61},
    "Cu(311)": {"base": "(100)", "step": "(111)", "terrace_width": 2, "step_height": 1, "a": 3.61},
    "Cu(511)": {"base": "(100)", "step": "(111)", "terrace_width": 3, "step_height": 1, "a": 3.61},
    
    # Au vicinal surfaces
    "Au(211)": {"base": "(111)", "step": "(100)", "terrace_width": 3, "step_height": 1, "a": 4.08},
    "Au(788)": {"base": "(111)", "step": "(111)", "terrace_width": 7, "step_height": 1, "a": 4.08},
    
    # Ni vicinal surfaces
    "Ni(211)": {"base": "(111)", "step": "(100)", "terrace_width": 3, "step_height": 1, "a": 3.52},
    "Ni(533)": {"base": "(111)", "step": "(100)", "terrace_width": 4, "step_height": 1, "a": 3.52},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_vicinal_surface(
    surface: str = "Pt(211)",
    n_terraces: int = 2,
    terrace_layers: int = 3,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate vicinal (stepped) surface.
    
    Args:
        surface: Vicinal surface name
        n_terraces: Number of terrace repeats
        terrace_layers: Substrate layers per terrace
        vacuum: Vacuum thickness
    
    Returns:
        Vicinal surface structure
    """
    if surface not in VICINAL_SURFACES:
        return {
            "success": False,
            "error": {"code": "INVALID_SURFACE", "message": f"Unknown surface '{surface}'",
                      "available": list(VICINAL_SURFACES.keys())}
        }
    
    info = VICINAL_SURFACES[surface]
    a = info["a"]
    terrace_width = info["terrace_width"]
    step_height = info["step_height"]
    element = surface.split("(")[0]
    
    # Build stepped structure
    a_terrace = a / np.sqrt(2)
    z_step = a / np.sqrt(3)
    
    # Total dimensions
    a_total = terrace_width * a_terrace * n_terraces
    b_total = a / np.sqrt(2) * 4  # Periodic in step direction
    c_total = (n_terraces * step_height) * z_step + terrace_layers * z_step + vacuum
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    
    for terrace_idx in range(n_terraces):
        x_offset = terrace_idx * terrace_width * a_terrace / a_total
        z_offset = terrace_idx * step_height * z_step / c_total
        
        # Build terrace
        for layer in range(terrace_layers):
            z_layer = (layer * z_step + vacuum / 2) / c_total + z_offset
            
            for i in range(terrace_width):
                for j in range(4):
                    x = (i * a_terrace / a_total) + x_offset + (layer % 2) * 0.5 * a_terrace / a_total
                    y = j / 4 + (layer % 2) * 0.125
                    
                    species.append(element)
                    coords.append([x % 1, y % 1, z_layer % 0.9])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate step density
    step_density = 1 / (terrace_width * a_terrace * 1e-8) / 1e6  # steps/µm
    miscut_angle = np.degrees(np.arctan(z_step / (terrace_width * a_terrace)))
    
    return {
        "success": True,
        "surface": surface,
        "element": element,
        "base_orientation": info["base"],
        "step_orientation": info["step"],
        "terrace_width_atoms": terrace_width,
        "terrace_width_angstrom": round(terrace_width * a_terrace, 2),
        "step_height_angstrom": round(z_step, 2),
        "miscut_angle_deg": round(miscut_angle, 2),
        "step_density_per_um": round(step_density, 1),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_stepped_surface_with_kinks(
    element: str = "Pt",
    base_facet: str = "111",
    terrace_width: int = 5,
    kink_spacing: int = 3,
    n_layers: int = 4,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate stepped surface with kinks.
    
    Args:
        element: Metal element
        base_facet: Base terrace orientation
        terrace_width: Atoms per terrace width
        kink_spacing: Kinks every N atoms along step edge
        n_layers: Substrate layers
        vacuum: Vacuum thickness
    
    Returns:
        Stepped surface with kinks
    """
    a_elem = {"Pt": 3.92, "Au": 4.08, "Cu": 3.61, "Ni": 3.52, "Pd": 3.89, "Ag": 4.09}.get(element, 4.0)
    
    a_surf = a_elem / np.sqrt(2)
    z_sep = a_elem / np.sqrt(3)
    
    # Build stepped + kinked supercell
    n_step_repeat = 2
    n_kink_repeat = max(2, kink_spacing)
    
    a_total = terrace_width * a_surf * n_step_repeat
    b_total = n_kink_repeat * a_surf
    c_total = n_layers * z_sep + z_sep + vacuum
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c_total
        
        for step in range(n_step_repeat):
            for i in range(terrace_width):
                for j in range(n_kink_repeat):
                    # Apply kink: shift x position at kink sites
                    kink_shift = 0
                    if step == 1 and i == terrace_width - 1 and j % kink_spacing == 0:
                        kink_shift = 0.5 / terrace_width / n_step_repeat
                    
                    x = (step * terrace_width + i + 0.5) / (terrace_width * n_step_repeat) + kink_shift
                    y = (j + 0.5) / n_kink_repeat
                    
                    species.append(element)
                    coords.append([x % 1, y % 1, z_layer])
    
    # Add step edge atoms (one layer higher)
    z_step = (n_layers * z_sep + vacuum / 2) / c_total
    for j in range(n_kink_repeat):
        x = 0.5 / n_step_repeat
        y = (j + 0.5) / n_kink_repeat
        
        # Skip kink positions for vacancy
        if j % kink_spacing != 0:
            species.append(element)
            coords.append([x, y, z_step])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "base_facet": base_facet,
        "terrace_width": terrace_width,
        "kink_spacing": kink_spacing,
        "n_kinks": n_kink_repeat // kink_spacing,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_roughened_surface(
    element: str = "Au",
    facet: str = "111",
    roughness_rms: float = 2.0,
    n_layers: int = 6,
    size: List[int] = [10, 10],
    vacuum: float = 15.0,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate roughened surface with controlled roughness.
    
    Args:
        element: Surface element
        facet: Base facet
        roughness_rms: RMS roughness in Angstrom
        n_layers: Base number of layers
        size: Surface size [nx, ny]
        vacuum: Vacuum thickness
        seed: Random seed
    
    Returns:
        Roughened surface
    """
    np.random.seed(seed)
    
    a_elem = {"Au": 4.08, "Pt": 3.92, "Cu": 3.61, "Ag": 4.09, "Pd": 3.89}.get(element, 4.0)
    a_surf = a_elem / np.sqrt(2)
    z_sep = a_elem / np.sqrt(3)
    
    nx, ny = size
    
    a_total = nx * a_surf
    b_total = ny * a_surf
    c_total = (n_layers + 3) * z_sep + vacuum  # Extra space for roughness
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    
    # Generate height map with correlation
    height_map = np.random.randn(nx, ny) * roughness_rms / z_sep
    
    # Smooth to add correlation
    from scipy.ndimage import gaussian_filter
    height_map = gaussian_filter(height_map, sigma=1.5)
    height_map = height_map * roughness_rms / (np.std(height_map) * z_sep + 0.001)
    
    for i in range(nx):
        for j in range(ny):
            local_height = int(round(height_map[i, j]))
            local_n_layers = n_layers + local_height
            
            for layer in range(max(1, local_n_layers)):
                z_layer = (layer * z_sep + vacuum / 2) / c_total
                
                x = (i + 0.5) / nx
                y = (j + 0.5) / ny
                
                species.append(element)
                coords.append([x, y, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "facet": facet,
        "roughness_rms_angstrom": roughness_rms,
        "size": size,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
