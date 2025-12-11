"""
defect/volume_defects.py - 3D Volume Defects

Generates volume defects:
- Voids (spherical, ellipsoidal, faceted)
- Precipitates (coherent, semi-coherent)
- Crack tips (mode I, II, III)
- Radiation damage cascades
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_void(
    element: str = "Fe",
    void_shape: str = "spherical",
    void_radius: float = 5.0,
    supercell_size: List[int] = [8, 8, 8],
    void_position: str = "center"
) -> Dict[str, Any]:
    """
    Generate structure with void.
    
    Args:
        element: Host material element
        void_shape: 'spherical', 'ellipsoidal', 'cubic', 'octahedral'
        void_radius: Void radius in Angstrom
        supercell_size: Supercell dimensions
        void_position: 'center' or 'random'
    
    Returns:
        Structure with void
    """
    # BCC lattice parameters
    a = {"Fe": 2.87, "W": 3.16, "Mo": 3.15, "Cr": 2.91, "V": 3.02, "Nb": 3.30}.get(element, 2.9)
    
    nx, ny, nz = supercell_size
    
    # Create supercell
    lattice = Lattice.cubic(a * nx)
    
    species = []
    coords = []
    
    # BCC positions
    base_positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    
    # Center of void
    if void_position == "center":
        void_center = np.array([0.5, 0.5, 0.5])
    else:
        np.random.seed(42)
        void_center = np.random.uniform(0.2, 0.8, 3)
    
    void_center_cart = void_center * np.array([a * nx, a * ny, a * nz])
    n_removed = 0
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in base_positions:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    
                    cart_pos = np.array([x, y, z]) * np.array([a * nx, a * ny, a * nz])
                    
                    # Check if inside void
                    diff = cart_pos - void_center_cart
                    
                    inside = False
                    if void_shape == "spherical":
                        inside = np.linalg.norm(diff) < void_radius
                    elif void_shape == "ellipsoidal":
                        # Ellipsoid with aspect ratio 2:1:1
                        scaled = diff / np.array([2, 1, 1])
                        inside = np.linalg.norm(scaled) < void_radius / 2
                    elif void_shape == "cubic":
                        inside = all(abs(d) < void_radius for d in diff)
                    elif void_shape == "octahedral":
                        inside = sum(abs(d) for d in diff) < void_radius * 1.5
                    
                    if inside:
                        n_removed += 1
                    else:
                        species.append(element)
                        coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate void volume
    if void_shape == "spherical":
        void_volume = 4/3 * np.pi * void_radius**3
    elif void_shape == "cubic":
        void_volume = (2 * void_radius)**3
    else:
        void_volume = 4/3 * np.pi * void_radius**3  # Approximate
    
    return {
        "success": True,
        "element": element,
        "void_shape": void_shape,
        "void_radius_angstrom": void_radius,
        "void_volume_angstrom3": round(void_volume, 1),
        "n_vacancies": n_removed,
        "n_atoms": len(structure),
        "supercell_size": supercell_size,
        "structure": structure_to_dict(structure)
    }


def generate_precipitate(
    matrix_element: str = "Al",
    precipitate_element: str = "Cu",
    precipitate_type: str = "coherent",
    precipitate_radius: float = 8.0,
    supercell_size: List[int] = [10, 10, 10]
) -> Dict[str, Any]:
    """
    Generate structure with precipitate.
    
    Args:
        matrix_element: Matrix material
        precipitate_element: Precipitate material
        precipitate_type: 'coherent', 'semi-coherent', 'incoherent'
        precipitate_radius: Precipitate radius
        supercell_size: Supercell size
    
    Returns:
        Structure with precipitate
    """
    # FCC lattice
    a_matrix = {"Al": 4.05, "Cu": 3.61, "Ni": 3.52, "Ag": 4.09, "Au": 4.08}.get(matrix_element, 4.0)
    a_precip = {"Cu": 3.61, "Al": 4.05, "Ni": 3.52, "Ti": 4.08}.get(precipitate_element, 3.8)
    
    nx, ny, nz = supercell_size
    
    lattice = Lattice.cubic(a_matrix * nx)
    
    species = []
    coords = []
    
    # FCC positions
    fcc_positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    precip_center = np.array([0.5, 0.5, 0.5]) * np.array([a_matrix * nx, a_matrix * ny, a_matrix * nz])
    
    n_matrix = 0
    n_precip = 0
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for fp in fcc_positions:
                    x = (i + fp[0]) / nx
                    y = (j + fp[1]) / ny
                    z = (k + fp[2]) / nz
                    
                    cart_pos = np.array([x, y, z]) * np.array([a_matrix * nx, a_matrix * ny, a_matrix * nz])
                    
                    dist = np.linalg.norm(cart_pos - precip_center)
                    
                    if dist < precipitate_radius:
                        species.append(precipitate_element)
                        n_precip += 1
                    else:
                        species.append(matrix_element)
                        n_matrix += 1
                    
                    coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate interface properties
    lattice_mismatch = (a_precip - a_matrix) / a_matrix * 100
    
    return {
        "success": True,
        "matrix_element": matrix_element,
        "precipitate_element": precipitate_element,
        "precipitate_type": precipitate_type,
        "precipitate_radius_angstrom": precipitate_radius,
        "lattice_mismatch_percent": round(lattice_mismatch, 2),
        "n_matrix_atoms": n_matrix,
        "n_precipitate_atoms": n_precip,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_crack_tip(
    element: str = "Fe",
    crack_mode: str = "I",
    crack_length: float = 10.0,
    supercell_size: List[int] = [20, 10, 4]
) -> Dict[str, Any]:
    """
    Generate structure with crack tip.
    
    Args:
        element: Material element
        crack_mode: 'I' (opening), 'II' (sliding), 'III' (tearing)
        crack_length: Crack length in Angstrom
        supercell_size: Supercell size [x, y, z]
    
    Returns:
        Structure with crack
    """
    a = {"Fe": 2.87, "Al": 4.05, "Cu": 3.61, "W": 3.16}.get(element, 3.0)
    
    nx, ny, nz = supercell_size
    
    lattice = Lattice.orthorhombic(a * nx, a * ny, a * nz)
    
    species = []
    coords = []
    
    # BCC positions
    base_positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    
    # Crack is along x, opening in y
    crack_tip_x = 0.3  # Fractional position of crack tip
    crack_y_center = 0.5
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in base_positions:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    
                    # Check if in crack region
                    in_crack = False
                    
                    if x < crack_tip_x:
                        # Before crack tip - apply crack opening
                        crack_half_width = crack_length / (a * nx) * 0.1  # Small opening
                        if abs(y - crack_y_center) < crack_half_width:
                            in_crack = True
                    
                    if not in_crack:
                        # Apply displacement field based on crack mode
                        if crack_mode == "I" and x < crack_tip_x:
                            # Opening mode - displaces in y
                            if y < crack_y_center:
                                y -= 0.01 * (crack_tip_x - x)
                            else:
                                y += 0.01 * (crack_tip_x - x)
                        
                        species.append(element)
                        coords.append([max(0, min(1, x)), max(0, min(1, y)), z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "crack_mode": crack_mode,
        "crack_mode_description": {"I": "opening", "II": "sliding", "III": "tearing"}[crack_mode],
        "crack_length_angstrom": crack_length,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_radiation_cascade(
    element: str = "Fe",
    pka_energy_keV: float = 10.0,
    supercell_size: List[int] = [15, 15, 15],
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate radiation damage cascade structure.
    
    Args:
        element: Material element
        pka_energy_keV: Primary knock-on atom energy in keV
        supercell_size: Supercell size
        seed: Random seed
    
    Returns:
        Structure with cascade damage
    """
    np.random.seed(seed)
    
    a = {"Fe": 2.87, "W": 3.16, "Cu": 3.61, "Ni": 3.52}.get(element, 3.0)
    
    nx, ny, nz = supercell_size
    n_atoms_perfect = nx * ny * nz * 2  # BCC
    
    # NRT model for displacement estimate
    E_d = 40  # eV, displacement threshold
    n_displacements = int(0.8 * pka_energy_keV * 1000 / (2 * E_d))
    n_displacements = min(n_displacements, n_atoms_perfect // 4)
    
    lattice = Lattice.cubic(a * nx)
    
    all_sites = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in [[0, 0, 0], [0.5, 0.5, 0.5]]:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    all_sites.append([x, y, z])
    
    # Select cascade center
    cascade_center = np.array([0.5, 0.5, 0.5])
    
    # Sort sites by distance from cascade center
    distances = [np.linalg.norm(np.array(s) - cascade_center) for s in all_sites]
    sorted_indices = np.argsort(distances)
    
    # Remove atoms near cascade (vacancies)
    vacancy_indices = set(sorted_indices[:n_displacements])
    
    species = []
    coords = []
    
    for i, site in enumerate(all_sites):
        if i not in vacancy_indices:
            species.append(element)
            coords.append(site)
    
    # Some displaced atoms become interstitials (simplified)
    n_interstitials = n_displacements // 2
    for _ in range(n_interstitials):
        # Random interstitial position near cascade
        pos = cascade_center + np.random.normal(0, 0.1, 3)
        pos = pos % 1  # Wrap to unit cell
        species.append(element)
        coords.append(list(pos))
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "pka_energy_keV": pka_energy_keV,
        "n_displacements_nrt": n_displacements,
        "n_vacancies": n_displacements,
        "n_interstitials": n_interstitials,
        "n_atoms_initial": n_atoms_perfect,
        "n_atoms_final": len(structure),
        "structure": structure_to_dict(structure)
    }
