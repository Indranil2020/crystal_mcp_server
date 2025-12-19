"""
photonic/crystals.py - Photonic Crystals

Comprehensive photonic crystal generation per structure_catalogue.md Category 13:
(iii) Photonic & phononic crystals with deliberately broken symmetries
Additional: All standard photonic crystal structures and bandgap engineering
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# 2D Photonic crystal database
PHOTONIC_2D_DATABASE = {
    "triangular_holes": {
        "symmetry": "hexagonal", "type": "air_in_dielectric",
        "bandgap": "TE", "dielectric_band": True,
        "materials": ["Si", "GaAs", "InP", "polymer"]
    },
    "square_holes": {
        "symmetry": "square", "type": "air_in_dielectric",
        "bandgap": "TE", "dielectric_band": True,
        "materials": ["Si", "SiO2"]
    },
    "triangular_rods": {
        "symmetry": "hexagonal", "type": "dielectric_in_air",
        "bandgap": "TM", "air_band": True,
        "materials": ["Al2O3", "Si"]
    },
    "square_rods": {
        "symmetry": "square", "type": "dielectric_in_air",
        "bandgap": "TM", "air_band": True,
        "materials": ["Si", "Ge"]
    },
    "honeycomb": {
        "symmetry": "hexagonal", "type": "air_in_dielectric",
        "bandgap": "complete", "materials": ["Si"]
    },
    "kagome": {
        "symmetry": "kagome", "type": "dielectric_network",
        "bandgap": "flat_band", "topological": True,
        "materials": ["Si", "graphene"]
    },
    "lieb": {
        "symmetry": "square", "type": "dielectric_network",
        "bandgap": "flat_band", "materials": ["Si"]
    },
}


# 3D Photonic crystal database
PHOTONIC_3D_DATABASE = {
    "diamond": {
        "symmetry": "FCC", "structure": "diamond",
        "bandgap": "complete_3D", "gap_percent": 15,
        "fabrication": "difficult", "materials": ["Si", "Ge"]
    },
    "inverse_opal": {
        "symmetry": "FCC", "structure": "inverse_FCC",
        "bandgap": "complete_3D", "gap_percent": 6,
        "fabrication": "self_assembly", "materials": ["Si", "TiO2", "polymer"]
    },
    "woodpile": {
        "symmetry": "FCC", "structure": "stacked_logs",
        "bandgap": "complete_3D", "gap_percent": 18,
        "fabrication": "layer_by_layer", "materials": ["Si", "GaAs"]
    },
    "yablonovite": {
        "symmetry": "FCC", "structure": "drilled_slab",
        "bandgap": "complete_3D", "gap_percent": 20,
        "fabrication": "drilling", "materials": ["Al2O3", "Si"]
    },
    "simple_cubic": {
        "symmetry": "SC", "structure": "spheres",
        "bandgap": "partial", "gap_percent": 3,
        "fabrication": "colloidal", "materials": ["polymer", "SiO2"]
    },
    "gyroid": {
        "symmetry": "BCC", "structure": "bicontinuous",
        "bandgap": "complete_3D", "gap_percent": 8,
        "fabrication": "block_copolymer", "materials": ["polymer", "Au"]
    },
}


# Defect types for bandgap engineering
PHOTONIC_DEFECTS = {
    "point_H1": {"type": "cavity", "mode": "monopole", "Q_typical": 1000},
    "point_L3": {"type": "cavity", "mode": "dipole", "Q_typical": 10000},
    "point_H0": {"type": "cavity", "mode": "monopole_modified", "Q_typical": 100000},
    "line_W1": {"type": "waveguide", "mode": "single", "dispersion": "slow_light"},
    "line_W3": {"type": "waveguide", "mode": "multimode", "bandwidth": "wide"},
    "heterostructure": {"type": "interface", "mode": "interface_state"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_2d_photonic_crystal(
    crystal_type: str = "triangular_holes",
    material: str = "Si",
    lattice_constant_nm: float = 400.0,
    hole_radius_fraction: float = 0.3,
    slab_thickness_nm: float = 220.0,
    supercell: List[int] = [5, 5]
) -> Dict[str, Any]:
    """
    Generate 2D photonic crystal slab.
    
    Args:
        crystal_type: Crystal type from database
        material: High-index material
        lattice_constant_nm: Lattice period in nm
        hole_radius_fraction: r/a ratio
        slab_thickness_nm: Slab thickness
        supercell: Supercell size
    
    Returns:
        Photonic crystal structure
    """
    if crystal_type not in PHOTONIC_2D_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown type '{crystal_type}'",
                      "available": list(PHOTONIC_2D_DATABASE.keys())}
        }
    
    info = PHOTONIC_2D_DATABASE[crystal_type]
    
    a = lattice_constant_nm * 10  # Convert to Angstrom
    r = hole_radius_fraction * a
    t = slab_thickness_nm * 10
    
    nx, ny = supercell
    
    if info["symmetry"] == "hexagonal":
        lattice = Lattice.hexagonal(a * nx, t + 200)
    else:
        lattice = Lattice.orthorhombic(a * nx, a * ny, t + 200)
    
    species = []
    coords = []
    
    # Generate slab with holes
    n_atoms_per_cell = 20
    
    for i in range(nx):
        for j in range(ny):
            # Cell center
            if info["symmetry"] == "hexagonal":
                cx = (i + (j % 2) * 0.5 + 0.5) / nx
                cy = (j * np.sqrt(3)/2 + 0.5) / ny
            else:
                cx = (i + 0.5) / nx
                cy = (j + 0.5) / ny
            
            # Add material atoms (simplified as points representing the dielectric)
            if "holes" in crystal_type or "air_in" in info["type"]:
                # Material everywhere except hole regions
                for k in range(n_atoms_per_cell):
                    angle = 2 * np.pi * k / n_atoms_per_cell
                    dist = (0.3 + 0.15) / nx  # Outside hole region
                    
                    x = cx + dist * np.cos(angle)
                    y = cy + dist * np.sin(angle)
                    z = 0.4 + 0.2 * (k % 3) / 3  # Within slab
                    
                    if 0 < x < 1 and 0 < y < 1:
                        species.append(material)
                        coords.append([x, y, z])
            else:
                # Rods: material only in rod positions
                for dz in range(5):
                    species.append(material)
                    coords.append([cx % 1, cy % 1, 0.4 + dz * 0.04])
    
    structure = Structure(lattice, species, coords)
    
    # Estimate bandgap frequency
    c = 3e8  # m/s
    n_eff = 2.5 if material == "Si" else 1.5
    f_midgap = c / (a * 1e-10 * n_eff)  # Hz
    wavelength_midgap = a * 1e-9 * n_eff / 1  # meters
    
    return {
        "success": True,
        "crystal_type": crystal_type,
        "dimension": "2D",
        "symmetry": info["symmetry"],
        "material": material,
        "lattice_constant_nm": lattice_constant_nm,
        "hole_radius_fraction": hole_radius_fraction,
        "slab_thickness_nm": slab_thickness_nm,
        "bandgap_polarization": info["bandgap"],
        "midgap_wavelength_nm": round(wavelength_midgap * 1e9, 1),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_3d_photonic_crystal(
    crystal_type: str = "woodpile",
    material: str = "Si",
    lattice_constant_um: float = 1.5,
    fill_fraction: float = 0.25,
    n_periods: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate 3D photonic crystal.
    
    Args:
        crystal_type: Crystal type from database
        material: High-index material
        lattice_constant_um: Lattice period in μm
        fill_fraction: Volume fraction of dielectric
        n_periods: Number of periods
    
    Returns:
        3D photonic crystal structure
    """
    if crystal_type not in PHOTONIC_3D_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown type",
                      "available": list(PHOTONIC_3D_DATABASE.keys())}
        }
    
    info = PHOTONIC_3D_DATABASE[crystal_type]
    
    a = lattice_constant_um * 10000  # Convert μm to Angstrom
    nx, ny, nz = n_periods
    
    lattice = Lattice.cubic(a * nx)
    
    species = []
    coords = []
    
    np.random.seed(42)
    
    if crystal_type == "woodpile":
        # Log-pile structure
        rod_width = fill_fraction * 0.3
        n_rods = 4
        
        for layer in range(nz * 4):
            z = (layer + 0.5) / (nz * 4)
            direction = layer % 2
            offset_shift = (layer // 2) % 2
            
            for rod in range(n_rods * max(nx, ny)):
                if direction == 0:
                    y = (rod + 0.5 + offset_shift * 0.5) / (n_rods * ny)
                    if y < 1:
                        for xi in range(8):
                            x = (xi + 0.5) / 8
                            species.append(material)
                            coords.append([x, y % 1, z])
                else:
                    x = (rod + 0.5 + offset_shift * 0.5) / (n_rods * nx)
                    if x < 1:
                        for yi in range(8):
                            y = (yi + 0.5) / 8
                            species.append(material)
                            coords.append([x % 1, y, z])
    
    elif crystal_type == "inverse_opal":
        # FCC spherical voids in dielectric
        fcc_sites = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        void_radius = 0.35
        
        # Fill space except at FCC sites
        n_sample = 500
        for _ in range(n_sample):
            x = np.random.random()
            y = np.random.random()
            z = np.random.random()
            
            # Check if point is outside all voids
            is_in_void = False
            for site in fcc_sites:
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            sx = site[0] + di
                            sy = site[1] + dj
                            sz = site[2] + dk
                            dist = np.sqrt((x - sx)**2 + (y - sy)**2 + (z - sz)**2)
                            if dist < void_radius:
                                is_in_void = True
                                break
            
            if not is_in_void:
                species.append(material)
                coords.append([x, y, z])
    
    elif crystal_type == "diamond":
        # Diamond structure
        basis = [
            [0, 0, 0], [0.25, 0.25, 0.25],
            [0.5, 0.5, 0], [0.75, 0.75, 0.25],
            [0.5, 0, 0.5], [0.75, 0.25, 0.75],
            [0, 0.5, 0.5], [0.25, 0.75, 0.75]
        ]
        
        for b in basis:
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        x = (b[0] + i) / nx
                        y = (b[1] + j) / ny
                        z = (b[2] + k) / nz
                        species.append(material)
                        coords.append([x % 1, y % 1, z % 1])
    
    else:
        # Simple cubic arrangement
        for i in range(nx * 3):
            for j in range(ny * 3):
                for k in range(nz * 3):
                    x = (i + 0.5) / (nx * 3)
                    y = (j + 0.5) / (ny * 3)
                    z = (k + 0.5) / (nz * 3)
                    species.append(material)
                    coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "crystal_type": crystal_type,
        "dimension": "3D",
        "symmetry": info["symmetry"],
        "material": material,
        "lattice_constant_um": lattice_constant_um,
        "fill_fraction": fill_fraction,
        "bandgap_type": info["bandgap"],
        "expected_gap_percent": info["gap_percent"],
        "fabrication_method": info["fabrication"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_photonic_defect(
    base_crystal: str = "triangular_holes",
    defect_type: str = "point_L3",
    material: str = "Si",
    lattice_constant_nm: float = 420.0
) -> Dict[str, Any]:
    """
    Generate photonic crystal with defect for cavity/waveguide.
    
    Args:
        base_crystal: Base photonic crystal type
        defect_type: Type of defect
        material: Dielectric material
        lattice_constant_nm: Lattice constant
    
    Returns:
        Photonic crystal with defect
    """
    if defect_type not in PHOTONIC_DEFECTS:
        return {
            "success": False,
            "error": {"code": "INVALID_DEFECT", "message": f"Unknown defect type",
                      "available": list(PHOTONIC_DEFECTS.keys())}
        }
    
    defect_info = PHOTONIC_DEFECTS[defect_type]
    
    # Generate base crystal
    base = generate_2d_photonic_crystal(base_crystal, material, lattice_constant_nm, supercell=[7, 7])
    
    if not base["success"]:
        return base
    
    atoms = list(base["structure"]["atoms"])
    
    # Apply defect
    center = [0.5, 0.5]
    
    if "point" in defect_type:
        # Remove atoms near center to create cavity
        if "L3" in defect_type:
            # Line of 3 holes
            defect_width = 0.15
            filtered = [a for a in atoms if abs(a["coords"][1] - 0.5) > 0.02 or 
                       abs(a["coords"][0] - 0.5) > defect_width]
        elif "H1" in defect_type:
            # Single hole
            defect_radius = 0.05
            filtered = [a for a in atoms if 
                       np.sqrt((a["coords"][0] - 0.5)**2 + (a["coords"][1] - 0.5)**2) > defect_radius]
        else:
            filtered = atoms
        atoms = filtered
    
    elif "line" in defect_type:
        # Remove row of holes for waveguide
        if "W1" in defect_type:
            atoms = [a for a in atoms if abs(a["coords"][1] - 0.5) > 0.03]
        else:
            atoms = [a for a in atoms if abs(a["coords"][1] - 0.5) > 0.08]
    
    return {
        "success": True,
        "base_crystal": base_crystal,
        "defect_type": defect_type,
        "defect_mode": defect_info["type"],
        "mode_character": defect_info.get("mode", ""),
        "typical_Q_factor": defect_info.get("Q_typical", 0),
        "material": material,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": base["structure"]["lattice"]}
    }
