"""
external_fields/magnetic_textures.py - Magnetic Textures

Generates magnetic texture structures:
- Skyrmions (Bloch, Néel)
- Antiskyrmions
- Bimerons, bobbers
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Skyrmion host materials
SKYRMION_DATABASE = {
    "MnSi": {"type": "bulk", "skyrmion_type": "Bloch", "lattice_a": 4.558, "Tc_K": 29},
    "FeGe": {"type": "bulk", "skyrmion_type": "Bloch", "lattice_a": 4.70, "Tc_K": 278},
    "Fe/Ir(111)": {"type": "interface", "skyrmion_type": "Néel", "lattice_a": 2.72, "size_nm": 1.0},
    "Pd/Fe/Ir(111)": {"type": "multilayer", "skyrmion_type": "Néel", "lattice_a": 2.72, "size_nm": 1.5},
    "Co/Pt": {"type": "multilayer", "skyrmion_type": "Néel", "lattice_a": 2.51, "DMI": "strong"},
    "Fe3GeTe2": {"type": "2D", "skyrmion_type": "Néel", "lattice_a": 3.99, "Tc_K": 220},
}


def generate_skyrmion(
    material: str = "FeGe",
    skyrmion_radius_nm: float = 5.0,
    supercell_size: List[int] = [20, 20, 1],
    helicity: float = 0.0,
    polarity: int = 1
) -> Dict[str, Any]:
    """
    Generate structure with skyrmion spin texture.
    
    Args:
        material: Host material from database
        skyrmion_radius_nm: Skyrmion radius in nanometers
        supercell_size: Supercell size
        helicity: Skyrmion helicity (0 = Néel, π/2 = Bloch)
        polarity: Core polarity (+1 or -1)
    
    Returns:
        Structure with skyrmion magnetic configuration
    """
    if material not in SKYRMION_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material '{material}'",
                      "available": list(SKYRMION_DATABASE.keys())}
        }
    
    info = SKYRMION_DATABASE[material]
    a = info["lattice_a"]
    
    nx, ny, nz = supercell_size
    
    lattice = Lattice.cubic(a * nx)
    
    species = []
    coords = []
    spins = []
    
    center_x = 0.5
    center_y = 0.5
    radius = skyrmion_radius_nm / (a * nx) * 10  # Convert to fractional
    
    # Get magnetic element
    mag_element = "Fe"
    if "Mn" in material:
        mag_element = "Mn"
    elif "Co" in material:
        mag_element = "Co"
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x = (i + 0.5) / nx
                y = (j + 0.5) / ny
                z = (k + 0.5) / nz
                
                species.append(mag_element)
                coords.append([x, y, z])
                
                # Calculate skyrmion spin profile
                dx = x - center_x
                dy = y - center_y
                r = np.sqrt(dx**2 + dy**2)
                
                if r < 0.001:
                    # Core
                    spin = [0, 0, polarity]
                else:
                    # Radial profile (simplified)
                    phi = np.arctan2(dy, dx)
                    theta = np.pi * (1 - np.exp(-r / radius))
                    
                    # Spin direction
                    if info["skyrmion_type"] == "Bloch":
                        spin_x = np.sin(theta) * np.sin(phi + helicity)
                        spin_y = np.sin(theta) * np.cos(phi + helicity)
                    else:  # Néel
                        spin_x = np.sin(theta) * np.cos(phi + helicity)
                        spin_y = np.sin(theta) * np.sin(phi + helicity)
                    
                    spin_z = polarity * np.cos(theta)
                    spin = [spin_x, spin_y, spin_z]
                
                spins.append(spin)
    
    structure = Structure(lattice, species, coords)
    
    # Add spin data to output
    atoms_with_spins = []
    for i, s in enumerate(structure):
        atoms_with_spins.append({
            "element": str(s.specie),
            "coords": list(s.frac_coords),
            "spin": spins[i]
        })
    
    return {
        "success": True,
        "material": material,
        "skyrmion_type": info["skyrmion_type"],
        "skyrmion_radius_nm": skyrmion_radius_nm,
        "polarity": polarity,
        "helicity_rad": helicity,
        "n_atoms": len(structure),
        "topological_charge": -polarity,
        "structure": {
            "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                        "matrix": lattice.matrix.tolist()},
            "atoms": atoms_with_spins
        }
    }


def generate_antiskyrmion(
    material: str = "Mn1.4Pt0.9Pd0.1Sn",
    size_nm: float = 100.0,
    supercell_size: List[int] = [20, 20, 1]
) -> Dict[str, Any]:
    """
    Generate antiskyrmion structure.
    
    Args:
        material: Host material
        size_nm: Antiskyrmion size
        supercell_size: Supercell size
    
    Returns:
        Antiskyrmion structure
    """
    # Antiskyrmions have different symmetry (D2d vs C∞)
    nx, ny, nz = supercell_size
    a = 3.5
    
    lattice = Lattice.cubic(a * nx)
    
    species = []
    coords = []
    spins = []
    
    center_x = 0.5
    center_y = 0.5
    radius = size_nm / (a * nx) * 10
    
    for i in range(nx):
        for j in range(ny):
            x = (i + 0.5) / nx
            y = (j + 0.5) / ny
            
            species.append("Mn")
            coords.append([x, y, 0.5])
            
            dx = x - center_x
            dy = y - center_y
            r = np.sqrt(dx**2 + dy**2)
            
            if r < 0.001:
                spin = [0, 0, 1]
            else:
                phi = np.arctan2(dy, dx)
                theta = np.pi * (1 - np.exp(-r / radius))
                
                # Antiskyrmion: opposite helicity in perpendicular directions
                spin_x = np.sin(theta) * np.cos(-phi)
                spin_y = np.sin(theta) * np.sin(-phi)
                spin_z = np.cos(theta)
                spin = [spin_x, spin_y, spin_z]
            
            spins.append(spin)
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "type": "antiskyrmion",
        "size_nm": size_nm,
        "topological_charge": +1,
        "n_atoms": len(structure),
        "structure": {"lattice": {"a": a * nx}, "n_atoms": len(structure)}
    }


def generate_magnetic_bobber(
    material: str = "FeGe",
    depth_nm: float = 50.0,
    supercell_size: List[int] = [10, 10, 20]
) -> Dict[str, Any]:
    """
    Generate magnetic bobber (truncated skyrmion tube).
    
    Args:
        material: Host material
        depth_nm: Bobber depth (extinction length)
        supercell_size: Supercell size
    
    Returns:
        Bobber structure
    """
    if material not in SKYRMION_DATABASE:
        return {"success": False, "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material"}}
    
    info = SKYRMION_DATABASE[material]
    a = info["lattice_a"]
    
    nx, ny, nz = supercell_size
    
    lattice = Lattice.orthorhombic(a * nx, a * ny, a * nz)
    
    species = []
    coords = []
    spins = []
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x = (i + 0.5) / nx
                y = (j + 0.5) / ny
                z = (k + 0.5) / nz
                
                species.append("Fe")
                coords.append([x, y, z])
                
                # Bobber profile - skyrmion radius decreases with depth
                z_phys = z * a * nz
                extinction = np.exp(-z_phys / (depth_nm / 10))
                
                dx = x - 0.5
                dy = y - 0.5
                r = np.sqrt(dx**2 + dy**2)
                
                radius = 0.1 * extinction
                
                if r < radius:
                    theta = np.pi * r / radius * extinction
                else:
                    theta = np.pi
                
                spin_z = np.cos(theta)
                spins.append([0, 0, spin_z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "type": "bobber",
        "depth_nm": depth_nm,
        "n_atoms": len(structure),
        "structure": {"n_atoms": len(structure)}
    }
