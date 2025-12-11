"""
external_fields/strain_fields.py - Strain and Deformation

Generates strained and deformed structures:
- Uniaxial/biaxial strain
- Bending/ripples
- Strain gradients
"""

from typing import Dict, Any, List, Optional, Union
import numpy as np
from pymatgen.core import Structure, Lattice


# Strain profile types
STRAIN_PROFILES = {
    "uniform": {"description": "Constant strain throughout"},
    "gradient_linear": {"description": "Linear gradient along z"},
    "gradient_exponential": {"description": "Exponential decay from surface"},
    "sinusoidal": {"description": "Sinusoidal ripple pattern"},
    "gaussian": {"description": "Gaussian bump/depression"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_strained_structure(
    element: str = "Si",
    strain_xx: float = 0.02,
    strain_yy: float = 0.0,
    strain_zz: float = 0.0,
    strain_type: str = "biaxial",
    supercell: List[int] = [3, 3, 3]
) -> Dict[str, Any]:
    """
    Generate uniformly strained structure.
    
    Args:
        element: Element
        strain_xx: Strain in x direction
        strain_yy: Strain in y (for biaxial, set equal to xx)
        strain_zz: Strain in z
        strain_type: 'uniaxial', 'biaxial', 'hydrostatic'
        supercell: Supercell size
    
    Returns:
        Strained structure
    """
    a = {"Si": 5.43, "Ge": 5.66, "GaAs": 5.65, "C": 3.57, "Fe": 2.87}.get(element, 5.0)
    
    if strain_type == "biaxial":
        strain_yy = strain_xx
        # Poisson effect
        strain_zz = -0.3 * (strain_xx + strain_yy)
    elif strain_type == "hydrostatic":
        strain_yy = strain_xx
        strain_zz = strain_xx
    
    # Apply strain to lattice
    a_strained = a * (1 + strain_xx)
    b_strained = a * (1 + strain_yy)
    c_strained = a * (1 + strain_zz)
    
    nx, ny, nz = supercell
    
    lattice = Lattice.orthorhombic(a_strained * nx, b_strained * ny, c_strained * nz)
    
    # FCC/diamond structure
    species = []
    coords = []
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # FCC positions
                for pos in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
                    x = (i + pos[0]) / nx
                    y = (j + pos[1]) / ny
                    z = (k + pos[2]) / nz
                    species.append(element)
                    coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "strain_type": strain_type,
        "strain_xx": strain_xx,
        "strain_yy": strain_yy,
        "strain_zz": strain_zz,
        "unstrained_a": a,
        "strained_a": round(a_strained, 4),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_rippled_structure(
    material: str = "graphene",
    amplitude_angstrom: float = 1.0,
    wavelength_angstrom: float = 20.0,
    direction: str = "x",
    size: List[int] = [10, 10]
) -> Dict[str, Any]:
    """
    Generate rippled 2D structure.
    
    Args:
        material: 2D material
        amplitude_angstrom: Ripple amplitude
        wavelength_angstrom: Ripple wavelength
        direction: Ripple direction ('x' or 'y')
        size: Supercell size
    
    Returns:
        Rippled structure
    """
    materials = {
        "graphene": {"element": "C", "a": 2.46},
        "h-BN": {"elements": ["B", "N"], "a": 2.51},
        "MoS2": {"elements": ["Mo", "S", "S"], "a": 3.16},
    }
    
    if material not in materials:
        return {"success": False, "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material"}}
    
    info = materials[material]
    a = info["a"]
    
    nx, ny = size
    
    lx = a * np.sqrt(3) * nx
    ly = a * ny
    lz = amplitude_angstrom * 2 + 10  # Vacuum
    
    lattice = Lattice.orthorhombic(lx, ly, lz)
    
    species = []
    coords = []
    
    elements = info.get("elements", [info["element"]])
    
    for i in range(nx):
        for j in range(ny):
            x_cart = a * np.sqrt(3) * (i + 0.5)
            y_cart = a * (j + 0.5)
            
            # Apply ripple
            if direction == "x":
                z_cart = amplitude_angstrom * np.sin(2 * np.pi * x_cart / wavelength_angstrom)
            else:
                z_cart = amplitude_angstrom * np.sin(2 * np.pi * y_cart / wavelength_angstrom)
            
            z_cart += lz / 2  # Center in box
            
            for elem in elements:
                species.append(elem)
                coords.append([x_cart / lx, y_cart / ly, z_cart / lz])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "amplitude_angstrom": amplitude_angstrom,
        "wavelength_angstrom": wavelength_angstrom,
        "direction": direction,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def apply_strain_field(
    structure_dict: Dict[str, Any],
    strain_profile: str = "gradient_linear",
    max_strain: float = 0.05,
    direction: str = "z"
) -> Dict[str, Any]:
    """
    Apply strain field to existing structure.
    
    Args:
        structure_dict: Structure dictionary
        strain_profile: Profile type from STRAIN_PROFILES
        max_strain: Maximum strain magnitude
        direction: Strain gradient direction
    
    Returns:
        Structure with applied strain field
    """
    if strain_profile not in STRAIN_PROFILES:
        return {"success": False, "error": {"code": "INVALID_PROFILE", "message": f"Unknown profile"}}
    
    atoms = structure_dict.get("structure", {}).get("atoms", [])
    lattice_info = structure_dict.get("structure", {}).get("lattice", {})
    
    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}
    
    new_atoms = []
    
    dir_idx = {"x": 0, "y": 1, "z": 2}[direction]
    
    for atom in atoms:
        pos = list(atom.get("coords", [0, 0, 0]))
        coord_val = pos[dir_idx]
        
        # Calculate local strain
        if strain_profile == "gradient_linear":
            local_strain = max_strain * coord_val
        elif strain_profile == "gradient_exponential":
            local_strain = max_strain * np.exp(-coord_val * 5)
        elif strain_profile == "sinusoidal":
            local_strain = max_strain * np.sin(2 * np.pi * coord_val)
        elif strain_profile == "gaussian":
            local_strain = max_strain * np.exp(-((coord_val - 0.5) ** 2) / 0.1)
        else:
            local_strain = max_strain
        
        # Apply strain to perpendicular directions
        for i in range(3):
            if i != dir_idx:
                pos[i] *= (1 + local_strain)
        
        new_atoms.append({**atom, "coords": pos, "local_strain": local_strain})
    
    return {
        "success": True,
        "strain_profile": strain_profile,
        "max_strain": max_strain,
        "direction": direction,
        "n_atoms": len(new_atoms),
        "structure": {"lattice": lattice_info, "atoms": new_atoms}
    }
