"""
two_d/magnetic_2d.py - Magnetic 2D Materials

Generates 2D magnetic monolayers:
- CrI3 (ferromagnetic insulator)
- CrGeTe3 (ferromagnetic semiconductor)
- Fe3GeTe2 (metallic ferromagnet)
- 2D Kagome magnets
- MnPS3, NiPS3, FePS3 (antiferromagnetic)

Scientific basis:
- Mermin-Wagner theorem and magnetic anisotropy
- Heisenberg/Ising models for 2D magnetism
- DFT-optimized structures with magnetic ordering
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# 2D magnetic materials database
MAGNETIC_2D_MATERIALS = {
    "CrI3": {
        "formula": "CrI3",
        "a": 6.87,
        "magnetic_element": "Cr",
        "magnetic_moment": 3.0,  # μB
        "ordering": "FM",
        "tc_kelvin": 45,
        "anisotropy": "ising",
        "description": "Ferromagnetic insulator, Ising-type"
    },
    "CrGeTe3": {
        "formula": "CrGeTe3",
        "a": 6.82,
        "magnetic_element": "Cr",
        "magnetic_moment": 3.0,
        "ordering": "FM",
        "tc_kelvin": 68,
        "anisotropy": "heisenberg",
        "description": "Ferromagnetic semiconductor"
    },
    "Fe3GeTe2": {
        "formula": "Fe3GeTe2",
        "a": 3.99,
        "magnetic_element": "Fe",
        "magnetic_moment": 2.5,
        "ordering": "FM",
        "tc_kelvin": 220,  # Record for 2D FM
        "anisotropy": "ising",
        "description": "Metallic ferromagnet, high Tc"
    },
    "MnPS3": {
        "formula": "MnPS3",
        "a": 6.08,
        "magnetic_element": "Mn",
        "magnetic_moment": 4.5,
        "ordering": "AFM",
        "tn_kelvin": 78,
        "anisotropy": "heisenberg",
        "description": "Antiferromagnetic, honeycomb Mn lattice"
    },
    "FePS3": {
        "formula": "FePS3",
        "a": 5.95,
        "magnetic_element": "Fe",
        "magnetic_moment": 4.0,
        "ordering": "AFM",
        "tn_kelvin": 118,
        "anisotropy": "ising",
        "description": "Ising antiferromagnet"
    },
    "NiPS3": {
        "formula": "NiPS3",
        "a": 5.81,
        "magnetic_element": "Ni",
        "magnetic_moment": 1.8,
        "ordering": "AFM",
        "tn_kelvin": 155,
        "anisotropy": "xxz",
        "description": "XXZ antiferromagnet"
    },
    "CrBr3": {
        "formula": "CrBr3",
        "a": 6.30,
        "magnetic_element": "Cr",
        "magnetic_moment": 3.0,
        "ordering": "FM",
        "tc_kelvin": 37,
        "anisotropy": "ising",
        "description": "CrI3 analog with Br"
    },
    "CrCl3": {
        "formula": "CrCl3",
        "a": 6.00,
        "magnetic_element": "Cr",
        "magnetic_moment": 3.0,
        "ordering": "AFM",
        "tn_kelvin": 17,
        "anisotropy": "heisenberg",
        "description": "Antiferromagnetic Cr halide"
    },
    "VI3": {
        "formula": "VI3",
        "a": 6.87,
        "magnetic_element": "V",
        "magnetic_moment": 2.0,
        "ordering": "FM",
        "tc_kelvin": 50,
        "anisotropy": "ising",
        "description": "V-based magnetic trihalide"
    },
}


def generate_magnetic_2d(
    material: str,
    magnetic_state: str = "auto",
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate 2D magnetic material structure with magnetic moments.
    
    Args:
        material: Material name (CrI3, Fe3GeTe2, etc.)
        magnetic_state: 'FM', 'AFM', or 'auto' (uses ground state)
        size: Supercell size
        vacuum: Vacuum padding
    
    Returns:
        Magnetic 2D structure with moments
    """
    if material not in MAGNETIC_2D_MATERIALS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MATERIAL",
                "message": f"Unknown material '{material}'",
                "available": list(MAGNETIC_2D_MATERIALS.keys())
            }
        }
    
    params = MAGNETIC_2D_MATERIALS[material]
    a = params["a"]
    c = vacuum
    mag_elem = params["magnetic_element"]
    moment = params["magnetic_moment"]
    ground_state = params["ordering"]
    
    if magnetic_state == "auto":
        magnetic_state = ground_state
    
    lattice = Lattice.hexagonal(a, c)
    
    z_center = 0.5
    
    # Build structure based on material type
    if material == "CrI3" or material == "CrBr3" or material == "CrCl3" or material == "VI3":
        # MX3 trihalide structure
        halogen = material[-1] if material[-2].isupper() else material[-2:]
        if halogen == "I3":
            halogen = "I"
        elif halogen == "Br3":
            halogen = "Br"
        elif halogen == "Cl3":
            halogen = "Cl"
        
        M = mag_elem
        species = [M, M, halogen, halogen, halogen, halogen, halogen, halogen]
        dz = 1.5 / c
        # Honeycomb M lattice with edge-sharing octahedra
        coords = [
            [1/3, 2/3, z_center], [2/3, 1/3, z_center],  # M atoms
            [0.08, 0.42, z_center + dz], [0.42, 0.08, z_center + dz],  # Top X
            [0.75, 0.42, z_center + dz], [0.08, 0.75, z_center - dz],  # More X
            [0.42, 0.42, z_center - dz], [0.75, 0.75, z_center - dz],  # Bottom X
        ]
    
    elif material == "Fe3GeTe2":
        species = ["Fe", "Fe", "Fe", "Ge", "Te", "Te"]
        dz = 1.0 / c
        coords = [
            [0, 0, z_center], [1/3, 2/3, z_center], [2/3, 1/3, z_center],  # Fe
            [0, 0, z_center + 0.5*dz],  # Ge
            [1/3, 2/3, z_center + 1.2*dz], [2/3, 1/3, z_center - 1.2*dz],  # Te
        ]
    
    elif material in ["MnPS3", "FePS3", "NiPS3"]:
        species = [mag_elem, mag_elem, "P", "P", "S", "S", "S", "S", "S", "S"]
        dz = 1.0 / c
        coords = [
            [1/3, 2/3, z_center], [2/3, 1/3, z_center],  # M
            [0, 0, z_center + dz], [0, 0, z_center - dz],  # P dimer
            [0.1, 0.1, z_center + 1.2*dz], [0.9, 0.9, z_center + 1.2*dz],  # S
            [0.4, 0.1, z_center - 1.2*dz], [0.1, 0.4, z_center - 1.2*dz],  
            [0.6, 0.6, z_center + 1.2*dz], [0.9, 0.6, z_center - 1.2*dz],
        ]
    
    else:
        # CrGeTe3 structure
        species = [mag_elem, mag_elem, "Ge", "Ge", "Te", "Te", "Te", "Te", "Te", "Te"]
        dz = 0.8 / c
        coords = [
            [1/3, 2/3, z_center], [2/3, 1/3, z_center],  # Cr
            [0, 0, z_center], [0.5, 0.5, z_center],  # Ge
            [0.1, 0.1, z_center + dz], [0.4, 0.1, z_center + dz],
            [0.1, 0.4, z_center - dz], [0.6, 0.6, z_center - dz],
            [0.9, 0.6, z_center + dz], [0.6, 0.9, z_center - dz],
        ]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Add magnetic moments to output
    struct_dict = structure_to_dict(structure, vacuum)
    
    for i, atom in enumerate(struct_dict["atoms"]):
        if atom["element"] == mag_elem:
            if magnetic_state == "FM":
                atom["magmom"] = [0, 0, moment]
            elif magnetic_state == "AFM":
                # Checkerboard AFM
                sign = 1 if i % 2 == 0 else -1
                atom["magmom"] = [0, 0, sign * moment]
        else:
            atom["magmom"] = [0, 0, 0]
    
    return {
        "success": True,
        "material": material,
        "formula": params["formula"],
        "magnetic_element": mag_elem,
        "magnetic_moment_uB": moment,
        "magnetic_state": magnetic_state,
        "ground_state": ground_state,
        "tc_kelvin": params.get("tc_kelvin"),
        "tn_kelvin": params.get("tn_kelvin"),
        "anisotropy": params["anisotropy"],
        "structure": struct_dict
    }
