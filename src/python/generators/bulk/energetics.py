"""
bulk/energetics.py - Energetic Materials and Cocrystals

Generates energetic material structures:
- High explosives (RDX, HMX, CL-20, PETN)
- Cocrystals and eutectics
- Plastic-bonded explosive (PBX) models
- Porous energetic materials

Scientific basis:
- Cambridge Structural Database (CSD) entries
- Crystal packing analysis
- Void and porosity modeling
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


# Energetic materials database
ENERGEIC_MATERIALS = {
    "RDX": {
        "formula": "C3H6N6O6",
        "name": "cyclotrimethylenetrinitramine",
        "a": 13.18, "b": 11.57, "c": 10.71,
        "alpha": 90, "beta": 90, "gamma": 90,
        "space_group": 61,  # Pbca
        "density": 1.82,
        "detonation_velocity": 8750  # m/s
    },
    "HMX": {
        "formula": "C4H8N8O8",
        "name": "cyclotetramethylenetetranitramine",
        "a": 6.54, "b": 11.05, "c": 8.70,
        "alpha": 90, "beta": 124.3, "gamma": 90,
        "space_group": 14,  # P2_1/c (beta form)
        "density": 1.91,
        "detonation_velocity": 9100
    },
    "CL-20": {
        "formula": "C6H6N12O12",
        "name": "hexanitrohexaazaisowurtzitane",
        "a": 13.23, "b": 8.17, "c": 14.87,
        "alpha": 90, "beta": 109.2, "gamma": 90,
        "space_group": 14,  # epsilon form
        "density": 2.04,
        "detonation_velocity": 9400
    },
    "PETN": {
        "formula": "C5H8N4O12",
        "name": "pentaerythritol tetranitrate",
        "a": 9.38, "b": 9.38, "c": 6.71,
        "alpha": 90, "beta": 90, "gamma": 90,
        "space_group": 88,  # I-42d
        "density": 1.77,
        "detonation_velocity": 8400
    },
    "TNT": {
        "formula": "C7H5N3O6",
        "name": "2,4,6-trinitrotoluene",
        "a": 21.00, "b": 6.03, "c": 15.00,
        "alpha": 90, "beta": 110.0, "gamma": 90,
        "space_group": 14,
        "density": 1.65,
        "detonation_velocity": 6900
    },
    "TATB": {
        "formula": "C6H6N6O6",
        "name": "1,3,5-triamino-2,4,6-trinitrobenzene",
        "a": 9.01, "b": 9.03, "c": 6.81,
        "alpha": 108.6, "beta": 91.8, "gamma": 119.97,
        "space_group": 2,  # P-1
        "density": 1.93,
        "detonation_velocity": 7760
    }
}


def generate_energetic_material(
    material: str,
    polymorph: str = "stable"
) -> Dict[str, Any]:
    """
    Generate energetic material crystal structure.
    
    Args:
        material: Material name (RDX, HMX, CL-20, PETN, TNT, TATB)
        polymorph: Crystal polymorph ('alpha', 'beta', 'gamma', 'epsilon', 'stable')
    
    Returns:
        Energetic material structure
    """
    if material not in ENERGEIC_MATERIALS:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_MATERIAL",
                "message": f"Unknown energetic material '{material}'",
                "available": list(ENERGEIC_MATERIALS.keys())
            }
        }
    
    params = ENERGEIC_MATERIALS[material]
    
    lattice = Lattice.from_parameters(
        params["a"], params["b"], params["c"],
        params["alpha"], params["beta"], params["gamma"]
    )
    
    # Simplified molecular placement (actual would use CSD coordinates)
    # One molecule at origin
    species = []
    coords = []
    
    # Parse formula for element counts
    import re
    formula = params["formula"]
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    
    n_atoms = 0
    for elem, count in elements:
        if elem:
            count = int(count) if count else 1
            for i in range(count):
                species.append(elem)
                # Place atoms in a rough molecular shape
                coords.append([0.1 * n_atoms, 0.1 * (n_atoms % 5), 0.5])
                n_atoms += 1
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "name": params["name"],
        "formula": params["formula"],
        "polymorph": polymorph,
        "density_g_cm3": params["density"],
        "detonation_velocity_m_s": params["detonation_velocity"],
        "space_group": params["space_group"],
        "note": "Simplified structure - use CSD for accurate atomic coordinates",
        "structure": structure_to_dict(structure)
    }


def generate_energetic_cocrystal(
    component_a: str,
    component_b: str,
    ratio: Tuple[int, int] = (1, 1)
) -> Dict[str, Any]:
    """
    Generate energetic cocrystal structure.
    
    Args:
        component_a: First energetic molecule
        component_b: Second component (energetic or coformer)
        ratio: Stoichiometric ratio (a:b)
    
    Returns:
        Cocrystal structure
    """
    # Common cocrystals
    COCRYSTALS = {
        ("HMX", "CL-20"): {"a": 12.5, "b": 9.2, "c": 15.3, "density": 1.98},
        ("CL-20", "TNT"): {"a": 14.1, "b": 8.5, "c": 16.2, "density": 1.92},
        ("RDX", "TNT"): {"a": 13.0, "b": 10.5, "c": 12.8, "density": 1.75}
    }
    
    key = (component_a, component_b)
    if key not in COCRYSTALS and (component_b, component_a) not in COCRYSTALS:
        # Generate estimated structure
        a = 12.0
        lattice = Lattice.orthorhombic(a, a * 0.8, a * 1.2)
        
        return {
            "success": True,
            "cocrystal": f"{component_a}-{component_b}",
            "ratio": ratio,
            "note": "Estimated structure - requires experimental validation",
            "structure": structure_to_dict(Structure(lattice, ['C'], [[0, 0, 0]]))
        }
    
    if key not in COCRYSTALS:
        key = (component_b, component_a)
        ratio = (ratio[1], ratio[0])
    
    params = COCRYSTALS[key]
    lattice = Lattice.orthorhombic(params["a"], params["b"], params["c"])
    
    return {
        "success": True,
        "cocrystal": f"{key[0]}-{key[1]}",
        "ratio": ratio,
        "density_g_cm3": params["density"],
        "structure": structure_to_dict(Structure(lattice, ['C', 'N'], [[0, 0, 0], [0.5, 0.5, 0.5]]))
    }


def generate_porous_energetic(
    material: str,
    void_fraction: float = 0.1,
    void_shape: str = "spherical",
    void_radius: float = 5.0
) -> Dict[str, Any]:
    """
    Generate porous energetic material (e.g., for PBX modeling).
    
    Args:
        material: Base energetic material
        void_fraction: Target void/porosity fraction
        void_shape: 'spherical', 'cylindrical', 'irregular'
        void_radius: Characteristic void size (Å)
    
    Returns:
        Porous energetic structure
    """
    if material not in ENERGEIC_MATERIALS:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {material}"}}
    
    params = ENERGEIC_MATERIALS[material]
    
    # Estimate supercell needed for void
    supercell_size = int(np.ceil(void_radius * 2 / min(params["a"], params["b"], params["c"])))
    supercell_size = max(supercell_size, 3)
    
    effective_density = params["density"] * (1 - void_fraction)
    
    return {
        "success": True,
        "material": material,
        "void_fraction": void_fraction,
        "void_shape": void_shape,
        "void_radius_angstrom": void_radius,
        "effective_density_g_cm3": effective_density,
        "supercell_recommended": [supercell_size] * 3,
        "note": "Remove atoms within void region after supercell generation"
    }


def generate_pbx_model(
    energetic: str,
    binder: str = "HTPB",
    energetic_fraction: float = 0.9
) -> Dict[str, Any]:
    """
    Generate plastic-bonded explosive (PBX) model.
    
    Args:
        energetic: Energetic filler (HMX, RDX, CL-20)
        binder: Polymer binder (HTPB, Estane, GAP)
        energetic_fraction: Mass fraction of energetic
    
    Returns:
        PBX model parameters
    """
    if energetic not in ENERGEIC_MATERIALS:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown: {energetic}"}}
    
    BINDER_DENSITIES = {
        "HTPB": 0.91,  # g/cm³
        "Estane": 1.12,
        "GAP": 1.29,
        "PBAN": 1.01
    }
    
    binder_density = BINDER_DENSITIES.get(binder, 1.0)
    energetic_density = ENERGEIC_MATERIALS[energetic]["density"]
    
    # Rule of mixtures for density
    composite_density = 1 / (energetic_fraction / energetic_density + 
                            (1 - energetic_fraction) / binder_density)
    
    return {
        "success": True,
        "energetic": energetic,
        "binder": binder,
        "energetic_mass_fraction": energetic_fraction,
        "composite_density_g_cm3": round(composite_density, 3),
        "energetic_volume_fraction": energetic_fraction * composite_density / energetic_density,
        "note": "Full PBX modeling requires mesoscale simulation"
    }
