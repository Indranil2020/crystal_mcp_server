"""
surface/alloys.py - Surface Alloying and De-alloying

Generates surface alloy structures:
- Pt-skin on Ni (core-shell)
- Cu-rich surfaces on CuAu
- Surface segregation profiles
- De-alloyed structures with porosity
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Surface alloy database
SURFACE_ALLOY_DATABASE = {
    "Pt-skin/Ni": {
        "core": "Ni", "skin": "Pt", "a_avg": 3.72, 
        "skin_layers": 1, "description": "Pt monolayer on Ni"
    },
    "Pt-skin/Ni3": {
        "core": "Ni", "skin": "Pt", "a_avg": 3.72,
        "skin_layers": 3, "description": "3-layer Pt-enriched surface"
    },
    "Pt-skeleton/Ni": {
        "core": "Ni", "skin": "Pt", "a_avg": 3.72,
        "skeleton": True, "description": "De-alloyed Pt skeleton"
    },
    "Pd-skin/Cu": {
        "core": "Cu", "skin": "Pd", "a_avg": 3.75,
        "skin_layers": 1, "description": "Pd monolayer on Cu"
    },
    "Au-skin/Ni": {
        "core": "Ni", "skin": "Au", "a_avg": 3.80,
        "skin_layers": 1, "description": "Au monolayer on Ni"
    },
    "Cu-rich/CuAu": {
        "bulk": ["Cu", "Au"], "surface": "Cu", "a_avg": 3.85,
        "segregation": True, "description": "Cu-segregated CuAu surface"
    },
    "Ag-skin/Pd": {
        "core": "Pd", "skin": "Ag", "a_avg": 3.97,
        "skin_layers": 1, "description": "Ag monolayer on Pd"
    },
    "Pt3Ni-skin": {
        "bulk": ["Pt", "Ni"], "ratio": [3, 1], "a_avg": 3.82,
        "skin_layers": 2, "description": "Pt3Ni with Pt-enriched surface"
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_surface_alloy(
    alloy_type: str = "Pt-skin/Ni",
    size: List[int] = [3, 3],
    n_layers: int = 6,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate surface alloy structure.
    
    Args:
        alloy_type: Type of surface alloy
        size: Surface size [nx, ny]
        n_layers: Total number of layers
        vacuum: Vacuum thickness
    
    Returns:
        Surface alloy structure
    """
    if alloy_type not in SURFACE_ALLOY_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ALLOY", "message": f"Unknown alloy '{alloy_type}'",
                      "available": list(SURFACE_ALLOY_DATABASE.keys())}
        }
    
    info = SURFACE_ALLOY_DATABASE[alloy_type]
    a = info["a_avg"]
    nx, ny = size
    
    a_surf = a / np.sqrt(2)
    z_sep = a / np.sqrt(3)
    
    a_total = nx * a_surf
    b_total = ny * a_surf
    c_total = n_layers * z_sep + vacuum
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    
    skin_layers = info.get("skin_layers", 1)
    core_elem = info.get("core", info.get("bulk", ["Ni"])[0] if isinstance(info.get("bulk"), list) else "Ni")
    skin_elem = info.get("skin", info.get("surface", "Pt"))
    
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c_total
        is_surface = (layer >= n_layers - skin_layers)
        
        for i in range(nx):
            for j in range(ny):
                x = (i + 0.5 + (layer % 2) * 0.5) / nx
                y = (j + 0.5 + (layer % 2) * 0.5) / ny
                
                if is_surface:
                    elem = skin_elem
                else:
                    elem = core_elem
                
                species.append(elem)
                coords.append([x % 1, y % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    n_skin = sum(1 for s in species if s == skin_elem)
    n_core = sum(1 for s in species if s == core_elem)
    
    return {
        "success": True,
        "alloy_type": alloy_type,
        "description": info["description"],
        "skin_element": skin_elem,
        "core_element": core_elem,
        "skin_layers": skin_layers,
        "n_skin_atoms": n_skin,
        "n_core_atoms": n_core,
        "surface_composition": f"{skin_elem}100",
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_dealloyed_surface(
    parent_alloy: str = "PtNi",
    composition: Dict[str, float] = {"Pt": 0.5, "Ni": 0.5},
    dealloying_depth: int = 3,
    porosity: float = 0.3,
    size: List[int] = [5, 5],
    n_layers: int = 8,
    vacuum: float = 15.0,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate de-alloyed porous surface.
    
    Args:
        parent_alloy: Parent alloy composition
        composition: Initial composition
        dealloying_depth: Depth of dealloying in layers
        porosity: Fraction of atoms removed
        size: Surface size
        n_layers: Total layers
        vacuum: Vacuum thickness
        seed: Random seed
    
    Returns:
        De-alloyed structure with porosity
    """
    np.random.seed(seed)
    
    elements = list(composition.keys())
    fractions = list(composition.values())
    
    # More noble element remains
    noble_order = {"Au": 1, "Pt": 2, "Pd": 3, "Ag": 4, "Cu": 5, "Ni": 6, "Co": 7, "Fe": 8}
    noble_elem = min(elements, key=lambda x: noble_order.get(x, 10))
    removed_elem = [e for e in elements if e != noble_elem][0] if len(elements) > 1 else elements[0]
    
    a = {"PtNi": 3.72, "PtCo": 3.70, "AuCu": 3.85, "PdCu": 3.78}.get(parent_alloy, 3.8)
    nx, ny = size
    
    a_surf = a / np.sqrt(2)
    z_sep = a / np.sqrt(3)
    
    a_total = nx * a_surf
    b_total = ny * a_surf
    c_total = n_layers * z_sep + vacuum
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c_total
        in_dealloyed_region = (layer >= n_layers - dealloying_depth)
        
        for i in range(nx):
            for j in range(ny):
                x = (i + 0.5 + (layer % 2) * 0.5) / nx
                y = (j + 0.5 + (layer % 2) * 0.5) / ny
                
                # Assign element based on composition
                elem = np.random.choice(elements, p=fractions)
                
                # In dealloyed region, remove less noble element with probability
                if in_dealloyed_region and elem == removed_elem:
                    if np.random.random() < porosity + (n_layers - 1 - layer) * 0.1:
                        continue  # Remove this atom (porosity)
                
                species.append(elem)
                coords.append([x % 1, y % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    # Calculate final composition
    final_comp = {e: sum(1 for s in species if s == e) / len(species) for e in elements}
    
    return {
        "success": True,
        "parent_alloy": parent_alloy,
        "initial_composition": composition,
        "final_composition": {k: round(v, 3) for k, v in final_comp.items()},
        "dealloying_depth_layers": dealloying_depth,
        "porosity": porosity,
        "removed_element": removed_elem,
        "noble_element": noble_elem,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_segregation_profile(
    alloy: str = "PtNi",
    bulk_composition: float = 0.5,
    surface_enrichment: float = 0.9,
    profile_type: str = "exponential",
    decay_length: float = 2.0,
    size: List[int] = [4, 4],
    n_layers: int = 10,
    vacuum: float = 15.0,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate alloy with segregation profile.
    
    Args:
        alloy: Binary alloy
        bulk_composition: Bulk composition of first element
        surface_enrichment: Surface composition of first element
        profile_type: 'exponential', 'linear', 'step'
        decay_length: Decay length in layers
        size: Surface size
        n_layers: Total layers
        vacuum: Vacuum thickness
        seed: Random seed
    
    Returns:
        Alloy with segregation profile
    """
    np.random.seed(seed)
    
    elements = [alloy[:2], alloy[2:]]
    if alloy == "PtNi":
        elements = ["Pt", "Ni"]
    elif alloy == "PtCo":
        elements = ["Pt", "Co"]
    elif alloy == "AuCu":
        elements = ["Au", "Cu"]
    elif alloy == "AgPd":
        elements = ["Ag", "Pd"]
    
    a = {"PtNi": 3.72, "PtCo": 3.70, "AuCu": 3.85, "AgPd": 3.98}.get(alloy, 3.8)
    nx, ny = size
    
    a_surf = a / np.sqrt(2)
    z_sep = a / np.sqrt(3)
    
    a_total = nx * a_surf
    b_total = ny * a_surf
    c_total = n_layers * z_sep + vacuum
    
    lattice = Lattice.orthorhombic(a_total, b_total, c_total)
    
    species = []
    coords = []
    layer_compositions = []
    
    for layer in range(n_layers):
        z_layer = (layer * z_sep + vacuum / 2) / c_total
        depth = n_layers - 1 - layer  # Distance from surface
        
        # Calculate local composition based on profile
        if profile_type == "exponential":
            local_comp = bulk_composition + (surface_enrichment - bulk_composition) * np.exp(-depth / decay_length)
        elif profile_type == "linear":
            local_comp = bulk_composition + (surface_enrichment - bulk_composition) * max(0, 1 - depth / (decay_length * 3))
        else:  # step
            local_comp = surface_enrichment if depth < decay_length else bulk_composition
        
        layer_compositions.append(round(local_comp, 3))
        
        for i in range(nx):
            for j in range(ny):
                x = (i + 0.5 + (layer % 2) * 0.5) / nx
                y = (j + 0.5 + (layer % 2) * 0.5) / ny
                
                elem = elements[0] if np.random.random() < local_comp else elements[1]
                
                species.append(elem)
                coords.append([x % 1, y % 1, z_layer])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "alloy": alloy,
        "elements": elements,
        "bulk_composition": bulk_composition,
        "surface_composition": surface_enrichment,
        "profile_type": profile_type,
        "decay_length_layers": decay_length,
        "layer_compositions": layer_compositions,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
