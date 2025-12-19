"""
twist/multilayers.py - Multilayer Twisted Structures

Comprehensive multilayer twisted structure generation:
- Trilayer graphene (alternating twist)
- Multilayer TMD stacks
- Mixed material heterostructures
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Multilayer system database
MULTILAYER_SYSTEMS = {
    # Graphene multilayers
    "ABA_graphene": {
        "n_layers": 3, "stacking": "ABA",
        "twist_angles": [0, 0], "Bernal": True,
        "description": "Bernal stacked trilayer"
    },
    "ABC_graphene": {
        "n_layers": 3, "stacking": "ABC",
        "twist_angles": [0, 0], "rhombohedral": True,
        "description": "Rhombohedral stacked trilayer"
    },
    "twisted_trilayer_ALT": {
        "n_layers": 3, "stacking": "twisted",
        "twist_angles": [1.57, -1.57], "alternating": True,
        "magic_angle": True, "FQHE": True,
        "description": "Alternating twist trilayer graphene"
    },
    "twisted_trilayer_MONO": {
        "n_layers": 3, "stacking": "twisted",
        "twist_angles": [1.08, 0], "mono_twist": True,
        "description": "Single twist trilayer"
    },
    "twisted_double_bilayer": {
        "n_layers": 4, "stacking": "AB_AB",
        "twist_angles": [0, 1.3, 0], 
        "correlated": True,
        "description": "Twisted double bilayer graphene"
    },
    
    # TMD multilayers
    "2H_MoS2_trilayer": {
        "n_layers": 3, "material": "MoS2",
        "stacking": "2H", "twist_angles": [0, 0],
        "description": "2H stacked MoS2 trilayer"
    },
    "3R_MoS2_trilayer": {
        "n_layers": 3, "material": "MoS2",
        "stacking": "3R", "twist_angles": [0, 0],
        "description": "3R stacked MoS2 trilayer"
    },
    "twisted_WSe2_multilayer": {
        "n_layers": 3, "material": "WSe2",
        "stacking": "twisted", "twist_angles": [2.0, 2.0],
        "moire_excitons": True
    },
    
    # Mixed heterostructures
    "graphene_hBN_graphene": {
        "n_layers": 3, "materials": ["graphene", "hBN", "graphene"],
        "twist_angles": [0, 0], "encapsulated": True,
        "description": "hBN encapsulated graphene"
    },
    "MoSe2_WSe2_MoSe2": {
        "n_layers": 3, "materials": ["MoSe2", "WSe2", "MoSe2"],
        "twist_angles": [0, 0], "type_II": True,
        "description": "TMD quantum well"
    },
}


# Interlayer distance database (Angstrom)
INTERLAYER_DISTANCES = {
    "graphene_graphene": 3.35,
    "graphene_hBN": 3.30,
    "MoS2_MoS2": 6.15,
    "WSe2_WSe2": 6.49,
    "MoSe2_WSe2": 6.35,
    "MoS2_WSe2": 6.30,
    "graphene_MoS2": 3.40,
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def get_interlayer_distance(material1: str, material2: str) -> float:
    """Get interlayer distance between two materials."""
    key = f"{material1}_{material2}"
    if key in INTERLAYER_DISTANCES:
        return INTERLAYER_DISTANCES[key]
    
    key_rev = f"{material2}_{material1}"
    if key_rev in INTERLAYER_DISTANCES:
        return INTERLAYER_DISTANCES[key_rev]
    
    # Default
    return 3.5


def generate_multilayer_stack(
    system: str = "twisted_trilayer_ALT",
    supercell: List[int] = [5, 5]
) -> Dict[str, Any]:
    """
    Generate multilayer twisted structure.
    
    Args:
        system: Multilayer system from database
        supercell: In-plane supercell
    
    Returns:
        Multilayer structure
    """
    if system not in MULTILAYER_SYSTEMS:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", "message": f"Unknown system",
                      "available": list(MULTILAYER_SYSTEMS.keys())}
        }
    
    info = MULTILAYER_SYSTEMS[system]
    n_layers = info["n_layers"]
    twist_angles = info["twist_angles"]
    
    # Determine material and lattice constant
    if "materials" in info:
        materials = info["materials"]
        a = 2.46  # Default graphene
    elif "material" in info:
        materials = [info["material"]] * n_layers
        a = 3.16 if "MoS2" in info["material"] else 3.28  # TMD
    else:
        materials = ["graphene"] * n_layers
        a = 2.46
    
    # Calculate total height
    total_d = 0
    for i in range(n_layers - 1):
        d = get_interlayer_distance(materials[i], materials[i+1])
        total_d += d
    
    nx, ny = supercell
    lattice = Lattice.hexagonal(a * nx, total_d + 20)
    
    species = []
    coords = []
    
    current_z = 10  # Start above bottom vacuum
    
    for layer_idx in range(n_layers):
        material = materials[layer_idx]
        
        # Get twist angle for this layer
        if layer_idx == 0:
            twist = 0
        else:
            twist = sum(twist_angles[:layer_idx])
        
        theta = np.radians(twist)
        
        z_frac = current_z / (total_d + 20)
        
        # Generate layer atoms
        if "graphene" in material.lower() or material == "C":
            # Graphene: 2 atoms per unit cell
            for i in range(nx):
                for j in range(ny):
                    for basis in [[0, 0], [1/3, 1/3]]:
                        x = (i + basis[0]) / nx
                        y = (j + basis[1]) / ny
                        
                        # Apply twist
                        x_rot = (x - 0.5) * np.cos(theta) - (y - 0.5) * np.sin(theta) + 0.5
                        y_rot = (x - 0.5) * np.sin(theta) + (y - 0.5) * np.cos(theta) + 0.5
                        
                        species.append("C")
                        coords.append([x_rot % 1, y_rot % 1, z_frac])
        
        elif "hBN" in material or "BN" in material:
            for i in range(nx):
                for j in range(ny):
                    for basis_idx, basis in enumerate([[0, 0], [1/3, 1/3]]):
                        x = (i + basis[0]) / nx
                        y = (j + basis[1]) / ny
                        
                        x_rot = (x - 0.5) * np.cos(theta) - (y - 0.5) * np.sin(theta) + 0.5
                        y_rot = (x - 0.5) * np.sin(theta) + (y - 0.5) * np.cos(theta) + 0.5
                        
                        species.append("B" if basis_idx == 0 else "N")
                        coords.append([x_rot % 1, y_rot % 1, z_frac])
        
        else:
            # TMD: M-X-X sandwich
            for i in range(nx):
                for j in range(ny):
                    x = (i + 1/3) / nx
                    y = (j + 1/3) / ny
                    
                    x_rot = (x - 0.5) * np.cos(theta) - (y - 0.5) * np.sin(theta) + 0.5
                    y_rot = (x - 0.5) * np.sin(theta) + (y - 0.5) * np.cos(theta) + 0.5
                    
                    # Metal center
                    metal = material[:2]
                    species.append(metal)
                    coords.append([x_rot % 1, y_rot % 1, z_frac])
                    
                    # Chalcogen above and below
                    chalc = material[-1] if material[-1] in ["S", "Se", "Te"] else "S"
                    delta_z = 0.02
                    species.append(chalc)
                    coords.append([x_rot % 1, y_rot % 1, z_frac - delta_z])
                    species.append(chalc)
                    coords.append([x_rot % 1, y_rot % 1, z_frac + delta_z])
        
        # Move to next layer
        if layer_idx < n_layers - 1:
            current_z += get_interlayer_distance(materials[layer_idx], materials[layer_idx + 1])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "system": system,
        "n_layers": n_layers,
        "materials": materials,
        "twist_angles_deg": twist_angles,
        "stacking": info.get("stacking", "unknown"),
        "is_magic_angle": info.get("magic_angle", False),
        "is_correlated": info.get("correlated", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_heterostructure_stack(
    materials: List[str],
    twist_angles: List[float],
    supercell: List[int] = [5, 5]
) -> Dict[str, Any]:
    """
    Generate custom heterostructure stack.
    
    Args:
        materials: List of materials for each layer
        twist_angles: Twist angles between layers
        supercell: In-plane supercell
    
    Returns:
        Custom heterostructure
    """
    if len(twist_angles) != len(materials) - 1:
        return {
            "success": False,
            "error": {"code": "LENGTH_MISMATCH", 
                      "message": f"Need {len(materials)-1} angles for {len(materials)} layers"}
        }
    
    # Create a custom system entry
    system_key = "_".join(materials)
    MULTILAYER_SYSTEMS[system_key] = {
        "n_layers": len(materials),
        "materials": materials,
        "twist_angles": twist_angles,
        "stacking": "custom"
    }
    
    result = generate_multilayer_stack(system_key, supercell)
    
    # Clean up
    del MULTILAYER_SYSTEMS[system_key]
    
    return result
