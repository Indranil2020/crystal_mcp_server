"""
two_d/he_2d.py - Heterostructure 2D Materials

Comprehensive 2D heterostructure generation:
- Vertical van der Waals heterostructures
- Lateral heterostructures with junctions
- Moiré heterostructures
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# 2D material database for heterostructures
MATERIAL_2D_DATABASE = {
    "graphene": {"a": 2.46, "element": "C", "type": "semimetal", "gap_eV": 0},
    "hBN": {"a": 2.50, "elements": ["B", "N"], "type": "insulator", "gap_eV": 5.9},
    "MoS2": {"a": 3.16, "elements": ["Mo", "S"], "type": "semiconductor", "gap_eV": 1.8},
    "MoSe2": {"a": 3.29, "elements": ["Mo", "Se"], "type": "semiconductor", "gap_eV": 1.5},
    "MoTe2": {"a": 3.52, "elements": ["Mo", "Te"], "type": "semiconductor", "gap_eV": 1.0},
    "WS2": {"a": 3.15, "elements": ["W", "S"], "type": "semiconductor", "gap_eV": 2.0},
    "WSe2": {"a": 3.28, "elements": ["W", "Se"], "type": "semiconductor", "gap_eV": 1.6},
    "WTe2": {"a": 3.48, "elements": ["W", "Te"], "type": "semimetal", "gap_eV": 0},
    "phosphorene": {"a": 3.31, "b": 4.37, "element": "P", "type": "semiconductor", "gap_eV": 2.0},
    "silicene": {"a": 3.87, "element": "Si", "type": "semimetal", "gap_eV": 0},
    "germanene": {"a": 4.06, "element": "Ge", "type": "semimetal", "gap_eV": 0},
    "stanene": {"a": 4.67, "element": "Sn", "type": "semimetal", "gap_eV": 0.1},
}


# Common heterostructure combinations
HETEROSTRUCTURE_SYSTEMS = {
    "graphene_hBN": {
        "layers": ["graphene", "hBN"], "mismatch_percent": 1.8,
        "band_alignment": "type_I", "encapsulation": True
    },
    "MoS2_WSe2": {
        "layers": ["MoS2", "WSe2"], "mismatch_percent": 3.8,
        "band_alignment": "type_II", "interlayer_excitons": True
    },
    "MoSe2_WSe2": {
        "layers": ["MoSe2", "WSe2"], "mismatch_percent": 0.3,
        "band_alignment": "type_II", "moire_trapping": True
    },
    "graphene_MoS2": {
        "layers": ["graphene", "MoS2"], 
        "band_alignment": "Schottky", "contact": True
    },
    "hBN_graphene_hBN": {
        "layers": ["hBN", "graphene", "hBN"],
        "encapsulation": True, "high_mobility": True
    },
    "WS2_MoS2": {
        "layers": ["WS2", "MoS2"], "mismatch_percent": 0.3,
        "band_alignment": "type_II"
    },
}


# Interlayer distances
INTERLAYER_DISTANCES = {
    ("graphene", "graphene"): 3.35,
    ("graphene", "hBN"): 3.30,
    ("graphene", "MoS2"): 3.40,
    ("hBN", "hBN"): 3.30,
    ("MoS2", "MoS2"): 6.15,
    ("MoS2", "WSe2"): 6.30,
    ("MoSe2", "WSe2"): 6.35,
    ("WS2", "MoS2"): 6.20,
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def get_interlayer_distance(mat1: str, mat2: str) -> float:
    """Get interlayer distance between two materials."""
    key = (mat1, mat2)
    if key in INTERLAYER_DISTANCES:
        return INTERLAYER_DISTANCES[key]
    key_rev = (mat2, mat1)
    if key_rev in INTERLAYER_DISTANCES:
        return INTERLAYER_DISTANCES[key_rev]
    return 3.5  # Default


def generate_vertical_heterostructure(
    materials: List[str],
    supercell: List[int] = [5, 5],
    twist_angles: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Generate vertical van der Waals heterostructure.
    
    Args:
        materials: List of 2D materials to stack
        supercell: In-plane supercell
        twist_angles: Twist angles between layers
    
    Returns:
        Heterostructure structure
    """
    for mat in materials:
        if mat not in MATERIAL_2D_DATABASE:
            return {
                "success": False,
                "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material {mat}",
                          "available": list(MATERIAL_2D_DATABASE.keys())}
            }
    
    if twist_angles is None:
        twist_angles = [0] * (len(materials) - 1)
    
    if len(twist_angles) != len(materials) - 1:
        twist_angles = [0] * (len(materials) - 1)
    
    # Calculate lattice constant (use average)
    a_values = [MATERIAL_2D_DATABASE[m]["a"] for m in materials]
    a_avg = np.mean(a_values)
    
    # Calculate total height
    total_height = 20  # Vacuum
    for i in range(len(materials) - 1):
        total_height += get_interlayer_distance(materials[i], materials[i + 1])
    
    nx, ny = supercell
    lattice = Lattice.hexagonal(a_avg * nx, total_height)
    
    species = []
    coords = []
    
    current_z = 10  # Start with vacuum
    
    for layer_idx, mat in enumerate(materials):
        info = MATERIAL_2D_DATABASE[mat]
        
        # Get twist angle
        cum_twist = sum(twist_angles[:layer_idx]) if layer_idx > 0 else 0
        theta = np.radians(cum_twist)
        
        z_frac = current_z / total_height
        
        # Generate layer
        if mat == "graphene":
            for i in range(nx):
                for j in range(ny):
                    for basis in [[0, 0], [1/3, 1/3]]:
                        x = (i + basis[0]) / nx
                        y = (j + basis[1]) / ny
                        
                        x_rot = (x - 0.5) * np.cos(theta) - (y - 0.5) * np.sin(theta) + 0.5
                        y_rot = (x - 0.5) * np.sin(theta) + (y - 0.5) * np.cos(theta) + 0.5
                        
                        species.append("C")
                        coords.append([x_rot % 1, y_rot % 1, z_frac])
        
        elif mat == "hBN":
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
            # TMD: Metal-Chalcogen-Chalcogen
            elems = info.get("elements", [mat[:2], mat[-1]])
            
            for i in range(nx):
                for j in range(ny):
                    x = (i + 1/3) / nx
                    y = (j + 1/3) / ny
                    
                    x_rot = (x - 0.5) * np.cos(theta) - (y - 0.5) * np.sin(theta) + 0.5
                    y_rot = (x - 0.5) * np.sin(theta) + (y - 0.5) * np.cos(theta) + 0.5
                    
                    species.append(elems[0])
                    coords.append([x_rot % 1, y_rot % 1, z_frac])
                    
                    delta_z = 1.5 / total_height
                    species.append(elems[1])
                    coords.append([x_rot % 1, y_rot % 1, z_frac - delta_z])
                    species.append(elems[1])
                    coords.append([x_rot % 1, y_rot % 1, z_frac + delta_z])
        
        # Move to next layer
        if layer_idx < len(materials) - 1:
            current_z += get_interlayer_distance(materials[layer_idx], materials[layer_idx + 1])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "materials": materials,
        "n_layers": len(materials),
        "twist_angles_deg": twist_angles,
        "lattice_constant_A": round(a_avg, 3),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_lateral_heterostructure(
    material1: str = "graphene",
    material2: str = "hBN",
    interface_type: str = "zigzag",
    width_nm: float = 5.0
) -> Dict[str, Any]:
    """
    Generate lateral heterostructure with in-plane junction.
    
    Args:
        material1: First material
        material2: Second material
        interface_type: 'zigzag' or 'armchair'
        width_nm: Width of each region in nm
    
    Returns:
        Lateral heterostructure
    """
    if material1 not in MATERIAL_2D_DATABASE or material2 not in MATERIAL_2D_DATABASE:
        return {"success": False, "error": {"code": "INVALID_MATERIAL", "message": "Unknown material"}}
    
    info1 = MATERIAL_2D_DATABASE[material1]
    info2 = MATERIAL_2D_DATABASE[material2]
    
    a = (info1["a"] + info2["a"]) / 2
    width_A = width_nm * 10
    
    n_cells = int(width_A / a)
    
    lattice = Lattice.orthorhombic(a * 2 * n_cells, a * np.sqrt(3) * 3, 20)
    
    species = []
    coords = []
    
    for i in range(2 * n_cells):
        for j in range(6):
            x = (i + 0.5 * (j % 2)) / (2 * n_cells)
            y = j / 6
            z = 0.5
            
            if i < n_cells:
                elem = info1.get("element", info1.get("elements", ["C"])[0])
            else:
                elem = info2.get("element", info2.get("elements", ["B"])[0])
            
            species.append(elem)
            coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material1": material1,
        "material2": material2,
        "interface_type": interface_type,
        "width_nm": width_nm,
        "junction_type": "lateral",
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_encapsulated_structure(
    active_layer: str = "graphene",
    encapsulant: str = "hBN",
    n_encap_layers: int = 3
) -> Dict[str, Any]:
    """
    Generate encapsulated 2D material (e.g., hBN/graphene/hBN).
    
    Args:
        active_layer: Central active layer
        encapsulant: Encapsulating material
        n_encap_layers: Number of encapsulant layers on each side
    
    Returns:
        Encapsulated structure
    """
    layers = [encapsulant] * n_encap_layers + [active_layer] + [encapsulant] * n_encap_layers
    
    result = generate_vertical_heterostructure(layers)
    
    if result["success"]:
        result["active_layer"] = active_layer
        result["encapsulant"] = encapsulant
        result["is_encapsulated"] = True
        result["n_encapsulant_layers_per_side"] = n_encap_layers
    
    return result
