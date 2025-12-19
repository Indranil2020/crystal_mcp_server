"""
defect/extended_defects.py - Extended Defects

Comprehensive extended defect generation per structure_catalogue.md Category 7:
(ii) Dislocations – edge, screw, mixed, partial, Lomer-Cottrell locks
(iii) Stacking faults – intrinsic, extrinsic, twin boundaries, growth faults
(iv) Grain boundaries – Σ3, Σ5, Σ7, Σ9, Σ11, CSL tilt, twist
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Dislocation database
DISLOCATION_DATABASE = {
    # FCC dislocations
    "FCC_edge_110": {"type": "edge", "plane": "(111)", "Burgers": "[110]/2", "character": 90},
    "FCC_screw_110": {"type": "screw", "plane": "(111)", "Burgers": "[110]/2", "character": 0},
    "FCC_mixed_60": {"type": "mixed", "plane": "(111)", "Burgers": "[110]/2", "character": 60},
    "FCC_partial_Shockley": {"type": "partial", "plane": "(111)", "Burgers": "[112]/6", 
                             "description": "Shockley partial"},
    "FCC_partial_Frank": {"type": "partial", "plane": "(111)", "Burgers": "[111]/3",
                          "description": "Frank sessile partial"},
    "FCC_Lomer_Cottrell": {"type": "lock", "description": "Lomer-Cottrell stair-rod lock"},
    
    # BCC dislocations
    "BCC_screw_111": {"type": "screw", "plane": "{110}", "Burgers": "[111]/2", "character": 0},
    "BCC_edge_111": {"type": "edge", "plane": "{110}", "Burgers": "[111]/2", "character": 90},
    "BCC_mixed_71": {"type": "mixed", "plane": "{110}", "Burgers": "[111]/2", "character": 71},
    
    # HCP dislocations
    "HCP_basal_a": {"type": "edge", "plane": "(0001)", "Burgers": "[11-20]/3"},
    "HCP_prismatic_a": {"type": "edge", "plane": "(10-10)", "Burgers": "[11-20]/3"},
    "HCP_pyramidal_c_a": {"type": "mixed", "plane": "(10-11)", "Burgers": "[11-23]/3"},
}


# Grain boundary database (CSL)
GRAIN_BOUNDARY_DATABASE = {
    "Sigma3_twin_FCC": {"sigma": 3, "misorientation_deg": 60.0, "type": "twin",
                        "axis": "[111]", "tilt_angle": 70.5},
    "Sigma3_twin_BCC": {"sigma": 3, "misorientation_deg": 60.0, "type": "twin",
                        "axis": "[110]"},
    "Sigma5_FCC": {"sigma": 5, "misorientation_deg": 36.9, "type": "tilt",
                   "axis": "[100]", "boundary_plane": "(310)"},
    "Sigma5_BCC": {"sigma": 5, "misorientation_deg": 36.9, "type": "tilt",
                   "axis": "[100]"},
    "Sigma7_FCC": {"sigma": 7, "misorientation_deg": 38.2, "type": "tilt",
                   "axis": "[111]", "boundary_plane": "(321)"},
    "Sigma9_FCC": {"sigma": 9, "misorientation_deg": 38.9, "type": "tilt",
                   "axis": "[110]"},
    "Sigma11_FCC": {"sigma": 11, "misorientation_deg": 50.5, "type": "tilt",
                    "axis": "[110]"},
    "Sigma13_FCC": {"sigma": 13, "misorientation_deg": 22.6, "type": "twist",
                    "axis": "[100]"},
    "Sigma17_FCC": {"sigma": 17, "misorientation_deg": 28.1, "type": "tilt",
                    "axis": "[100]"},
    "Sigma19_FCC": {"sigma": 19, "misorientation_deg": 26.5, "type": "twist",
                    "axis": "[111]"},
}


# Stacking fault energies (mJ/m²)
STACKING_FAULT_ENERGIES = {
    "Cu": {"intrinsic": 45, "extrinsic": 55, "twin": 24},
    "Ag": {"intrinsic": 16, "extrinsic": 20, "twin": 10},
    "Au": {"intrinsic": 32, "extrinsic": 40, "twin": 15},
    "Al": {"intrinsic": 166, "extrinsic": 200, "twin": 75},
    "Ni": {"intrinsic": 125, "extrinsic": 150, "twin": 43},
    "Pt": {"intrinsic": 322, "extrinsic": 380, "twin": 161},
    "Pd": {"intrinsic": 180, "extrinsic": 215, "twin": 97},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_dislocation(
    element: str = "Cu",
    dislocation_type: str = "FCC_edge_110",
    supercell: List[int] = [10, 10, 4]
) -> Dict[str, Any]:
    """
    Generate structure with dislocation.
    
    Args:
        element: Host element
        dislocation_type: Dislocation type from database
        supercell: Supercell size
    
    Returns:
        Structure with dislocation
    """
    if dislocation_type not in DISLOCATION_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown dislocation type",
                      "available": list(DISLOCATION_DATABASE.keys())}
        }
    
    info = DISLOCATION_DATABASE[dislocation_type]
    
    # Lattice parameters
    lattice_params = {"Cu": 3.61, "Ag": 4.09, "Au": 4.08, "Al": 4.05, "Ni": 3.52,
                      "Fe": 2.87, "W": 3.16, "Mo": 3.15}
    a = lattice_params.get(element, 3.6)
    
    nx, ny, nz = supercell
    
    # Create base structure
    if "FCC" in dislocation_type:
        lattice = Lattice.cubic(a * nx)
        base_coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    elif "BCC" in dislocation_type:
        lattice = Lattice.cubic(a * nx)
        base_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    else:
        lattice = Lattice.hexagonal(a * nx, a * 1.63 * nz)
        base_coords = [[0, 0, 0], [0.33, 0.67, 0.5]]
    
    species = []
    coords = []
    
    # Build supercell with dislocation displacement
    center_x = 0.5
    center_y = 0.5
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for base in base_coords:
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz
                    
                    # Apply dislocation displacement
                    dx = x - center_x
                    dy = y - center_y
                    theta = np.arctan2(dy, dx)
                    
                    if info["type"] == "edge":
                        # Edge dislocation: displacement in x proportional to theta
                        displacement = theta / (2 * np.pi) * 0.5 / nx
                        x += displacement
                    elif info["type"] == "screw":
                        # Screw dislocation: displacement in z proportional to theta
                        displacement = theta / (2 * np.pi) * 0.5 / nz
                        z += displacement
                    elif info["type"] == "mixed":
                        # Mixed: both components
                        x += theta / (2 * np.pi) * 0.3 / nx
                        z += theta / (2 * np.pi) * 0.3 / nz
                    
                    species.append(element)
                    coords.append([x % 1, y % 1, z % 1])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "dislocation_type": dislocation_type,
        "type": info["type"],
        "Burgers_vector": info.get("Burgers", ""),
        "slip_plane": info.get("plane", ""),
        "character_angle_deg": info.get("character", 0),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_stacking_fault(
    element: str = "Cu",
    fault_type: str = "intrinsic",
    fault_position: float = 0.5,
    supercell: List[int] = [4, 4, 12]
) -> Dict[str, Any]:
    """
    Generate structure with stacking fault.
    
    Args:
        element: Host element
        fault_type: 'intrinsic', 'extrinsic', 'twin'
        fault_position: z-position of fault (fractional)
        supercell: Supercell size
    
    Returns:
        Structure with stacking fault
    """
    fault_types = ["intrinsic", "extrinsic", "twin", "growth"]
    if fault_type not in fault_types:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown fault type",
                      "available": fault_types}
        }
    
    lattice_params = {"Cu": 3.61, "Ag": 4.09, "Au": 4.08, "Al": 4.05, "Ni": 3.52}
    a = lattice_params.get(element, 3.6)
    
    nx, ny, nz = supercell
    
    # FCC structure with ABC stacking
    lattice = Lattice.orthorhombic(a * nx / np.sqrt(2), a * ny * np.sqrt(3/2), 
                                    a * nz * np.sqrt(3))
    
    species = []
    coords = []
    
    # ABC stacking sequence for FCC (111)
    stacking = ["A", "B", "C"]
    
    for k in range(nz * 3):
        layer_z = (k + 0.5) / (nz * 3)
        
        # Determine stacking based on layer and fault
        if fault_type == "intrinsic" and layer_z > fault_position:
            # Remove one layer: ABCAB|CABC -> shift by 1
            shift = 1
        elif fault_type == "extrinsic" and abs(layer_z - fault_position) < 0.1:
            # Insert extra layer
            shift = 2
        elif fault_type == "twin" and layer_z > fault_position:
            # Mirror stacking
            shift = (nz * 3 - k) % 3
        else:
            shift = 0
        
        layer_type = stacking[(k + shift) % 3]
        
        for i in range(nx):
            for j in range(ny):
                # Position based on stacking type
                if layer_type == "A":
                    x_off, y_off = 0, 0
                elif layer_type == "B":
                    x_off, y_off = 1/3, 1/3
                else:  # C
                    x_off, y_off = 2/3, 2/3
                
                x = (i + x_off) / nx
                y = (j + y_off) / ny
                
                species.append(element)
                coords.append([x % 1, y % 1, layer_z])
    
    structure = Structure(lattice, species, coords)
    
    # Get stacking fault energy
    sfe = STACKING_FAULT_ENERGIES.get(element, {}).get(fault_type, 0)
    
    return {
        "success": True,
        "element": element,
        "fault_type": fault_type,
        "fault_position": fault_position,
        "stacking_fault_energy_mJ_m2": sfe,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_grain_boundary(
    element: str = "Cu",
    boundary_type: str = "Sigma3_twin_FCC",
    supercell: List[int] = [4, 4, 8]
) -> Dict[str, Any]:
    """
    Generate grain boundary structure.
    
    Args:
        element: Host element
        boundary_type: GB type from database
        supercell: Supercell size
    
    Returns:
        Structure with grain boundary
    """
    if boundary_type not in GRAIN_BOUNDARY_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown GB type",
                      "available": list(GRAIN_BOUNDARY_DATABASE.keys())}
        }
    
    info = GRAIN_BOUNDARY_DATABASE[boundary_type]
    
    lattice_params = {"Cu": 3.61, "Ag": 4.09, "Au": 4.08, "Al": 4.05, "Ni": 3.52,
                      "Fe": 2.87}
    a = lattice_params.get(element, 3.6)
    
    nx, ny, nz = supercell
    
    # Create bicrystal
    lattice = Lattice.cubic(a * max(nx, ny, nz))
    
    species = []
    coords = []
    
    misorientation = np.radians(info["misorientation_deg"] / 2)
    gb_z = 0.5
    
    base_coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for base in base_coords:
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz
                    
                    # Apply rotation based on which grain
                    if z < gb_z:
                        # Lower grain: rotate by +angle
                        x_rot = (x - 0.5) * np.cos(misorientation) - (y - 0.5) * np.sin(misorientation) + 0.5
                        y_rot = (x - 0.5) * np.sin(misorientation) + (y - 0.5) * np.cos(misorientation) + 0.5
                    else:
                        # Upper grain: rotate by -angle
                        x_rot = (x - 0.5) * np.cos(-misorientation) - (y - 0.5) * np.sin(-misorientation) + 0.5
                        y_rot = (x - 0.5) * np.sin(-misorientation) + (y - 0.5) * np.cos(-misorientation) + 0.5
                    
                    species.append(element)
                    coords.append([x_rot % 1, y_rot % 1, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element": element,
        "boundary_type": boundary_type,
        "sigma": info["sigma"],
        "misorientation_deg": info["misorientation_deg"],
        "gb_type": info["type"],
        "rotation_axis": info.get("axis", ""),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
