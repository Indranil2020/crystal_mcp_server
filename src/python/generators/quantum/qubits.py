"""
quantum/qubits.py - Qubit Material Structures

Comprehensive qubit material generation for quantum computing:
- Superconducting qubits (transmon, flux, charge)
- Spin qubits (Si:P, NV centers, quantum dots)
- Topological qubits (Majorana platforms)
- Trapped ion qubit materials
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Qubit platform database
QUBIT_PLATFORMS = {
    # Superconducting qubits
    "transmon_Al": {
        "type": "superconducting", "material": "Al", "Tc_K": 1.2,
        "junction": "Al/AlOx/Al", "T1_us": 100, "T2_us": 50,
        "description": "Transmon qubit with Al Josephson junction"
    },
    "transmon_Nb": {
        "type": "superconducting", "material": "Nb", "Tc_K": 9.3,
        "junction": "Nb/AlOx/Nb", "T1_us": 200, "T2_us": 100,
        "description": "Transmon qubit with Nb technology"
    },
    "flux_qubit": {
        "type": "superconducting", "material": "Al", 
        "n_junctions": 3, "persistent_current": True,
        "description": "Three-junction flux qubit"
    },
    "fluxonium": {
        "type": "superconducting", "material": "Al",
        "superinductor": True, "T1_us": 1000,
        "description": "Fluxonium with superinductor"
    },
    
    # Spin qubits
    "Si_P_donor": {
        "type": "spin", "host": "Si", "dopant": "P",
        "T2_ms": 0.5, "single_electron": True,
        "description": "Phosphorus donor in Si"
    },
    "Si_Ge_QD": {
        "type": "spin", "host": "SiGe", 
        "T2_ms": 0.02, "electrostatically_defined": True,
        "description": "Si/SiGe quantum dot spin qubit"
    },
    "NV_diamond": {
        "type": "spin", "host": "diamond", "defect": "NV",
        "T2_ms": 2.0, "room_temp": True, "optical_readout": True,
        "description": "Nitrogen-vacancy center in diamond"
    },
    "SiV_diamond": {
        "type": "spin", "host": "diamond", "defect": "SiV",
        "T2_us": 100, "optical_readout": True,
        "description": "Silicon-vacancy center in diamond"
    },
    
    # Topological qubits
    "InSb_Al_nanowire": {
        "type": "topological", "semiconductor": "InSb", "superconductor": "Al",
        "g_factor": 50, "SOC_strong": True, "Majorana": True,
        "description": "InSb/Al nanowire Majorana platform"
    },
    "InAs_Al_nanowire": {
        "type": "topological", "semiconductor": "InAs", "superconductor": "Al",
        "g_factor": 15, "SOC_strong": True, "Majorana": True,
        "description": "InAs/Al nanowire Majorana platform"
    },
    
    # Trapped ion
    "Ca40_ion": {
        "type": "trapped_ion", "ion": "Ca+", "mass": 40,
        "optical_transition": "397nm", "T2_s": 1.0,
        "description": "Ca-40 trapped ion qubit"
    },
    "Yb171_ion": {
        "type": "trapped_ion", "ion": "Yb+", "mass": 171,
        "hyperfine_qubit": True, "T2_s": 10.0,
        "description": "Yb-171 hyperfine qubit"
    },
}


# NV center and defect database
DEFECT_QUBITS = {
    "NV_minus": {
        "defect": "N-V", "charge": -1, "host": "diamond",
        "spin": 1, "ZFS_GHz": 2.87, "optical_ZPL_nm": 637,
        "T1_ms": 6, "T2_ms": 2, "room_temp": True
    },
    "SiV_minus": {
        "defect": "Si-V", "charge": -1, "host": "diamond",
        "spin": 0.5, "ZFS_GHz": 48, "optical_ZPL_nm": 737,
        "T2_us": 100, "mK_required": True
    },
    "GeV_minus": {
        "defect": "Ge-V", "charge": -1, "host": "diamond",
        "spin": 0.5, "ZFS_GHz": 152, "optical_ZPL_nm": 602
    },
    "SnV_minus": {
        "defect": "Sn-V", "charge": -1, "host": "diamond",
        "spin": 0.5, "ZFS_GHz": 850, "optical_ZPL_nm": 619
    },
    "divacancy_SiC": {
        "defect": "VV", "host": "SiC", "polytype": "4H",
        "spin": 1, "ZFS_GHz": 1.3, "optical_ZPL_nm": 1078,
        "room_temp_partially": True
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


def generate_qubit_structure(
    platform: str = "NV_diamond",
    supercell: List[int] = [3, 3, 3]
) -> Dict[str, Any]:
    """
    Generate qubit material structure.
    
    Args:
        platform: Qubit platform from database
        supercell: Supercell dimensions
    
    Returns:
        Qubit structure
    """
    # Check both databases
    if platform in QUBIT_PLATFORMS:
        info = QUBIT_PLATFORMS[platform]
        is_defect = False
    elif platform in DEFECT_QUBITS:
        info = DEFECT_QUBITS[platform]
        is_defect = True
    else:
        all_available = list(QUBIT_PLATFORMS.keys()) + list(DEFECT_QUBITS.keys())
        return {
            "success": False,
            "error": {"code": "INVALID_PLATFORM", "message": f"Unknown platform",
                      "available": all_available}
        }
    
    host = info.get("host", info.get("material", "Si"))
    
    # Create host structure
    lattice_params = {
        "diamond": 3.567, "Si": 5.43, "SiGe": 5.5, "SiC": 3.08,
        "Al": 4.05, "Nb": 3.30, "InSb": 6.48, "InAs": 6.06
    }
    
    a = lattice_params.get(host, 5.0)
    
    species = []
    coords = []
    
    nx, ny, nz = supercell
    
    if host == "diamond":
        lattice = Lattice.cubic(a * nx)
        
        # Diamond structure
        basis = [
            [0, 0, 0], [0.25, 0.25, 0.25],
            [0.5, 0.5, 0], [0.75, 0.75, 0.25],
            [0.5, 0, 0.5], [0.75, 0.25, 0.75],
            [0, 0.5, 0.5], [0.25, 0.75, 0.75]
        ]
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for b in basis:
                        x = (i + b[0]) / nx
                        y = (j + b[1]) / ny
                        z = (k + b[2]) / nz
                        species.append("C")
                        coords.append([x % 1, y % 1, z % 1])
        
        # Add defect at center if defect qubit
        if is_defect:
            center_idx = len(species) // 2
            defect = info.get("defect", "NV")
            
            if "NV" in platform or "NV" in defect:
                # Replace one C with N, remove neighbor for vacancy
                species[center_idx] = "N"
                if center_idx + 1 < len(species):
                    del species[center_idx + 1]
                    del coords[center_idx + 1]
            elif "SiV" in platform or "Si" in defect:
                # Replace two Cs with Si
                species[center_idx] = "Si"
                if center_idx + 1 < len(species):
                    del species[center_idx + 1]
                    del coords[center_idx + 1]
    
    elif host in ["Si", "SiGe"]:
        lattice = Lattice.cubic(a * nx)
        
        basis = [[0, 0, 0], [0.25, 0.25, 0.25]]
        
        for i in range(nx * 2):
            for j in range(ny * 2):
                for k in range(nz * 2):
                    for b in basis:
                        x = (i + b[0]) / (nx * 2)
                        y = (j + b[1]) / (ny * 2)
                        z = (k + b[2]) / (nz * 2)
                        species.append("Si")
                        coords.append([x % 1, y % 1, z % 1])
        
        # Add P dopant if Si:P
        if "P" in platform or info.get("dopant") == "P":
            center_idx = len(species) // 2
            species[center_idx] = "P"
    
    elif host in ["InSb", "InAs"]:
        lattice = Lattice.cubic(a)
        species = ["In", "In", "In", "In", "Sb" if "Sb" in host else "As"] * nx
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.25, 0.25, 0.25]]
        coords = coords * nx
    
    else:
        lattice = Lattice.cubic(a * nx)
        species = [host[:2]] * (nx * ny * nz * 4)
        coords = [[i/nx, j/ny, k/nz] for i in range(nx) for j in range(ny) for k in range(nz) for _ in range(4)]
    
    if species and coords:
        structure = Structure(lattice, species[:len(coords)], coords[:len(species)])
    else:
        structure = Structure(lattice, ["C"], [[0, 0, 0]])
    
    result = {
        "success": True,
        "platform": platform,
        "qubit_type": info.get("type", "unknown"),
        "host_material": host,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
    
    # Add platform-specific properties
    if "T1" in str(info):
        result["T1"] = info.get("T1_us", info.get("T1_ms", 0) * 1000)
    if "T2" in str(info):
        result["T2"] = info.get("T2_us", info.get("T2_ms", 0) * 1000)
    if is_defect:
        result["ZFS_GHz"] = info.get("ZFS_GHz", 0)
        result["optical_ZPL_nm"] = info.get("optical_ZPL_nm", 0)
        result["room_temp_operation"] = info.get("room_temp", False)
    if info.get("Majorana"):
        result["supports_Majorana"] = True
    
    return result


def generate_nv_center(
    supercell: List[int] = [4, 4, 4],
    n_nv_centers: int = 1
) -> Dict[str, Any]:
    """
    Generate diamond with NV center(s).
    
    Args:
        supercell: Supercell size
        n_nv_centers: Number of NV centers
    
    Returns:
        Diamond structure with NV centers
    """
    result = generate_qubit_structure("NV_minus", supercell)
    
    if result["success"]:
        result["n_nv_centers"] = n_nv_centers
        result["qubit_type"] = "spin"
        result["room_temp_operation"] = True
        result["optical_addressable"] = True
    
    return result


def generate_spin_qubit_array(
    platform: str = "Si_P_donor",
    n_qubits: int = 4,
    spacing_nm: float = 20.0
) -> Dict[str, Any]:
    """
    Generate array of spin qubits.
    
    Args:
        platform: Spin qubit platform
        n_qubits: Number of qubits in array
        spacing_nm: Qubit spacing in nm
    
    Returns:
        Spin qubit array structure
    """
    if platform not in QUBIT_PLATFORMS:
        return {"success": False, "error": {"code": "INVALID_PLATFORM", "message": f"Unknown platform"}}
    
    info = QUBIT_PLATFORMS[platform]
    host = info.get("host", "Si")
    
    a = 5.43 if host == "Si" else 5.5  # Si or SiGe
    spacing_A = spacing_nm * 10
    
    # Calculate supercell size
    n_cells = int(spacing_A / a) * n_qubits + 5
    
    lattice = Lattice.cubic(a * n_cells)
    
    species = []
    coords = []
    
    # Build Si host
    for i in range(n_cells):
        for j in range(3):
            for k in range(3):
                species.append("Si")
                coords.append([i / n_cells, j / 3, k / 3])
    
    # Add P dopants as qubits
    for q in range(n_qubits):
        x_pos = (q + 1) * spacing_A / (a * n_cells)
        species.append("P")
        coords.append([x_pos, 0.5, 0.5])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "platform": platform,
        "n_qubits": n_qubits,
        "qubit_spacing_nm": spacing_nm,
        "array_type": "linear",
        "host_material": host,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
