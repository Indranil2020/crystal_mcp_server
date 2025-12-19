"""
catalyst/clusters.py - Supported Metal Clusters

Comprehensive supported cluster generation per structure_catalogue.md Category 10:
(iii) Oxide-supported clusters (Pt4, Ir4, Rh6, Au20) on CeO2, TiO2, ZnO, Al2O3, FeOx
(iv) Zeolite-confined clusters (Cu-CHA, Fe-MFI, Co-AFX)
(vii) Electric-field & gate-tunable catalysis configurations
(viii) Photocatalyst heterojunctions (type-I, type-II, Z-scheme, S-scheme)
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Metal cluster database
CLUSTER_DATABASE = {
    # Pt clusters
    "Pt4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.5, "edge_length": 2.77},
    "Pt6": {"n_atoms": 6, "geometry": "octahedron", "diameter_A": 5.5, "edge_length": 2.77},
    "Pt13": {"n_atoms": 13, "geometry": "icosahedron", "diameter_A": 7.0, "edge_length": 2.77},
    "Pt38": {"n_atoms": 38, "geometry": "truncated_octahedron", "diameter_A": 10.0},
    "Pt55": {"n_atoms": 55, "geometry": "icosahedron", "diameter_A": 11.0},
    # Ir clusters
    "Ir4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.3, "edge_length": 2.71},
    "Ir6": {"n_atoms": 6, "geometry": "octahedron", "diameter_A": 5.2, "edge_length": 2.71},
    # Rh clusters
    "Rh4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.4, "edge_length": 2.69},
    "Rh6": {"n_atoms": 6, "geometry": "octahedron", "diameter_A": 5.3, "edge_length": 2.69},
    "Rh13": {"n_atoms": 13, "geometry": "icosahedron", "diameter_A": 6.8, "edge_length": 2.69},
    # Au clusters
    "Au4": {"n_atoms": 4, "geometry": "rhombus", "diameter_A": 4.8, "planar": True},
    "Au6": {"n_atoms": 6, "geometry": "triangle", "diameter_A": 5.8, "planar": True},
    "Au13": {"n_atoms": 13, "geometry": "planar", "diameter_A": 8.0, "planar": True},
    "Au20": {"n_atoms": 20, "geometry": "tetrahedron", "diameter_A": 9.5},
    "Au55": {"n_atoms": 55, "geometry": "icosahedron", "diameter_A": 14.0},
    # Pd clusters
    "Pd4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.4},
    "Pd6": {"n_atoms": 6, "geometry": "octahedron", "diameter_A": 5.3},
    "Pd13": {"n_atoms": 13, "geometry": "icosahedron", "diameter_A": 7.0},
    # Ru clusters
    "Ru4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.2},
    "Ru6": {"n_atoms": 6, "geometry": "octahedron", "diameter_A": 5.1},
    # Ni clusters
    "Ni4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.0},
    "Ni13": {"n_atoms": 13, "geometry": "icosahedron", "diameter_A": 6.5},
    # Co clusters
    "Co4": {"n_atoms": 4, "geometry": "tetrahedron", "diameter_A": 4.0},
    "Co13": {"n_atoms": 13, "geometry": "icosahedron", "diameter_A": 6.5},
}


# Oxide support database
OXIDE_SUPPORTS = {
    "CeO2": {"a": 5.41, "structure": "fluorite", "reducible": True, 
             "surfaces": ["111", "110", "100"], "oxygen_vacancy_energy": 1.5},
    "TiO2_rutile": {"a": 4.59, "c": 2.96, "structure": "rutile", "reducible": True,
                   "surfaces": ["110", "101", "100"]},
    "TiO2_anatase": {"a": 3.78, "c": 9.51, "structure": "anatase", "reducible": True,
                    "surfaces": ["101", "001", "100"]},
    "ZnO": {"a": 3.25, "c": 5.21, "structure": "wurtzite", "reducible": True,
            "surfaces": ["0001", "10-10"]},
    "Al2O3": {"a": 4.76, "c": 12.99, "structure": "corundum", "reducible": False,
              "surfaces": ["0001", "1-102"]},
    "MgO": {"a": 4.21, "structure": "rocksalt", "reducible": False,
            "surfaces": ["100", "110"]},
    "SiO2": {"a": 4.91, "structure": "quartz", "reducible": False,
             "surfaces": ["0001", "10-10"]},
    "FeOx_magnetite": {"a": 8.40, "structure": "spinel", "reducible": True,
                       "surfaces": ["111", "100", "110"]},
    "FeOx_hematite": {"a": 5.04, "c": 13.75, "structure": "corundum", "reducible": True,
                      "surfaces": ["0001", "10-12"]},
}


# Zeolite frameworks
ZEOLITE_FRAMEWORKS = {
    "CHA": {"name": "Chabazite", "rings": [8, 6, 4], "cage": "t-she", 
            "pore_size_A": 3.8, "Si/Al": 12},
    "MFI": {"name": "ZSM-5", "rings": [10, 6, 5, 4], "channels": "3D",
            "pore_size_A": 5.5, "Si/Al": 25},
    "FAU": {"name": "Faujasite", "rings": [12, 6, 4], "supercage": True,
            "pore_size_A": 7.4, "Si/Al": 2.5},
    "AFX": {"name": "SAPO-56", "rings": [8, 6, 4], "cage": "aft",
            "pore_size_A": 3.6, "Si/Al": 10},
    "BEA": {"name": "Beta", "rings": [12, 6, 5, 4], "channels": "3D",
            "pore_size_A": 6.6, "Si/Al": 20},
    "MOR": {"name": "Mordenite", "rings": [12, 8, 5, 4], "channels": "1D",
            "pore_size_A": 6.5, "Si/Al": 10},
    "LTA": {"name": "Zeolite A", "rings": [8, 6, 4], "cage": "alpha",
            "pore_size_A": 4.2, "Si/Al": 1},
    "FER": {"name": "Ferrierite", "rings": [10, 8, 6, 5], "channels": "2D",
            "pore_size_A": 5.4, "Si/Al": 8},
}


# Photocatalyst heterojunctions
HETEROJUNCTIONS = {
    "type_I": {"description": "Straddling gap", "charge_separation": "poor",
               "examples": [("TiO2", "CdS"), ("ZnO", "CdSe")]},
    "type_II": {"description": "Staggered gap", "charge_separation": "good",
                "examples": [("TiO2", "g-C3N4"), ("WO3", "BiVO4")]},
    "Z_scheme": {"description": "Direct Z-scheme", "charge_separation": "excellent",
                 "examples": [("WO3", "g-C3N4"), ("Fe2O3", "Cu2O")]},
    "S_scheme": {"description": "Step-scheme", "charge_separation": "excellent",
                 "examples": [("TiO2", "g-C3N4"), ("CdS", "WO3")]},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_cluster(
    cluster_type: str = "Pt4",
    return_atoms_only: bool = False
) -> Dict[str, Any]:
    """
    Generate isolated metal cluster.
    
    Args:
        cluster_type: Cluster name from database
        return_atoms_only: Return only atom list
    
    Returns:
        Cluster structure
    """
    if cluster_type not in CLUSTER_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_CLUSTER", "message": f"Unknown cluster '{cluster_type}'",
                      "available": list(CLUSTER_DATABASE.keys())}
        }
    
    info = CLUSTER_DATABASE[cluster_type]
    n_atoms = info["n_atoms"]
    geometry = info["geometry"]
    
    # Extract element
    element = ""
    for char in cluster_type:
        if char.isalpha():
            element += char
        else:
            break
    
    atoms = []
    
    if geometry == "tetrahedron":
        # Regular tetrahedron
        edge = info.get("edge_length", 2.7)
        h = edge * np.sqrt(2/3)
        atoms = [
            [0, 0, 0],
            [edge, 0, 0],
            [edge/2, edge * np.sqrt(3)/2, 0],
            [edge/2, edge * np.sqrt(3)/6, h]
        ]
    elif geometry == "octahedron":
        edge = info.get("edge_length", 2.7)
        d = edge / np.sqrt(2)
        atoms = [
            [d, 0, 0], [-d, 0, 0],
            [0, d, 0], [0, -d, 0],
            [0, 0, d], [0, 0, -d]
        ]
    elif geometry == "icosahedron":
        # Mackay icosahedron
        r = info.get("diameter_A", 7.0) / 2
        phi = (1 + np.sqrt(5)) / 2  # Golden ratio
        
        # 12 vertices + 1 center
        atoms = [[0, 0, 0]]  # Center
        for i in range(12):
            theta = np.arccos(-1 + 2*i/11)
            phi_angle = np.pi * (1 + np.sqrt(5)) * i
            x = r * np.sin(theta) * np.cos(phi_angle)
            y = r * np.sin(theta) * np.sin(phi_angle)
            z = r * np.cos(theta)
            atoms.append([x, y, z])
    elif geometry in ["rhombus", "triangle", "planar"]:
        # Planar Au clusters
        n = info["n_atoms"]
        r = info.get("diameter_A", 5.0) / 2
        for i in range(n):
            angle = 2 * np.pi * i / n
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            atoms.append([x, y, 0])
    else:
        # Generic spherical arrangement
        n = info["n_atoms"]
        r = info.get("diameter_A", 5.0) / 2
        for i in range(n):
            phi = np.arccos(1 - 2*(i+0.5)/n)
            theta = np.pi * (1 + np.sqrt(5)) * i
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)
            atoms.append([x, y, z])
    
    if return_atoms_only:
        return {"atoms": atoms, "element": element}
    
    # Put in large box
    box_size = info.get("diameter_A", 10) + 15
    lattice = Lattice.cubic(box_size)
    
    species = [element] * len(atoms)
    coords = []
    for atom in atoms:
        x = (atom[0] + box_size/2) / box_size
        y = (atom[1] + box_size/2) / box_size
        z = (atom[2] + box_size/2) / box_size
        coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "cluster_type": cluster_type,
        "element": element,
        "n_atoms": len(atoms),
        "geometry": geometry,
        "diameter_A": info.get("diameter_A", 0),
        "is_planar": info.get("planar", False),
        "structure": structure_to_dict(structure)
    }


def generate_supported_cluster(
    cluster: str = "Pt4",
    support: str = "CeO2",
    surface: str = "111",
    supercell: List[int] = [3, 3, 4],
    adsorption_site: str = "O_vacancy"
) -> Dict[str, Any]:
    """
    Generate oxide-supported metal cluster.
    
    Args:
        cluster: Cluster type
        support: Oxide support
        surface: Surface Miller index
        supercell: Supercell size
        adsorption_site: 'atop', 'bridge', 'hollow', 'O_vacancy'
    
    Returns:
        Supported cluster structure
    """
    if cluster not in CLUSTER_DATABASE:
        return {"success": False, "error": {"code": "INVALID_CLUSTER", "message": f"Unknown cluster"}}
    
    if support not in OXIDE_SUPPORTS:
        return {"success": False, "error": {"code": "INVALID_SUPPORT", "message": f"Unknown support",
                "available": list(OXIDE_SUPPORTS.keys())}}
    
    cluster_info = CLUSTER_DATABASE[cluster]
    support_info = OXIDE_SUPPORTS[support]
    
    # Generate cluster atoms
    cluster_result = generate_cluster(cluster, return_atoms_only=True)
    cluster_atoms = cluster_result["atoms"]
    element = cluster_result["element"]
    
    # Generate support slab
    a = support_info["a"]
    nx, ny, nz = supercell
    
    lattice = Lattice.cubic(a * max(nx, ny))
    
    species = []
    coords = []
    
    # Simplified oxide structure
    if support_info["structure"] == "fluorite":
        # CeO2 fluorite
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = (i + 0.0) / nx
                    y = (j + 0.0) / ny
                    z = (k + 0.0) / nz * 0.4  # Slab
                    species.append("Ce")
                    coords.append([x, y, z])
                    
                    # O atoms
                    for ox, oy, oz in [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]:
                        oxi = (i + ox) / nx
                        oyi = (j + oy) / ny
                        ozi = (k + oz) / nz * 0.4
                        species.append("O")
                        coords.append([oxi, oyi, ozi])
    elif "rutile" in support_info["structure"]:
        # TiO2 rutile
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = i / nx
                    y = j / ny
                    z = k / nz * 0.4
                    species.append("Ti")
                    coords.append([x, y, z])
                    species.extend(["O", "O"])
                    coords.extend([[x + 0.3/nx, y + 0.3/ny, z], [x - 0.3/nx, y - 0.3/ny, z]])
    else:
        # Generic oxide
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = i / nx
                    y = j / ny
                    z = k / nz * 0.4
                    species.append(support.split("_")[0].replace("Ox", "").replace("2", "").replace("3", "") or "M")
                    coords.append([x, y, z])
                    species.append("O")
                    coords.append([x + 0.25/nx, y + 0.25/ny, z])
    
    # Find highest z point for cluster placement
    z_max = max(c[2] for c in coords)
    cluster_z_base = z_max + 0.05
    
    # Add cluster atoms
    cluster_diameter = cluster_info.get("diameter_A", 5) / (a * max(nx, ny))
    
    for atom in cluster_atoms[:cluster_info["n_atoms"]]:
        scaled_x = 0.5 + atom[0] / (a * max(nx, ny))
        scaled_y = 0.5 + atom[1] / (a * max(nx, ny))
        scaled_z = cluster_z_base + atom[2] / (a * max(nx, ny))
        
        species.append(element)
        coords.append([scaled_x, scaled_y, scaled_z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "cluster": cluster,
        "support": support,
        "surface": surface,
        "adsorption_site": adsorption_site,
        "cluster_geometry": cluster_info["geometry"],
        "support_reducible": support_info.get("reducible", False),
        "n_cluster_atoms": cluster_info["n_atoms"],
        "n_total_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_zeolite_cluster(
    metal: str = "Cu",
    zeolite: str = "CHA",
    cluster_size: int = 2,
    location: str = "cage"
) -> Dict[str, Any]:
    """
    Generate zeolite-confined metal cluster.
    
    Args:
        metal: Metal element
        zeolite: Zeolite framework type
        cluster_size: Number of metal atoms
        location: 'cage', 'channel', 'intersection'
    
    Returns:
        Zeolite-confined cluster structure
    """
    if zeolite not in ZEOLITE_FRAMEWORKS:
        return {"success": False, "error": {"code": "INVALID_ZEOLITE", 
                "message": f"Unknown zeolite", "available": list(ZEOLITE_FRAMEWORKS.keys())}}
    
    info = ZEOLITE_FRAMEWORKS[zeolite]
    
    # Simplified zeolite structure
    a = 13.0  # Typical CHA
    lattice = Lattice.rhombohedral(a, 94)
    
    species = []
    coords = []
    
    # T-sites (Si/Al)
    n_t_sites = 36  # CHA has 36 T-sites per unit cell
    np.random.seed(42)
    
    for i in range(n_t_sites):
        theta = np.pi * i / n_t_sites
        phi = 2 * np.pi * np.sqrt(5) * i
        r = 0.3 + 0.2 * np.random.random()
        
        x = 0.5 + r * np.sin(theta) * np.cos(phi)
        y = 0.5 + r * np.sin(theta) * np.sin(phi)
        z = 0.5 + r * np.cos(theta)
        
        # Si or Al
        if i % int(info["Si/Al"] + 1) == 0:
            species.append("Al")
        else:
            species.append("Si")
        coords.append([x % 1, y % 1, z % 1])
        
        # Bridging O atoms
        for _ in range(2):
            ox = x + np.random.uniform(-0.05, 0.05)
            oy = y + np.random.uniform(-0.05, 0.05)
            oz = z + np.random.uniform(-0.05, 0.05)
            species.append("O")
            coords.append([ox % 1, oy % 1, oz % 1])
    
    # Add metal cluster at cage center
    if location == "cage":
        center = [0.5, 0.5, 0.5]
    elif location == "channel":
        center = [0.25, 0.25, 0.5]
    else:
        center = [0.5, 0.5, 0.25]
    
    # Metal cluster
    cluster_spacing = 0.03
    for i in range(cluster_size):
        mx = center[0] + cluster_spacing * np.cos(2 * np.pi * i / cluster_size)
        my = center[1] + cluster_spacing * np.sin(2 * np.pi * i / cluster_size)
        mz = center[2]
        species.append(metal)
        coords.append([mx, my, mz])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "cluster_size": cluster_size,
        "zeolite": zeolite,
        "zeolite_name": info["name"],
        "pore_size_A": info["pore_size_A"],
        "location": location,
        "ring_sizes": info["rings"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_photocatalyst_heterojunction(
    junction_type: str = "type_II",
    materials: Optional[Tuple[str, str]] = None,
    interface_orientation: str = "(001)"
) -> Dict[str, Any]:
    """
    Generate photocatalyst heterojunction.
    
    Args:
        junction_type: 'type_I', 'type_II', 'Z_scheme', 'S_scheme'
        materials: Tuple of (material_1, material_2)
        interface_orientation: Interface plane
    
    Returns:
        Heterojunction structure
    """
    if junction_type not in HETEROJUNCTIONS:
        return {"success": False, "error": {"code": "INVALID_JUNCTION",
                "message": f"Unknown junction type", "available": list(HETEROJUNCTIONS.keys())}}
    
    info = HETEROJUNCTIONS[junction_type]
    
    if materials is None:
        materials = info["examples"][0]
    
    mat1, mat2 = materials
    
    # Simplified heterojunction (stacked slabs)
    a = 10.0
    lattice = Lattice.orthorhombic(a, a, 30)
    
    species = []
    coords = []
    
    # Material 1 slab
    for i in range(3):
        for j in range(3):
            for k in range(2):
                x = (i + 0.5) / 3
                y = (j + 0.5) / 3
                z = (k + 0.5) / 6 * 0.4  # Lower slab
                
                if "TiO2" in mat1:
                    species.extend(["Ti", "O", "O"])
                    coords.extend([[x, y, z], [x+0.03, y, z], [x, y+0.03, z]])
                elif "g-C3N4" in mat1:
                    species.extend(["C", "N"])
                    coords.extend([[x, y, z], [x+0.02, y+0.02, z]])
                else:
                    species.append(mat1[0])
                    coords.append([x, y, z])
    
    # Interface region
    interface_z = 0.45
    
    # Material 2 slab
    for i in range(3):
        for j in range(3):
            for k in range(2):
                x = (i + 0.5) / 3
                y = (j + 0.5) / 3
                z = interface_z + (k + 0.5) / 6 * 0.4  # Upper slab
                
                if "TiO2" in mat2:
                    species.extend(["Ti", "O", "O"])
                    coords.extend([[x, y, z], [x+0.03, y, z], [x, y+0.03, z]])
                elif "g-C3N4" in mat2:
                    species.extend(["C", "N"])
                    coords.extend([[x, y, z], [x+0.02, y+0.02, z]])
                else:
                    species.append(mat2[0])
                    coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "junction_type": junction_type,
        "description": info["description"],
        "charge_separation": info["charge_separation"],
        "material_1": mat1,
        "material_2": mat2,
        "interface_orientation": interface_orientation,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
