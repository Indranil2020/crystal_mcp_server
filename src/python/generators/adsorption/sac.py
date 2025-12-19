"""
adsorption/sac.py - Single Atom Catalysts on Supports

Comprehensive SAC generation per structure_catalogue.md Category 6 & 10:
(ii) Single-atom catalysts on various supports (graphene, MN4, oxides)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# SAC metal database
SAC_METALS = {
    "Fe": {"atomic_number": 26, "common_oxidation": [2, 3], "spin": "high", "applications": ["ORR", "NRR"]},
    "Co": {"atomic_number": 27, "common_oxidation": [2, 3], "spin": "low", "applications": ["ORR", "HER"]},
    "Ni": {"atomic_number": 28, "common_oxidation": [2], "spin": "low", "applications": ["CO2RR", "HER"]},
    "Cu": {"atomic_number": 29, "common_oxidation": [1, 2], "spin": "none", "applications": ["CO2RR"]},
    "Mn": {"atomic_number": 25, "common_oxidation": [2, 3, 4], "spin": "high", "applications": ["OER", "ORR"]},
    "Cr": {"atomic_number": 24, "common_oxidation": [3], "spin": "high", "applications": ["NRR"]},
    "V": {"atomic_number": 23, "common_oxidation": [3, 4, 5], "spin": "variable", "applications": ["NRR"]},
    "Zn": {"atomic_number": 30, "common_oxidation": [2], "spin": "none", "applications": ["CO2RR"]},
    "Pt": {"atomic_number": 78, "common_oxidation": [2, 4], "spin": "none", "applications": ["HER", "ORR"]},
    "Pd": {"atomic_number": 46, "common_oxidation": [2], "spin": "none", "applications": ["hydrogenation"]},
    "Rh": {"atomic_number": 45, "common_oxidation": [3], "spin": "low", "applications": ["hydroformylation"]},
    "Ir": {"atomic_number": 77, "common_oxidation": [3, 4], "spin": "low", "applications": ["OER"]},
    "Ru": {"atomic_number": 44, "common_oxidation": [2, 3], "spin": "variable", "applications": ["NRR", "OER"]},
    "Mo": {"atomic_number": 42, "common_oxidation": [4, 6], "spin": "none", "applications": ["HER"]},
    "W": {"atomic_number": 74, "common_oxidation": [4, 6], "spin": "none", "applications": ["HER"]},
}


# Support database
SAC_SUPPORTS = {
    # Carbon-based
    "graphene": {"type": "2D", "a": 2.46, "defect_sites": ["vacancy", "edge", "dopant"]},
    "N_graphene": {"type": "2D", "a": 2.46, "N_content": 0.1, "defect_sites": ["pyridinic", "graphitic"]},
    "graphene_oxide": {"type": "2D", "a": 2.46, "O_content": 0.3, "defect_sites": ["epoxy", "hydroxyl"]},
    "carbon_nitride": {"type": "2D", "a": 7.14, "formula": "C3N4", "defect_sites": ["N_cavity"]},
    
    # Oxide supports
    "CeO2_111": {"type": "oxide", "a": 5.41, "reducible": True, "surface": "(111)"},
    "TiO2_101": {"type": "oxide", "a": 3.78, "reducible": True, "surface": "(101)"},
    "Al2O3_0001": {"type": "oxide", "a": 4.76, "reducible": False, "surface": "(0001)"},
    "MgO_100": {"type": "oxide", "a": 4.21, "reducible": False, "surface": "(100)"},
    "ZnO_0001": {"type": "oxide", "a": 3.25, "reducible": True, "surface": "(0001)"},
    
    # 2D materials
    "MoS2": {"type": "2D_TMD", "a": 3.16, "defect_sites": ["S_vacancy", "Mo_edge"]},
    "BN": {"type": "2D", "a": 2.50, "defect_sites": ["B_vacancy", "N_vacancy"]},
}


# Coordination motifs
COORDINATION_MOTIFS = {
    "MN4": {"n_ligands": 4, "ligand": "N", "geometry": "square_planar", "spin_state": "variable"},
    "MN3C": {"n_ligands": 4, "ligand": ["N", "N", "N", "C"], "geometry": "square_planar"},
    "MN2C2": {"n_ligands": 4, "ligand": ["N", "N", "C", "C"], "geometry": "square_planar"},
    "MO4": {"n_ligands": 4, "ligand": "O", "geometry": "square_planar"},
    "MN2O2": {"n_ligands": 4, "ligand": ["N", "N", "O", "O"], "geometry": "square_planar"},
    "MS4": {"n_ligands": 4, "ligand": "S", "geometry": "tetrahedral"},
    "MN3": {"n_ligands": 3, "ligand": "N", "geometry": "trigonal"},
    "MC3": {"n_ligands": 3, "ligand": "C", "geometry": "trigonal"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_sac_graphene(
    metal: str = "Fe",
    coordination: str = "MN4",
    supercell: List[int] = [5, 5],
    n_sac_sites: int = 1
) -> Dict[str, Any]:
    """
    Generate SAC on graphene-based support.
    
    Args:
        metal: Metal atom
        coordination: Coordination motif
        supercell: Supercell size
        n_sac_sites: Number of SAC sites
    
    Returns:
        SAC structure on graphene
    """
    if metal not in SAC_METALS:
        return {
            "success": False,
            "error": {"code": "INVALID_METAL", "message": f"Unknown metal '{metal}'",
                      "available": list(SAC_METALS.keys())}
        }
    
    if coordination not in COORDINATION_MOTIFS:
        return {
            "success": False,
            "error": {"code": "INVALID_COORDINATION", "message": f"Unknown coordination",
                      "available": list(COORDINATION_MOTIFS.keys())}
        }
    
    metal_info = SAC_METALS[metal]
    coord_info = COORDINATION_MOTIFS[coordination]
    
    a = 2.46  # Graphene lattice constant
    c = 20.0  # Vacuum
    
    nx, ny = supercell
    
    lattice = Lattice.hexagonal(a * nx, c)
    
    species = []
    coords = []
    
    # Generate graphene lattice
    for i in range(nx):
        for j in range(ny):
            # Two atoms per unit cell
            for basis_idx, basis in enumerate([[0, 0], [1/3, 1/3]]):
                x = (i + basis[0]) / nx
                y = (j + basis[1]) / ny
                z = 0.5
                
                species.append("C")
                coords.append([x, y, z])
    
    # Create SAC site at center
    center_x = 0.5
    center_y = 0.5
    
    # Remove carbon atoms near SAC site
    ligands = coord_info["ligand"]
    if isinstance(ligands, str):
        ligands = [ligands] * coord_info["n_ligands"]
    
    # Filter out carbons near center (create vacancy)
    vacancy_radius = 0.08
    kept_species = []
    kept_coords = []
    
    for sp, co in zip(species, coords):
        dist = np.sqrt((co[0] - center_x)**2 + (co[1] - center_y)**2)
        if dist > vacancy_radius:
            kept_species.append(sp)
            kept_coords.append(co)
    
    species = kept_species
    coords = kept_coords
    
    # Add metal center
    species.append(metal)
    coords.append([center_x, center_y, 0.5])
    
    # Add coordinating atoms
    n_ligands = coord_info["n_ligands"]
    ligand_distance = 0.04  # Fractional
    
    for i, lig in enumerate(ligands):
        angle = 2 * np.pi * i / n_ligands
        x = center_x + ligand_distance * np.cos(angle)
        y = center_y + ligand_distance * np.sin(angle)
        species.append(lig)
        coords.append([x, y, 0.5])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "coordination": coordination,
        "geometry": coord_info["geometry"],
        "support": "graphene",
        "n_ligands": n_ligands,
        "ligand_types": list(set(ligands)),
        "applications": metal_info["applications"],
        "spin_state": metal_info["spin"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_sac_oxide(
    metal: str = "Pt",
    support: str = "CeO2_111",
    supercell: List[int] = [3, 3, 3],
    anchoring_site: str = "O_vacancy"
) -> Dict[str, Any]:
    """
    Generate SAC on oxide support.
    
    Args:
        metal: Metal atom
        support: Oxide support type
        supercell: Supercell size
        anchoring_site: 'O_vacancy', 'surface_O', 'step_edge'
    
    Returns:
        SAC on oxide structure
    """
    if metal not in SAC_METALS:
        return {"success": False, "error": {"code": "INVALID_METAL", "message": f"Unknown metal"}}
    
    if support not in SAC_SUPPORTS:
        return {"success": False, "error": {"code": "INVALID_SUPPORT", "message": f"Unknown support",
                "available": [s for s in SAC_SUPPORTS if SAC_SUPPORTS[s]["type"] == "oxide"]}}
    
    support_info = SAC_SUPPORTS[support]
    metal_info = SAC_METALS[metal]
    
    a = support_info["a"]
    nx, ny, nz = supercell
    
    lattice = Lattice.cubic(a * max(nx, ny))
    
    species = []
    coords = []
    
    # Generate oxide slab
    if "CeO2" in support:
        # Fluorite structure
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = i / nx
                    y = j / ny
                    z = k / nz * 0.4  # Slab
                    
                    species.append("Ce")
                    coords.append([x, y, z])
                    
                    # O atoms
                    for ox, oy, oz in [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25]]:
                        species.append("O")
                        coords.append([(x + ox/nx) % 1, (y + oy/ny) % 1, z + oz/nz * 0.4])
    else:
        # Generic oxide
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = i / nx
                    y = j / ny
                    z = k / nz * 0.4
                    
                    metal_elem = support.split("O")[0].split("_")[0].replace("2", "").replace("3", "")
                    species.append(metal_elem if metal_elem else "M")
                    coords.append([x, y, z])
                    species.append("O")
                    coords.append([(x + 0.25/nx) % 1, (y + 0.25/ny) % 1, z])
    
    # Find surface and add SAC
    z_max = max(c[2] for c in coords)
    
    # Add metal at surface
    species.append(metal)
    coords.append([0.5, 0.5, z_max + 0.05])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "support": support,
        "anchoring_site": anchoring_site,
        "is_reducible": support_info.get("reducible", False),
        "applications": metal_info["applications"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_sac_tmd(
    metal: str = "Co",
    support: str = "MoS2",
    defect_type: str = "S_vacancy",
    supercell: List[int] = [4, 4]
) -> Dict[str, Any]:
    """
    Generate SAC on 2D TMD support.
    
    Args:
        metal: Metal atom
        support: TMD support (MoS2, WS2, etc.)
        defect_type: Type of defect site
        supercell: Supercell size
    
    Returns:
        SAC on TMD structure
    """
    if metal not in SAC_METALS:
        return {"success": False, "error": {"code": "INVALID_METAL", "message": f"Unknown metal"}}
    
    metal_info = SAC_METALS[metal]
    
    # MoS2 parameters
    a = 3.16
    c = 20.0
    
    nx, ny = supercell
    
    lattice = Lattice.hexagonal(a * nx, c)
    
    species = []
    coords = []
    
    # Generate MoS2 monolayer
    for i in range(nx):
        for j in range(ny):
            x = (i + 1/3) / nx
            y = (j + 1/3) / ny
            
            # Mo
            species.append("Mo")
            coords.append([x, y, 0.5])
            
            # S above and below
            species.append("S")
            coords.append([x, y, 0.5 + 0.08])
            species.append("S")
            coords.append([x, y, 0.5 - 0.08])
    
    # Create defect at center and add SAC
    center_idx = len(species) // 2
    if defect_type == "S_vacancy":
        # Remove one S atom near center
        del species[center_idx + 1]
        del coords[center_idx + 1]
    
    # Add metal at defect site
    species.append(metal)
    coords.append([0.5, 0.5, 0.5])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "support": support,
        "defect_type": defect_type,
        "applications": metal_info["applications"],
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
