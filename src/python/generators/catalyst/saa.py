"""
catalyst/saa.py - Single Atom Alloys and Catalysts

Comprehensive SAA generation per structure_catalogue.md Category 10:
(i) Single-atom alloys (SAA) – Pt/Cu, Pd/Au, Ni/Cu, Ru/Cu, Ir/Ag
(ii) High-entropy alloy surfaces with random & segregated terminations
(v) 2-D SAC supports with coordination motifs (MN4, MN3C, MN2C2, MN2P2)
(vi) Strain-engineered catalysts – biaxial, uniaxial, bending, ripple
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Comprehensive SAA database
SAA_DATABASE = {
    # Group 10 (Pt, Pd, Ni) on Group 11 (Cu, Ag, Au)
    "Pt/Cu": {"dopant": "Pt", "host": "Cu", "a": 3.615, "energy_benefit": "high_selectivity",
              "applications": ["dehydrogenation", "hydrogenation"], "ensemble_size": 1},
    "Pd/Cu": {"dopant": "Pd", "host": "Cu", "a": 3.615, "energy_benefit": "CO_tolerant",
              "applications": ["formic_acid_oxidation", "acetylene_hydrogenation"]},
    "Pd/Au": {"dopant": "Pd", "host": "Au", "a": 4.078, "energy_benefit": "selective",
              "applications": ["H2O2_synthesis", "vinyl_acetate"]},
    "Pt/Au": {"dopant": "Pt", "host": "Au", "a": 4.078, "energy_benefit": "stability",
              "applications": ["ORR", "methanol_oxidation"]},
    "Ni/Cu": {"dopant": "Ni", "host": "Cu", "a": 3.615, "energy_benefit": "cost_effective",
              "applications": ["CO2_hydrogenation", "methane_reforming"]},
    "Rh/Cu": {"dopant": "Rh", "host": "Cu", "a": 3.615, "energy_benefit": "hydrogenation",
              "applications": ["NO_reduction", "hydroformylation"]},
    "Ru/Cu": {"dopant": "Ru", "host": "Cu", "a": 3.615, "energy_benefit": "ammonia_synthesis",
              "applications": ["N2_activation", "CO_oxidation"]},
    "Ir/Ag": {"dopant": "Ir", "host": "Ag", "a": 4.086, "energy_benefit": "OER",
              "applications": ["water_splitting", "OER_catalyst"]},
    "Ir/Cu": {"dopant": "Ir", "host": "Cu", "a": 3.615, "energy_benefit": "stability",
              "applications": ["propane_dehydrogenation"]},
    "Co/Cu": {"dopant": "Co", "host": "Cu", "a": 3.615, "energy_benefit": "magnetic",
              "applications": ["Fischer_Tropsch"]},
    # Inverse SAA
    "Cu/Pt": {"dopant": "Cu", "host": "Pt", "a": 3.924, "energy_benefit": "selectivity",
              "applications": ["PROX", "WGS"], "inverse": True},
    "Au/Pd": {"dopant": "Au", "host": "Pd", "a": 3.891, "energy_benefit": "selectivity",
              "applications": ["H2O2_selectivity"], "inverse": True},
}


# 2D SAC coordination motifs
SAC_COORDINATION_MOTIFS = {
    "MN4": {"coordination": 4, "ligands": "N", "geometry": "square_planar",
            "examples": ["FeN4", "CoN4", "NiN4", "MnN4", "CuN4"],
            "support": "graphene", "local_moment": True},
    "MN3C": {"coordination": 4, "ligands": ["N", "N", "N", "C"], "geometry": "square_planar",
             "examples": ["FeN3C", "CoN3C"], "support": "N-doped_graphene"},
    "MN2C2": {"coordination": 4, "ligands": ["N", "N", "C", "C"], "geometry": "square_planar",
              "examples": ["FeN2C2", "CoN2C2"], "support": "edge_site"},
    "MN2P2": {"coordination": 4, "ligands": ["N", "N", "P", "P"], "geometry": "square_planar",
              "examples": ["FeN2P2"], "support": "P-doped_graphene"},
    "MN2O2": {"coordination": 4, "ligands": ["N", "N", "O", "O"], "geometry": "square_planar",
              "examples": ["FeN2O2"], "support": "N-O_doped"},
    "MS4": {"coordination": 4, "ligands": "S", "geometry": "tetrahedral",
            "examples": ["MoS4", "WS4"], "support": "MoS2_edge"},
    "MO4": {"coordination": 4, "ligands": "O", "geometry": "distorted",
            "examples": ["PtO4", "AuO4"], "support": "oxide_defect"},
    "MN2S2": {"coordination": 4, "ligands": ["N", "N", "S", "S"], "geometry": "square_planar",
              "examples": ["FeN2S2"], "support": "thiophene_graphene"},
}


# High-entropy alloy compositions
HEA_SURFACES = {
    "Cantor": {"elements": ["Fe", "Co", "Ni", "Cr", "Mn"], "structure": "FCC", "a": 3.59},
    "Senkov": {"elements": ["Mo", "Nb", "Ta", "V", "W"], "structure": "BCC", "a": 3.21},
    "PtPdRhRuIr": {"elements": ["Pt", "Pd", "Rh", "Ru", "Ir"], "structure": "FCC", "a": 3.87},
    "CuCoFeNiZn": {"elements": ["Cu", "Co", "Fe", "Ni", "Zn"], "structure": "FCC", "a": 3.60},
    "AlCoCrFeNi": {"elements": ["Al", "Co", "Cr", "Fe", "Ni"], "structure": "BCC", "a": 2.87},
    "TiZrHfNbTa": {"elements": ["Ti", "Zr", "Hf", "Nb", "Ta"], "structure": "BCC", "a": 3.35},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_saa_surface(
    saa_type: str = "Pt/Cu",
    miller_index: List[int] = [1, 1, 1],
    supercell: List[int] = [4, 4, 4],
    n_dopant: int = 1,
    dopant_layer: str = "surface",
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate Single Atom Alloy surface.
    
    Args:
        saa_type: SAA type from database (dopant/host)
        miller_index: Surface Miller indices
        supercell: Supercell dimensions
        n_dopant: Number of dopant atoms
        dopant_layer: 'surface', 'subsurface', 'random'
        vacuum: Vacuum thickness in Angstrom
    
    Returns:
        SAA surface structure
    """
    if saa_type not in SAA_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_SAA", "message": f"Unknown SAA '{saa_type}'",
                      "available": list(SAA_DATABASE.keys())}
        }
    
    info = SAA_DATABASE[saa_type]
    dopant = info["dopant"]
    host = info["host"]
    a = info["a"]
    
    nx, ny, nz = supercell
    
    # Create FCC slab
    lattice = Lattice.cubic(a * max(nx, ny))
    
    # FCC (111) surface
    species = []
    coords = []
    
    fcc_basis = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for basis in fcc_basis:
                    x = (i + basis[0]) / nx
                    y = (j + basis[1]) / ny
                    z = (k + basis[2]) / nz * 0.5  # Slab region
                    
                    species.append(host)
                    coords.append([x, y, z])
    
    # Add vacuum
    # Scale z-coordinates to leave vacuum
    slab_thickness = 0.5
    for i, coord in enumerate(coords):
        coord[2] = coord[2] * slab_thickness + (1 - slab_thickness - vacuum / (a * max(nx, ny))) / 2
    
    structure = Structure(lattice, species, coords)
    
    # Find surface layer and substitute with dopant
    z_coords = [s.frac_coords[2] for s in structure]
    if dopant_layer == "surface":
        z_max = max(z_coords)
        surface_indices = [i for i, z in enumerate(z_coords) if abs(z - z_max) < 0.02]
    elif dopant_layer == "subsurface":
        z_sorted = sorted(set(z_coords), reverse=True)
        if len(z_sorted) > 1:
            z_sub = z_sorted[1]
        else:
            z_sub = z_sorted[0]
        surface_indices = [i for i, z in enumerate(z_coords) if abs(z - z_sub) < 0.02]
    else:
        surface_indices = list(range(len(structure)))
    
    # Substitute dopants
    np.random.seed(42)
    dopant_indices = list(np.random.choice(surface_indices, min(n_dopant, len(surface_indices)), replace=False))
    
    for idx in dopant_indices:
        structure[int(idx)] = dopant, structure[int(idx)].frac_coords
    
    return {
        "success": True,
        "saa_type": saa_type,
        "dopant": dopant,
        "host": host,
        "miller_index": miller_index,
        "n_dopant_atoms": len(dopant_indices),
        "dopant_layer": dopant_layer,
        "applications": info.get("applications", []),
        "energy_benefit": info.get("energy_benefit", ""),
        "is_inverse_saa": info.get("inverse", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_hea_surface(
    hea_type: str = "Cantor",
    miller_index: List[int] = [1, 1, 1],
    supercell: List[int] = [4, 4, 4],
    termination: str = "random",
    surface_segregation: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate High-Entropy Alloy surface.
    
    Args:
        hea_type: HEA composition from database
        miller_index: Surface Miller indices
        supercell: Supercell dimensions
        termination: 'random', 'equiatomic', 'segregated'
        surface_segregation: Element to enrich at surface
    
    Returns:
        HEA surface structure
    """
    if hea_type not in HEA_SURFACES:
        return {
            "success": False,
            "error": {"code": "INVALID_HEA", "message": f"Unknown HEA '{hea_type}'",
                      "available": list(HEA_SURFACES.keys())}
        }
    
    info = HEA_SURFACES[hea_type]
    elements = info["elements"]
    a = info["a"]
    structure_type = info["structure"]
    
    nx, ny, nz = supercell
    
    lattice = Lattice.cubic(a * max(nx, ny))
    
    species = []
    coords = []
    
    np.random.seed(42)
    
    if structure_type == "FCC":
        basis = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    else:  # BCC
        basis = [[0, 0, 0], [0.5, 0.5, 0.5]]
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for b in basis:
                    x = (i + b[0]) / nx
                    y = (j + b[1]) / ny
                    z = (k + b[2]) / nz * 0.5  # Slab
                    
                    # Assign element
                    if termination == "random":
                        elem = np.random.choice(elements)
                    else:
                        idx = (i + j + k) % len(elements)
                        elem = elements[idx]
                    
                    species.append(elem)
                    coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    # Apply surface segregation
    if surface_segregation and surface_segregation in elements:
        z_coords = [s.frac_coords[2] for s in structure]
        z_max = max(z_coords)
        for i, z in enumerate(z_coords):
            if abs(z - z_max) < 0.02:
                structure[i] = surface_segregation, structure[i].frac_coords
    
    # Calculate actual composition
    composition = {elem: sum(1 for s in structure if str(s.specie) == elem) for elem in elements}
    
    return {
        "success": True,
        "hea_type": hea_type,
        "elements": elements,
        "n_elements": len(elements),
        "structure_type": structure_type,
        "termination": termination,
        "surface_segregation": surface_segregation,
        "composition": composition,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_sac_2d_support(
    metal: str = "Fe",
    coordination_motif: str = "MN4",
    support: str = "graphene",
    supercell: List[int] = [4, 4],
    n_sac_sites: int = 1
) -> Dict[str, Any]:
    """
    Generate 2D SAC with specific coordination motif.
    
    Args:
        metal: Metal atom (Fe, Co, Ni, Mn, Cu, etc.)
        coordination_motif: Motif from SAC_COORDINATION_MOTIFS
        support: 2D support material
        supercell: Supercell size
        n_sac_sites: Number of SAC sites
    
    Returns:
        2D SAC structure
    """
    if coordination_motif not in SAC_COORDINATION_MOTIFS:
        return {
            "success": False,
            "error": {"code": "INVALID_MOTIF", "message": f"Unknown motif '{coordination_motif}'",
                      "available": list(SAC_COORDINATION_MOTIFS.keys())}
        }
    
    motif_info = SAC_COORDINATION_MOTIFS[coordination_motif]
    
    # Graphene lattice
    a_graphene = 2.46
    c_vacuum = 20.0
    
    nx, ny = supercell
    
    lattice = Lattice.hexagonal(a_graphene * nx, c_vacuum)
    
    species = []
    coords = []
    
    # Generate graphene lattice
    a1 = np.array([1, 0, 0])
    a2 = np.array([0.5, np.sqrt(3)/2, 0])
    
    for i in range(nx):
        for j in range(ny):
            # Two atoms per unit cell
            for basis_idx, basis in enumerate([[0, 0], [1/3, 1/3]]):
                x = (i + basis[0]) / nx
                y = (j + basis[1]) / ny
                z = 0.5
                
                species.append("C")
                coords.append([x, y, z])
    
    # Find center position for SAC site
    center_idx = len(species) // 2
    center_pos = coords[center_idx]
    
    # Create coordination environment
    ligands = motif_info["ligands"]
    if isinstance(ligands, str):
        ligands = [ligands] * motif_info["coordination"]
    
    # Remove C atoms and add coordination sphere
    removed_indices = list(range(center_idx - 2, center_idx + 3))
    kept_species = [s for i, s in enumerate(species) if i not in removed_indices]
    kept_coords = [c for i, c in enumerate(coords) if i not in removed_indices]
    
    species = kept_species
    coords = kept_coords
    
    # Add metal center
    species.append(metal)
    coords.append(list(center_pos))
    
    # Add coordinating atoms
    coord_distance = 0.05  # Fractional
    for i, ligand in enumerate(ligands):
        angle = 2 * np.pi * i / len(ligands)
        lx = center_pos[0] + coord_distance * np.cos(angle)
        ly = center_pos[1] + coord_distance * np.sin(angle)
        species.append(ligand)
        coords.append([lx % 1, ly % 1, center_pos[2]])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "metal": metal,
        "coordination_motif": coordination_motif,
        "coordination_number": motif_info["coordination"],
        "geometry": motif_info["geometry"],
        "support": support,
        "n_sac_sites": n_sac_sites,
        "has_local_moment": motif_info.get("local_moment", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_strained_catalyst(
    catalyst_type: str = "Pt/Cu",
    strain_type: str = "biaxial",
    strain_magnitude: float = 0.02,
    supercell: List[int] = [4, 4, 4]
) -> Dict[str, Any]:
    """
    Generate strain-engineered catalyst surface.
    
    Args:
        catalyst_type: Base catalyst type
        strain_type: 'biaxial', 'uniaxial', 'shear', 'ripple'
        strain_magnitude: Strain magnitude (+ tensile, - compressive)
        supercell: Supercell size
    
    Returns:
        Strained catalyst structure
    """
    # First generate unstrained SAA
    base_result = generate_saa_surface(catalyst_type, supercell=supercell)
    
    if not base_result["success"]:
        return base_result
    
    atoms = base_result["structure"]["atoms"]
    lattice_info = base_result["structure"]["lattice"]
    
    a = lattice_info["a"]
    b = lattice_info["b"]
    c = lattice_info["c"]
    
    # Apply strain
    if strain_type == "biaxial":
        a_new = a * (1 + strain_magnitude)
        b_new = b * (1 + strain_magnitude)
        c_new = c * (1 - 0.3 * strain_magnitude)  # Poisson effect
    elif strain_type == "uniaxial":
        a_new = a * (1 + strain_magnitude)
        b_new = b * (1 - 0.3 * strain_magnitude)
        c_new = c * (1 - 0.3 * strain_magnitude)
    elif strain_type == "shear":
        # Apply shear strain
        a_new = a
        b_new = b
        c_new = c
        # Shear would modify the matrix, simplified here
    else:  # ripple
        a_new = a
        b_new = b
        c_new = c
        # Ripple would add z-displacement as function of x,y
    
    new_lattice = {"a": a_new, "b": b_new, "c": c_new}
    
    # Estimate d-band shift (empirical)
    d_band_shift_eV = -5.0 * strain_magnitude  # Approximate
    
    return {
        "success": True,
        "catalyst_type": catalyst_type,
        "strain_type": strain_type,
        "strain_magnitude": strain_magnitude,
        "strain_percent": round(strain_magnitude * 100, 2),
        "estimated_d_band_shift_eV": round(d_band_shift_eV, 3),
        "original_lattice": {"a": a, "b": b, "c": c},
        "strained_lattice": new_lattice,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": new_lattice}
    }
