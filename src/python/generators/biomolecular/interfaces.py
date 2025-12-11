"""
biomolecular/interfaces.py - Protein-Surface Interfaces

Generates biomolecule-surface interfaces:
- Protein on Ti, silica, gold
- DNA origami on 2D materials
- Peptide-semiconductor
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Surface protein database
SURFACE_PROTEINS = {
    "fibronectin": {
        "size_A": [100, 50, 30], "domains": ["RGD", "heparin_binding"],
        "surfaces": ["Ti", "gold", "silica"]
    },
    "lysozyme": {
        "size_A": [45, 30, 30], "molecular_weight": 14300,
        "surfaces": ["silica", "mica", "graphene"]
    },
    "IgG": {
        "size_A": [140, 80, 40], "molecular_weight": 150000,
        "surfaces": ["gold", "silicon", "polymer"]
    },
    "BSA": {
        "size_A": [80, 80, 30], "molecular_weight": 66500,
        "surfaces": ["silica", "Ti", "polymer"]
    },
    "collagen": {
        "size_A": [300, 15, 15], "triple_helix": True,
        "surfaces": ["hydroxyapatite", "Ti"]
    },
    "RGD_peptide": {
        "size_A": [15, 10, 10], "sequence": "RGD",
        "surfaces": ["Ti", "gold", "silica"]
    },
}


def generate_protein_surface(
    protein: str = "fibronectin",
    surface: str = "Ti",
    surface_size: List[int] = [100, 100],
    orientation: str = "flat",
    coverage: float = 0.3
) -> Dict[str, Any]:
    """
    Generate protein-surface interface.
    
    Args:
        protein: Protein type
        surface: Surface material
        surface_size: Surface dimensions in Angstrom
        orientation: 'flat', 'upright', 'random'
        coverage: Surface coverage fraction
    
    Returns:
        Protein-surface structure
    """
    if protein not in SURFACE_PROTEINS:
        return {
            "success": False,
            "error": {"code": "INVALID_PROTEIN", "message": f"Unknown protein '{protein}'",
                      "available": list(SURFACE_PROTEINS.keys())}
        }
    
    np.random.seed(42)
    info = SURFACE_PROTEINS[protein]
    protein_size = info["size_A"]
    
    lx, ly = surface_size
    lz = max(protein_size) + 50
    
    atoms = []
    
    # Generate surface
    surface_params = {
        "Ti": {"a": 2.95, "element": "Ti"},
        "gold": {"a": 2.88, "element": "Au"},
        "silica": {"a": 4.91, "elements": ["Si", "O"]},
        "graphene": {"a": 2.46, "element": "C"},
    }
    
    if surface in surface_params:
        params = surface_params[surface]
        a = params["a"]
        elem = params.get("element", params.get("elements", ["X"])[0])
        
        nx = int(lx / a)
        ny = int(ly / a)
        
        for i in range(nx):
            for j in range(ny):
                x = i * a
                y = j * a
                atoms.append({"element": elem, "cartesian": [x, y, 0], "layer": "surface"})
                
                if surface == "silica":
                    atoms.append({"element": "O", "cartesian": [x + a/2, y + a/2, 1.6], "layer": "surface"})
    
    # Add proteins
    n_proteins = max(1, int(coverage * lx * ly / (protein_size[0] * protein_size[1])))
    
    for p in range(n_proteins):
        # Random position
        px = np.random.uniform(protein_size[0] / 2, lx - protein_size[0] / 2)
        py = np.random.uniform(protein_size[1] / 2, ly - protein_size[1] / 2)
        
        if orientation == "flat":
            pz = 5  # Distance from surface
            width, depth, height = protein_size
        elif orientation == "upright":
            pz = 5
            height, depth, width = protein_size
        else:
            pz = 5
            # Random rotation
            width, depth, height = np.random.permutation(protein_size)
        
        # Simplified protein as ellipsoid of backbone atoms
        n_atoms_protein = int(info.get("molecular_weight", 10000) / 110)  # ~110 Da per residue
        
        for i in range(min(n_atoms_protein, 100)):
            # Random position within ellipsoid
            theta = np.random.uniform(0, 2 * np.pi)
            phi = np.random.uniform(0, np.pi)
            r = np.random.uniform(0, 1) ** (1/3)  # Uniform in volume
            
            x = px + r * width / 2 * np.sin(phi) * np.cos(theta)
            y = py + r * depth / 2 * np.sin(phi) * np.sin(theta)
            z = pz + height / 2 + r * height / 2 * np.cos(phi)
            
            # Alternate backbone atoms
            elem = ["N", "C", "C", "O"][i % 4]
            atoms.append({"element": elem, "cartesian": [x, y, z], "molecule": protein, "protein_id": p})
    
    return {
        "success": True,
        "protein": protein,
        "surface": surface,
        "orientation": orientation,
        "n_proteins": n_proteins,
        "surface_coverage": coverage,
        "surface_size_A": surface_size,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_dna_surface(
    surface: str = "graphene",
    dna_type: str = "origami_plate",
    surface_size: List[int] = [200, 200]
) -> Dict[str, Any]:
    """
    Generate DNA origami on 2D surface.
    
    Args:
        surface: Surface material (graphene, MoS2)
        dna_type: DNA structure type
        surface_size: Surface size
    
    Returns:
        DNA-surface structure
    """
    dna_types = {
        "origami_plate": {"shape": "rectangular", "size": [70, 100, 2]},
        "origami_rod": {"shape": "cylindrical", "length": 200, "diameter": 20},
        "origami_barrel": {"shape": "barrel", "length": 100, "diameter": 30},
    }
    
    if dna_type not in dna_types:
        return {"success": False, "error": {"code": "INVALID_TYPE", "message": f"Unknown DNA type"}}
    
    np.random.seed(42)
    info = dna_types[dna_type]
    
    lx, ly = surface_size
    lz = 50
    
    atoms = []
    
    # Generate surface
    if surface == "graphene":
        a = 2.46
        for i in range(int(lx / a)):
            for j in range(int(ly / a)):
                x = i * a
                y = j * a
                atoms.append({"element": "C", "cartesian": [x, y, 0], "layer": "surface"})
    elif surface == "MoS2":
        a = 3.16
        for i in range(int(lx / a)):
            for j in range(int(ly / a)):
                x = i * a
                y = j * a
                atoms.append({"element": "Mo", "cartesian": [x, y, 0], "layer": "surface"})
                atoms.append({"element": "S", "cartesian": [x, y, 1.6], "layer": "surface"})
                atoms.append({"element": "S", "cartesian": [x, y, -1.6], "layer": "surface"})
    
    # Add DNA origami
    center_x = lx / 2
    center_y = ly / 2
    
    if info["shape"] == "rectangular":
        width, length, height = info["size"]
        # Generate grid of dsDNA
        for i in range(int(width / 2)):
            for j in range(int(length / 3.4)):
                x = center_x - width / 2 + i * 2
                y = center_y - length / 2 + j * 3.4
                z = 5  # Distance from surface
                
                # DNA backbone
                atoms.append({"element": "P", "cartesian": [x, y, z], "molecule": "DNA"})
                atoms.append({"element": "O", "cartesian": [x + 0.5, y, z], "molecule": "DNA"})
                atoms.append({"element": "C", "cartesian": [x, y + 0.5, z], "molecule": "DNA"})
    
    return {
        "success": True,
        "surface": surface,
        "dna_type": dna_type,
        "shape": info["shape"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_peptide_surface(
    peptide: str = "RGD",
    surface: str = "TiO2",
    surface_size: List[int] = [50, 50],
    coverage_nmol_nm2: float = 1.0
) -> Dict[str, Any]:
    """
    Generate peptide-semiconductor interface.
    
    Args:
        peptide: Peptide sequence
        surface: Surface material
        surface_size: Surface size in Angstrom
        coverage_nmol_nm2: Coverage in molecules per nm²
    
    Returns:
        Peptide-surface structure
    """
    # Common peptide sequences
    peptide_sequences = {
        "RGD": ["Arg", "Gly", "Asp"],
        "RGDS": ["Arg", "Gly", "Asp", "Ser"],
        "diphenylalanine": ["Phe", "Phe"],
        "KLVFF": ["Lys", "Leu", "Val", "Phe", "Phe"],
    }
    
    if peptide not in peptide_sequences:
        return {"success": False, "error": {"code": "INVALID_PEPTIDE", "message": f"Unknown peptide"}}
    
    np.random.seed(42)
    sequence = peptide_sequences[peptide]
    
    lx, ly = surface_size
    lz = 30
    
    atoms = []
    
    # Generate surface
    if surface == "TiO2":
        a = 4.59
        for i in range(int(lx / a)):
            for j in range(int(ly / a)):
                x = i * a
                y = j * a
                atoms.append({"element": "Ti", "cartesian": [x, y, 0], "layer": "surface"})
                atoms.append({"element": "O", "cartesian": [x + a/2, y, 0.9], "layer": "surface"})
                atoms.append({"element": "O", "cartesian": [x, y + a/2, 0.9], "layer": "surface"})
    elif surface == "ZnO":
        a = 3.25
        for i in range(int(lx / a)):
            for j in range(int(ly / a)):
                x = i * a
                y = j * a
                atoms.append({"element": "Zn", "cartesian": [x, y, 0], "layer": "surface"})
                atoms.append({"element": "O", "cartesian": [x, y, 1.9], "layer": "surface"})
    
    # Calculate number of peptides
    area_nm2 = (lx * ly) / 100  # Convert Å² to nm²
    n_peptides = int(coverage_nmol_nm2 * area_nm2)
    
    for p in range(n_peptides):
        px = np.random.uniform(5, lx - 5)
        py = np.random.uniform(5, ly - 5)
        pz = 5  # Distance from surface
        
        # Add peptide backbone
        for i, res in enumerate(sequence):
            z = pz + i * 3.8  # Residue spacing
            
            # Backbone atoms
            atoms.append({"element": "N", "cartesian": [px, py, z], "molecule": peptide, "peptide_id": p})
            atoms.append({"element": "C", "cartesian": [px + 0.5, py + 0.5, z + 1.3], "molecule": peptide})
            atoms.append({"element": "C", "cartesian": [px, py, z + 2.3], "molecule": peptide})
            atoms.append({"element": "O", "cartesian": [px + 1.2, py, z + 2.8], "molecule": peptide})
    
    return {
        "success": True,
        "peptide": peptide,
        "sequence": sequence,
        "surface": surface,
        "n_peptides": n_peptides,
        "coverage_nmol_nm2": coverage_nmol_nm2,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
