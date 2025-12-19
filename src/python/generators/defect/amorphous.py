"""
defect/amorphous.py - Amorphous Structure Generation

Comprehensive amorphous structure generation:
- Melt-quench simulation
- Random network models
- Amorphous alloys and glasses
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Amorphous material database
AMORPHOUS_DATABASE = {
    # Oxide glasses
    "SiO2": {"type": "network_glass", "network_former": "Si", 
             "coordination": 4, "density_gcm3": 2.20, "Tg_K": 1473},
    "GeO2": {"type": "network_glass", "network_former": "Ge", 
             "coordination": 4, "density_gcm3": 3.65, "Tg_K": 823},
    "B2O3": {"type": "network_glass", "network_former": "B", 
             "coordination": 3, "density_gcm3": 1.84, "Tg_K": 530},
    "P2O5": {"type": "network_glass", "network_former": "P", 
             "coordination": 4, "density_gcm3": 2.39},
    "Na2O_SiO2": {"type": "modified_glass", "modifier": "Na", 
                  "network_former": "Si", "density_gcm3": 2.40},
    "CaO_SiO2": {"type": "modified_glass", "modifier": "Ca",
                 "network_former": "Si"},
    
    # Metallic glasses
    "Zr_Cu": {"type": "metallic_glass", "composition": {"Zr": 0.5, "Cu": 0.5},
              "density_gcm3": 6.8, "Tg_K": 670},
    "Zr_Cu_Ni_Al": {"type": "bulk_metallic_glass", 
                    "composition": {"Zr": 0.55, "Cu": 0.30, "Ni": 0.05, "Al": 0.10},
                    "density_gcm3": 6.7, "Tg_K": 673, "GFA": "high"},
    "Pd_Ni_P": {"type": "metallic_glass",
                "composition": {"Pd": 0.40, "Ni": 0.40, "P": 0.20},
                "density_gcm3": 9.4, "Tg_K": 590},
    "Fe_B": {"type": "metallic_glass", 
             "composition": {"Fe": 0.80, "B": 0.20},
             "magnetic": True, "soft_magnetic": True},
    
    # Chalcogenide glasses
    "As2S3": {"type": "chalcogenide_glass", "elements": ["As", "S"],
              "density_gcm3": 3.20, "infrared_transparent": True},
    "As2Se3": {"type": "chalcogenide_glass", "elements": ["As", "Se"],
               "density_gcm3": 4.60},
    "GeSe2": {"type": "chalcogenide_glass", "elements": ["Ge", "Se"],
              "phase_change": True},
    "Ge2Sb2Te5": {"type": "phase_change", "PCM": True,
                  "elements": ["Ge", "Sb", "Te"],
                  "crystallization_fast": True},
    
    # Amorphous semiconductors
    "a_Si": {"type": "amorphous_semiconductor", "element": "Si",
             "density_gcm3": 2.29, "dangling_bonds": True},
    "a_Si_H": {"type": "hydrogenated_aSi", "elements": ["Si", "H"],
               "H_content_percent": 10, "solar_cell": True},
    "a_Ge": {"type": "amorphous_semiconductor", "element": "Ge",
             "density_gcm3": 5.1},
    "a_C": {"type": "amorphous_carbon", "element": "C",
            "sp3_fraction": 0.3, "DLC": True},
}


# Coordination number database
COORDINATION_NUMBERS = {
    "Si": 4, "Ge": 4, "B": 3, "P": 4, "As": 3,
    "O": 2, "S": 2, "Se": 2, "Te": 2,
    "Na": 6, "Ca": 6, "Mg": 6,
    "Zr": 12, "Cu": 12, "Ni": 12, "Al": 12, "Fe": 12,
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_amorphous_structure(
    material: str = "SiO2",
    n_atoms: int = 100,
    density_gcm3: Optional[float] = None,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate amorphous structure using random placement with constraints.
    
    Args:
        material: Amorphous material from database
        n_atoms: Number of atoms
        density_gcm3: Target density (uses database default if None)
        seed: Random seed
    
    Returns:
        Amorphous structure
    """
    if material not in AMORPHOUS_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material",
                      "available": list(AMORPHOUS_DATABASE.keys())}
        }
    
    info = AMORPHOUS_DATABASE[material]
    np.random.seed(seed)
    
    # Get density
    if density_gcm3 is None:
        density_gcm3 = info.get("density_gcm3", 2.5)
    
    # Calculate box size from density
    composition = info.get("composition", {})
    if not composition:
        if "elements" in info:
            n_elem = len(info["elements"])
            composition = {e: 1/n_elem for e in info["elements"]}
        elif "element" in info:
            composition = {info["element"]: 1.0}
        else:
            # Extract from formula
            if "SiO2" in material:
                composition = {"Si": 0.33, "O": 0.67}
            elif "GeO2" in material:
                composition = {"Ge": 0.33, "O": 0.67}
            elif "B2O3" in material:
                composition = {"B": 0.40, "O": 0.60}
            elif "As2S3" in material:
                composition = {"As": 0.40, "S": 0.60}
            else:
                composition = {"Si": 1.0}
    
    # Estimate atomic volume
    avg_mass = sum(comp * COORDINATION_NUMBERS.get(elem, 12) for elem, comp in composition.items())
    volume_A3 = (avg_mass * n_atoms * 1.66054) / density_gcm3  # Approximate
    box_size = volume_A3 ** (1/3)
    
    lattice = Lattice.cubic(box_size)
    
    species = []
    coords = []
    
    # Generate species based on composition
    for elem, frac in composition.items():
        n_elem = int(n_atoms * frac)
        species.extend([elem] * n_elem)
    
    # Pad to n_atoms
    while len(species) < n_atoms:
        species.append(list(composition.keys())[0])
    species = species[:n_atoms]
    
    # Generate random positions with minimum distance constraint
    min_dist = 1.5 / box_size
    
    for i in range(n_atoms):
        for attempt in range(100):
            pos = np.random.random(3)
            
            # Check minimum distance
            valid = True
            for existing in coords:
                dist = np.sqrt(sum((pos[j] - existing[j])**2 for j in range(3)))
                if dist < min_dist:
                    valid = False
                    break
            
            if valid:
                coords.append(list(pos))
                break
        else:
            # If can't place, just add anyway
            coords.append(list(np.random.random(3)))
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "material": material,
        "type": info["type"],
        "n_atoms": len(structure),
        "target_density_gcm3": density_gcm3,
        "box_size_A": round(box_size, 2),
        "Tg_K": info.get("Tg_K", 0),
        "composition": composition,
        "random_seed": seed,
        "structure": structure_to_dict(structure)
    }


def generate_melt_quench(
    elements: List[str],
    composition: List[float],
    n_atoms: int = 100,
    quench_rate: str = "fast"
) -> Dict[str, Any]:
    """
    Generate amorphous structure using simulated melt-quench.
    
    Args:
        elements: List of elements
        composition: Composition fractions
        n_atoms: Total atoms
        quench_rate: 'fast', 'medium', 'slow'
    
    Returns:
        Melt-quenched structure
    """
    if len(elements) != len(composition):
        return {"success": False, "error": {"code": "MISMATCH", "message": "Elements and composition length mismatch"}}
    
    if abs(sum(composition) - 1.0) > 0.01:
        return {"success": False, "error": {"code": "COMPOSITION", "message": "Composition must sum to 1"}}
    
    np.random.seed(42)
    
    # Estimate box size
    box_size = (n_atoms * 15) ** (1/3)
    
    lattice = Lattice.cubic(box_size)
    
    species = []
    coords = []
    
    # Build species list
    for elem, frac in zip(elements, composition):
        n = int(n_atoms * frac)
        species.extend([elem] * n)
    
    while len(species) < n_atoms:
        species.append(elements[0])
    
    # Random positions
    for _ in range(len(species)):
        coords.append(list(np.random.random(3)))
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "elements": elements,
        "composition": dict(zip(elements, composition)),
        "n_atoms": len(structure),
        "quench_rate": quench_rate,
        "simulated": True,
        "structure": structure_to_dict(structure)
    }


def generate_amorphous_interface(
    crystalline: str = "Si",
    amorphous: str = "a_Si",
    crystal_layers: int = 5,
    amorphous_thickness_A: float = 20.0
) -> Dict[str, Any]:
    """
    Generate crystalline-amorphous interface.
    
    Args:
        crystalline: Crystalline phase
        amorphous: Amorphous phase
        crystal_layers: Number of crystal layers
        amorphous_thickness_A: Thickness of amorphous region
    
    Returns:
        Interface structure
    """
    a_Si = 5.43
    interlayer = a_Si / 4
    
    total_c = crystal_layers * interlayer + amorphous_thickness_A + 10
    
    lattice = Lattice.orthorhombic(a_Si * 2, a_Si * 2, total_c)
    
    species = []
    coords = []
    
    # Crystal region
    for layer in range(crystal_layers):
        z = (layer * interlayer) / total_c
        for i in range(4):
            for j in range(4):
                x = (i + 0.5 * (layer % 2)) / 4
                y = (j + 0.5 * (layer % 2)) / 4
                species.append("Si")
                coords.append([x, y, z])
    
    # Amorphous region
    crystal_z = crystal_layers * interlayer / total_c
    n_amorphous = int(amorphous_thickness_A**3 * 0.05 / 20)
    
    np.random.seed(42)
    for _ in range(n_amorphous):
        x = np.random.random()
        y = np.random.random()
        z = crystal_z + np.random.random() * (amorphous_thickness_A / total_c)
        species.append("Si")
        coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "crystalline": crystalline,
        "amorphous": amorphous,
        "crystal_layers": crystal_layers,
        "amorphous_thickness_A": amorphous_thickness_A,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
