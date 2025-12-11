"""
defect/amorphous.py - Amorphous and Glassy Structures

Generates amorphous/disordered structures:
- a-Si, a-Ge, a-C (amorphous semiconductors)
- a-SiO2 (silica glass)
- a-IGZO (amorphous oxide semiconductors)
- Metallic glasses
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Amorphous material database
AMORPHOUS_DATABASE = {
    "a-Si": {
        "elements": {"Si": 1.0}, "density": 2.3, "coordination": 4,
        "bond_length": 2.35, "description": "Amorphous silicon"
    },
    "a-Ge": {
        "elements": {"Ge": 1.0}, "density": 5.3, "coordination": 4,
        "bond_length": 2.45, "description": "Amorphous germanium"
    },
    "a-C": {
        "elements": {"C": 1.0}, "density": 2.0, "coordination": 3.5,
        "bond_length": 1.54, "description": "Amorphous carbon (diamond-like)"
    },
    "a-SiO2": {
        "elements": {"Si": 0.33, "O": 0.67}, "density": 2.2, "coordination": 4,
        "bond_length": 1.61, "description": "Fused silica"
    },
    "a-IGZO": {
        "elements": {"In": 0.4, "Ga": 0.2, "Zn": 0.2, "O": 0.2}, "density": 6.4,
        "bond_length": 2.0, "description": "Amorphous InGaZnO"
    },
    "a-GST": {
        "elements": {"Ge": 0.225, "Sb": 0.225, "Te": 0.55}, "density": 6.4,
        "bond_length": 2.8, "description": "Amorphous Ge2Sb2Te5 (phase change)"
    },
    "a-Li3PS4": {
        "elements": {"Li": 0.375, "P": 0.125, "S": 0.5}, "density": 1.9,
        "bond_length": 2.0, "description": "Amorphous lithium thiophosphate"
    },
    # Metallic glasses
    "Zr-Cu-Al_MG": {
        "elements": {"Zr": 0.55, "Cu": 0.30, "Al": 0.15}, "density": 6.5,
        "bond_length": 2.8, "description": "Zr-Cu-Al metallic glass"
    },
    "Pd-Ni-P_MG": {
        "elements": {"Pd": 0.4, "Ni": 0.4, "P": 0.2}, "density": 9.5,
        "bond_length": 2.5, "description": "Pd-Ni-P metallic glass"
    },
    "Fe-B_MG": {
        "elements": {"Fe": 0.8, "B": 0.2}, "density": 7.3,
        "bond_length": 2.5, "description": "Fe-B amorphous alloy"
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


def generate_amorphous_structure(
    material: str = "a-Si",
    n_atoms: int = 64,
    density: Optional[float] = None,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate amorphous structure using random packing.
    
    Args:
        material: Material type from database
        n_atoms: Number of atoms
        density: Override density (g/cm³)
        seed: Random seed
    
    Returns:
        Amorphous structure
    """
    if material not in AMORPHOUS_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": f"Unknown material '{material}'",
                      "available": list(AMORPHOUS_DATABASE.keys())}
        }
    
    np.random.seed(seed)
    info = AMORPHOUS_DATABASE[material]
    
    # Use provided or database density
    rho = density if density else info["density"]
    
    # Calculate box size from density
    # Approximate average atomic mass
    elements = info["elements"]
    avg_mass = sum(
        frac * {"Si": 28, "Ge": 73, "C": 12, "O": 16, "In": 115, "Ga": 70, "Zn": 65,
                "Te": 128, "Sb": 122, "Li": 7, "P": 31, "S": 32, "Zr": 91, "Cu": 64,
                "Al": 27, "Pd": 106, "Ni": 59, "Fe": 56, "B": 11}.get(elem, 50)
        for elem, frac in elements.items()
    )
    
    # Volume in Å³
    volume = n_atoms * avg_mass / (rho * 0.6022)  # 0.6022 converts g/cm³ to amu/Å³
    box_size = volume ** (1/3)
    
    lattice = Lattice.cubic(box_size)
    
    # Generate atom types based on composition
    atom_types = []
    for elem, frac in elements.items():
        n_elem = max(1, int(n_atoms * frac))
        atom_types.extend([elem] * n_elem)
    
    # Adjust to exact number
    while len(atom_types) < n_atoms:
        atom_types.append(list(elements.keys())[0])
    atom_types = atom_types[:n_atoms]
    np.random.shuffle(atom_types)
    
    # Random positions with minimum distance constraint
    min_dist = info["bond_length"] * 0.8
    positions = []
    max_attempts = 1000
    
    for i in range(n_atoms):
        for attempt in range(max_attempts):
            pos = np.random.uniform(0, 1, 3)
            
            # Check minimum distance to existing atoms
            valid = True
            for existing in positions:
                diff = pos - existing
                # Periodic boundary
                diff = diff - np.round(diff)
                dist = np.linalg.norm(diff * box_size)
                if dist < min_dist:
                    valid = False
                    break
            
            if valid:
                positions.append(pos)
                break
        else:
            # If couldn't place with constraint, place randomly
            positions.append(np.random.uniform(0, 1, 3))
    
    structure = Structure(lattice, atom_types, positions)
    
    return {
        "success": True,
        "material": material,
        "description": info["description"],
        "n_atoms": n_atoms,
        "density_g_cm3": rho,
        "box_size_angstrom": round(box_size, 2),
        "composition": {k: round(v * n_atoms) for k, v in elements.items()},
        "structure": structure_to_dict(structure)
    }


def generate_melt_quenched_amorphous(
    material: str = "a-Si",
    n_atoms: int = 64,
    quench_rate: str = "fast",
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate amorphous structure simulating melt-quench process.
    
    Args:
        material: Material type
        n_atoms: Number of atoms
        quench_rate: 'fast', 'medium', 'slow'
        seed: Random seed
    
    Returns:
        Melt-quenched amorphous structure
    """
    # For now, use random structure with disorder parameter
    result = generate_amorphous_structure(material, n_atoms, seed=seed)
    
    if result["success"]:
        # Add quench-specific metadata
        result["generation_method"] = "melt_quench"
        result["quench_rate"] = quench_rate
        
        # Disorder increases with faster quench
        disorder_factor = {"fast": 1.2, "medium": 1.0, "slow": 0.8}[quench_rate]
        result["disorder_factor"] = disorder_factor
    
    return result
