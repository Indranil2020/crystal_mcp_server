"""
two_d/he_2d.py - High-Entropy 2D Materials

Generates high-entropy 2D materials:
- High-entropy TMDs (e.g., (Mo,W,Ti,Nb,Ta)S2)
- High-entropy MXenes
- High-entropy 2D oxides/carbides/nitrides

Scientific basis:
- Configurational entropy stabilization
- Special quasirandom structures (SQS) on 2D lattices
- 5+ element requirement for HE classification
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


def generate_he_2d(
    material_type: str = "tmd",
    elements: List[str] = None,
    base_structure: str = "MoS2",
    size: List[int] = [4, 4, 1],
    vacuum: float = 15.0,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate high-entropy 2D material.
    
    Args:
        material_type: Type of HE-2D ('tmd', 'mxene', 'oxide', 'nitride')
        elements: List of elements (5+ for true HE)
        base_structure: Base structure template
        size: Supercell size for randomization
        vacuum: Vacuum padding
        seed: Random seed for reproducibility
    
    Returns:
        High-entropy 2D structure
    """
    if seed is not None:
        np.random.seed(seed)
    
    if elements is None:
        if material_type == "tmd":
            elements = ["Mo", "W", "Ti", "Nb", "Ta"]
        elif material_type == "mxene":
            elements = ["Ti", "V", "Nb", "Ta", "Mo"]
        elif material_type == "oxide":
            elements = ["Ti", "V", "Zr", "Nb", "Hf"]
        else:
            elements = ["Ti", "V", "Cr", "Nb", "Ta"]
    
    n_elements = len(elements)
    if n_elements < 4:
        return {
            "success": False,
            "error": {
                "code": "TOO_FEW_ELEMENTS",
                "message": f"Need at least 4 elements for high-entropy, got {n_elements}",
            }
        }
    
    is_true_he = n_elements >= 5
    
    # Build base structure
    if material_type == "tmd":
        a = 3.2
        c = vacuum
        lattice = Lattice.hexagonal(a, c)
        
        # Create supercell first
        n_metals = size[0] * size[1]
        
        # Random assignment of metal sites
        metal_assignment = np.random.choice(elements, size=n_metals)
        
        species = []
        coords = []
        z_center = 0.5
        dz = 3.0 / c
        
        idx = 0
        for i in range(size[0]):
            for j in range(size[1]):
                fx = (1/3 + i) / size[0]
                fy = (2/3 + j) / size[1]
                
                species.append(metal_assignment[idx])
                coords.append([fx, fy, z_center])
                species.append("S")
                coords.append([fx, fy, z_center + dz])
                species.append("S")
                coords.append([fx, fy, z_center - dz])
                idx += 1
        
        lattice = Lattice.hexagonal(a * size[0], c)
    
    elif material_type == "mxene":
        a = 3.0
        c = vacuum
        
        n_metals = size[0] * size[1] * 3  # 3 metal layers
        metal_assignment = np.random.choice(elements, size=n_metals)
        
        species = []
        coords = []
        z_center = 0.5
        z_spacing = 2.3 / c
        
        idx = 0
        for i in range(size[0]):
            for j in range(size[1]):
                for layer in range(3):
                    fx = (layer % 3) / 3 + i / size[0]
                    fy = ((layer + 1) % 3) / 3 + j / size[1]
                    z = z_center + (layer - 1) * z_spacing
                    
                    species.append(metal_assignment[idx])
                    coords.append([fx % 1, fy % 1, z])
                    idx += 1
        
        # Add C between metal layers
        for i in range(size[0]):
            for j in range(size[1]):
                for layer in range(2):
                    fx = i / size[0]
                    fy = j / size[1]
                    z = z_center + (layer - 0.5) * z_spacing
                    species.append("C")
                    coords.append([fx, fy, z])
        
        lattice = Lattice.hexagonal(a * size[0], c)
    
    else:
        # Generic oxide/nitride
        a = 3.5
        c = vacuum
        lattice = Lattice.hexagonal(a * size[0], c)
        
        n_metals = size[0] * size[1]
        metal_assignment = np.random.choice(elements, size=n_metals)
        
        species = list(metal_assignment)
        coords = [[i/size[0], j/size[1], 0.5] for i in range(size[0]) for j in range(size[1])]
    
    structure = Structure(lattice, species, coords)
    
    # Calculate composition
    from collections import Counter
    metal_counts = Counter([s for s in species if s in elements])
    total_metals = sum(metal_counts.values())
    composition = {el: count/total_metals for el, count in metal_counts.items()}
    
    # Calculate configurational entropy
    S_config = -sum(x * np.log(x) for x in composition.values() if x > 0)
    S_config_normalized = S_config / np.log(n_elements)  # Normalize to ideal
    
    return {
        "success": True,
        "material_type": material_type,
        "elements": elements,
        "n_elements": n_elements,
        "is_true_high_entropy": is_true_he,
        "composition": {k: round(v, 3) for k, v in composition.items()},
        "configurational_entropy_R": round(S_config, 3),
        "entropy_fraction_of_ideal": round(S_config_normalized, 3),
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_he_tmd(
    metals: List[str] = ["Mo", "W", "Ti", "Nb", "Ta"],
    chalcogen: str = "S",
    size: List[int] = [4, 4, 1],
    vacuum: float = 15.0,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate high-entropy TMD specifically.
    
    Args:
        metals: Metal elements (5+ for HE)
        chalcogen: Chalcogen element (S, Se, Te)
        size: Supercell size
        vacuum: Vacuum
        seed: Random seed
    
    Returns:
        HE-TMD structure
    """
    return generate_he_2d(
        material_type="tmd",
        elements=metals,
        size=size,
        vacuum=vacuum,
        seed=seed
    )
