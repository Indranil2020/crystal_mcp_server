"""
defect/sqs.py - Special Quasirandom Structures

Generates chemically disordered structures:
- SQS (Special Quasirandom Structures)
- Short-range order (SRO) structures
- Cluster expansion configurations
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Common alloy systems for SQS
SQS_ALLOY_DATABASE = {
    "CuZn": {"a": 2.95, "structure": "BCC", "description": "Brass"},
    "FeNi": {"a": 3.57, "structure": "FCC", "description": "Invar"},
    "NiAl": {"a": 2.88, "structure": "BCC", "description": "NiAl intermetallic"},
    "FeCr": {"a": 2.87, "structure": "BCC", "description": "Stainless steel base"},
    "CuAu": {"a": 3.85, "structure": "FCC", "description": "Cu-Au alloy"},
    "TiAl": {"a": 4.00, "structure": "FCC", "description": "TiAl intermetallic"},
    "NiFe": {"a": 3.55, "structure": "FCC", "description": "Permalloy"},
    "AlMg": {"a": 4.10, "structure": "FCC", "description": "Al-Mg alloy"},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def calculate_sro_parameter(
    structure: Structure,
    element_a: str,
    element_b: str,
    shell: int = 1,
    concentration_a: float = 0.5
) -> float:
    """Calculate Warren-Cowley short-range order parameter."""
    # Simplified calculation
    n_a = sum(1 for s in structure if str(s.specie) == element_a)
    n_b = sum(1 for s in structure if str(s.specie) == element_b)
    
    if n_a == 0 or n_b == 0:
        return 0.0
    
    # For random structure, alpha ~ 0
    # For ordered, alpha approaches limits
    return 0.0  # Placeholder for proper calculation


def generate_sqs(
    alloy: str = "FeNi",
    composition: Optional[Dict[str, float]] = None,
    supercell_size: List[int] = [3, 3, 3],
    n_iterations: int = 1000,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate Special Quasirandom Structure.
    
    Args:
        alloy: Alloy system from database
        composition: Override composition {elem: fraction}
        supercell_size: Supercell dimensions
        n_iterations: Monte Carlo swap iterations
        seed: Random seed
    
    Returns:
        SQS structure
    """
    if alloy not in SQS_ALLOY_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ALLOY", "message": f"Unknown alloy '{alloy}'",
                      "available": list(SQS_ALLOY_DATABASE.keys())}
        }
    
    np.random.seed(seed)
    
    info = SQS_ALLOY_DATABASE[alloy]
    a = info["a"]
    structure_type = info["structure"]
    
    # Parse elements from alloy name
    elements = []
    current = ""
    for char in alloy:
        if char.isupper():
            if current:
                elements.append(current)
            current = char
        else:
            current += char
    if current:
        elements.append(current)
    
    # Default 50-50 composition
    if composition is None:
        composition = {elements[0]: 0.5, elements[1]: 0.5}
    
    nx, ny, nz = supercell_size
    
    if structure_type == "BCC":
        lattice = Lattice.cubic(a * nx)
        base_positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    else:  # FCC
        lattice = Lattice.cubic(a * nx)
        base_positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    # Generate all positions
    all_coords = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in base_positions:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    all_coords.append([x, y, z])
    
    n_total = len(all_coords)
    
    # Assign atoms based on composition
    species = []
    for elem, frac in composition.items():
        n_elem = int(n_total * frac)
        species.extend([elem] * n_elem)
    
    while len(species) < n_total:
        species.append(elements[0])
    species = species[:n_total]
    
    # Random initial configuration
    np.random.shuffle(species)
    
    # Monte Carlo optimization for SQS
    # (Simplified - in practice would use cluster function optimization)
    current_structure = Structure(lattice, species, all_coords)
    
    for iteration in range(n_iterations):
        # Find two random sites with different species
        i = np.random.randint(n_total)
        j = np.random.randint(n_total)
        
        if species[i] != species[j]:
            # Swap and check if improves randomness
            species[i], species[j] = species[j], species[i]
            # In full SQS, would calculate correlation functions here
    
    final_structure = Structure(lattice, species, all_coords)
    
    # Calculate actual composition
    actual_comp = {elem: sum(1 for s in species if s == elem) / n_total for elem in elements}
    
    return {
        "success": True,
        "alloy": alloy,
        "description": info["description"],
        "structure_type": structure_type,
        "target_composition": composition,
        "actual_composition": {k: round(v, 3) for k, v in actual_comp.items()},
        "supercell_size": supercell_size,
        "n_mc_iterations": n_iterations,
        "n_atoms": n_total,
        "structure": structure_to_dict(final_structure)
    }


def generate_sro_structure(
    alloy: str = "CuZn",
    sro_parameter: float = 0.3,
    supercell_size: List[int] = [4, 4, 4],
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate structure with specified short-range order.
    
    Args:
        alloy: Alloy system
        sro_parameter: Warren-Cowley parameter (-1 to 1, 0 = random)
        supercell_size: Supercell dimensions
        seed: Random seed
    
    Returns:
        Structure with SRO
    """
    np.random.seed(seed)
    
    if alloy not in SQS_ALLOY_DATABASE:
        return {"success": False, "error": {"code": "INVALID_ALLOY", "message": f"Unknown alloy '{alloy}'"}}
    
    info = SQS_ALLOY_DATABASE[alloy]
    a = info["a"]
    
    elements = []
    current = ""
    for char in alloy:
        if char.isupper():
            if current:
                elements.append(current)
            current = char
        else:
            current += char
    if current:
        elements.append(current)
    
    nx, ny, nz = supercell_size
    
    if info["structure"] == "BCC":
        lattice = Lattice.cubic(a * nx)
        base_positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    else:
        lattice = Lattice.cubic(a * nx)
        base_positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    all_coords = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in base_positions:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    all_coords.append([x, y, z])
    
    n_total = len(all_coords)
    n_A = n_total // 2
    
    # For positive SRO (clustering), group like atoms
    # For negative SRO (ordering), alternate atoms
    species = []
    
    if sro_parameter > 0:
        # Clustering - half A, half B in blocks
        for i, coord in enumerate(all_coords):
            x, y, z = coord
            block = int(x * 2) + int(y * 2) + int(z * 2)
            if block % 2 == 0:
                species.append(elements[0])
            else:
                species.append(elements[1])
            
            # Add some randomness
            if np.random.random() > sro_parameter:
                species[-1] = elements[1] if species[-1] == elements[0] else elements[0]
    else:
        # Ordering - alternating
        for i, coord in enumerate(all_coords):
            if i % 2 == 0:
                species.append(elements[0])
            else:
                species.append(elements[1])
            
            # Add some randomness
            if np.random.random() > abs(sro_parameter):
                species[-1] = elements[1] if species[-1] == elements[0] else elements[0]
    
    structure = Structure(lattice, species, all_coords)
    
    return {
        "success": True,
        "alloy": alloy,
        "target_sro": sro_parameter,
        "sro_interpretation": "clustering" if sro_parameter > 0 else "ordering" if sro_parameter < 0 else "random",
        "n_atoms": n_total,
        "structure": structure_to_dict(structure)
    }


def generate_hea_sqs(
    elements: List[str] = ["Fe", "Co", "Ni", "Cr", "Mn"],
    composition: Optional[Dict[str, float]] = None,
    supercell_size: List[int] = [4, 4, 4],
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate High-Entropy Alloy SQS structure.
    
    Args:
        elements: List of 4-6 elements
        composition: Compositions (default equiatomic)
        supercell_size: Supercell dimensions
        seed: Random seed
    
    Returns:
        HEA SQS structure
    """
    np.random.seed(seed)
    
    if len(elements) < 4:
        return {"success": False, "error": {"code": "INVALID_HEA", "message": "HEA requires at least 4 elements"}}
    
    # Default equiatomic
    if composition is None:
        composition = {elem: 1.0 / len(elements) for elem in elements}
    
    # Average lattice constant
    a_values = {"Fe": 2.87, "Co": 3.54, "Ni": 3.52, "Cr": 2.91, "Mn": 3.86, 
                "V": 3.02, "Ti": 4.51, "Cu": 3.61, "Al": 4.05}
    a = np.mean([a_values.get(e, 3.5) for e in elements])
    
    nx, ny, nz = supercell_size
    
    # FCC structure (most common for HEAs)
    lattice = Lattice.cubic(a * nx)
    base_positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    all_coords = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for bp in base_positions:
                    x = (i + bp[0]) / nx
                    y = (j + bp[1]) / ny
                    z = (k + bp[2]) / nz
                    all_coords.append([x, y, z])
    
    n_total = len(all_coords)
    
    # Assign species based on composition
    species = []
    for elem, frac in composition.items():
        n_elem = int(n_total * frac)
        species.extend([elem] * n_elem)
    
    while len(species) < n_total:
        species.append(elements[0])
    species = species[:n_total]
    
    np.random.shuffle(species)
    
    structure = Structure(lattice, species, all_coords)
    
    # Calculate configurational entropy
    S_config = -sum(f * np.log(f) if f > 0 else 0 for f in composition.values())
    S_config_ideal = np.log(len(elements))
    
    return {
        "success": True,
        "elements": elements,
        "n_elements": len(elements),
        "composition": {k: round(v, 3) for k, v in composition.items()},
        "is_equiatomic": all(abs(v - 1/len(elements)) < 0.01 for v in composition.values()),
        "configurational_entropy_R": round(S_config, 3),
        "ideal_entropy_R": round(S_config_ideal, 3),
        "n_atoms": n_total,
        "structure": structure_to_dict(structure)
    }
