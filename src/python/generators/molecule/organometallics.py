"""
molecule/organometallics.py - Organometallic Compounds

Generates organometallic molecular structures:
- Metallocenes (ferrocene, cobaltocene, etc.)
- Piano-stool complexes
- Carbenes (NHC, Grubbs)
- Sandwich and half-sandwich complexes
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Metallocene database
METALLOCENE_DATABASE = {
    "ferrocene": {
        "metal": "Fe", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.66, "ring_ring_distance": 3.32,
        "eclipsed": True, "spin": 0, "oxidation": 2
    },
    "cobaltocene": {
        "metal": "Co", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.73, "ring_ring_distance": 3.46,
        "eclipsed": True, "spin": 0.5, "oxidation": 2
    },
    "nickelocene": {
        "metal": "Ni", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.83, "ring_ring_distance": 3.66,
        "eclipsed": True, "spin": 1, "oxidation": 2
    },
    "ruthenocene": {
        "metal": "Ru", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.82, "ring_ring_distance": 3.64,
        "eclipsed": True, "spin": 0, "oxidation": 2
    },
    "osmocene": {
        "metal": "Os", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.85, "ring_ring_distance": 3.70,
        "eclipsed": True, "spin": 0, "oxidation": 2
    },
    "chromocene": {
        "metal": "Cr", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.80, "ring_ring_distance": 3.60,
        "eclipsed": False, "spin": 1, "oxidation": 2
    },
    "manganocene": {
        "metal": "Mn", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.77, "ring_ring_distance": 3.54,
        "eclipsed": False, "spin": 2.5, "oxidation": 2
    },
    "vanadocene": {
        "metal": "V", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 1.92, "ring_ring_distance": 3.84,
        "eclipsed": True, "spin": 1.5, "oxidation": 2
    },
    "titanocene": {
        "metal": "Ti", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 2.06, "ring_ring_distance": 4.12,
        "eclipsed": True, "spin": 1, "oxidation": 2
    },
    "zirconocene": {
        "metal": "Zr", "ligand": "Cp", "n_ligands": 2,
        "metal_ligand_distance": 2.21, "ring_ring_distance": 4.42,
        "eclipsed": True, "spin": 0, "oxidation": 2
    },
}


# Piano-stool complexes
PIANO_STOOL_DATABASE = {
    "CpFe(CO)2I": {
        "metal": "Fe", "ring": "Cp", "carbonyls": 2, "halide": "I",
        "description": "Iron piano-stool"
    },
    "CpMn(CO)3": {
        "metal": "Mn", "ring": "Cp", "carbonyls": 3, "halide": None,
        "description": "Cymantrene"
    },
    "CpRe(CO)3": {
        "metal": "Re", "ring": "Cp", "carbonyls": 3, "halide": None,
        "description": "Rhenium cymantrenyl"
    },
    "Cp*Rh(CO)2": {
        "metal": "Rh", "ring": "Cp*", "carbonyls": 2, "halide": None,
        "description": "Pentamethyl rhodium dicarbonyl"
    },
    "Cp*Ir(CO)2": {
        "metal": "Ir", "ring": "Cp*", "carbonyls": 2, "halide": None,
        "description": "Pentamethyl iridium dicarbonyl"
    },
}


# Carbene catalysts
CARBENE_CATALYSTS = {
    "Grubbs-1": {
        "metal": "Ru", "ligands": ["PCy3", "PCy3", "Cl", "Cl", "CHPh"],
        "generation": 1, "description": "First-generation Grubbs catalyst"
    },
    "Grubbs-2": {
        "metal": "Ru", "ligands": ["NHC", "PCy3", "Cl", "Cl", "CHPh"],
        "generation": 2, "description": "Second-generation Grubbs catalyst"
    },
    "Hoveyda-Grubbs-2": {
        "metal": "Ru", "ligands": ["NHC", "ether-chelate", "Cl", "Cl"],
        "generation": 2, "description": "Hoveyda-Grubbs catalyst"
    },
    "Schrock": {
        "metal": "Mo", "ligands": ["NAr", "OR", "OR", "CHR"],
        "generation": 0, "description": "Schrock carbene catalyst"
    },
}


def generate_cyclopentadienyl_ring(radius: float = 1.21, z: float = 0, rotation: float = 0) -> List[Dict[str, Any]]:
    """Generate Cp ring atom positions."""
    atoms = []
    for i in range(5):
        angle = 2 * np.pi * i / 5 + np.radians(rotation)
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        atoms.append({"element": "C", "cartesian": [x, y, z]})
    return atoms


def generate_metallocene(
    name: str = "ferrocene",
    conformation: str = "eclipsed",
    add_hydrogens: bool = True
) -> Dict[str, Any]:
    """
    Generate metallocene structure.
    
    Args:
        name: Metallocene name
        conformation: 'eclipsed' or 'staggered'
        add_hydrogens: Whether to add H atoms to Cp rings
    
    Returns:
        Metallocene structure
    """
    if name not in METALLOCENE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_NAME", "message": f"Unknown metallocene '{name}'",
                      "available": list(METALLOCENE_DATABASE.keys())}
        }
    
    info = METALLOCENE_DATABASE[name]
    metal = info["metal"]
    half_dist = info["ring_ring_distance"] / 2
    
    atoms = []
    
    # Metal center
    atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    
    # Generate Cp rings
    cp_radius = 1.21  # C-C distance in Cp
    
    # Top ring
    rotation_top = 0
    atoms.extend(generate_cyclopentadienyl_ring(cp_radius, half_dist, rotation_top))
    
    # Bottom ring
    rotation_bottom = 36 if conformation == "staggered" else 0  # 36° rotation for staggered
    atoms.extend(generate_cyclopentadienyl_ring(cp_radius, -half_dist, rotation_bottom))
    
    # Add hydrogens if requested
    if add_hydrogens:
        h_radius = 2.18  # C-H projected radius
        
        for i in range(5):
            # Top ring H
            angle = 2 * np.pi * i / 5
            x = h_radius * np.cos(angle)
            y = h_radius * np.sin(angle)
            atoms.append({"element": "H", "cartesian": [x, y, half_dist]})
            
            # Bottom ring H
            angle_bottom = angle + np.radians(36 if conformation == "staggered" else 0)
            x = h_radius * np.cos(angle_bottom)
            y = h_radius * np.sin(angle_bottom)
            atoms.append({"element": "H", "cartesian": [x, y, -half_dist]})
    
    return {
        "success": True,
        "name": name,
        "metal": metal,
        "conformation": conformation,
        "spin": info["spin"],
        "oxidation_state": info["oxidation"],
        "metal_ring_distance_angstrom": info["metal_ligand_distance"],
        "ring_ring_distance_angstrom": info["ring_ring_distance"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_piano_stool(
    name: str = "CpMn(CO)3",
    add_hydrogens: bool = True
) -> Dict[str, Any]:
    """
    Generate piano-stool complex.
    
    Args:
        name: Complex name
        add_hydrogens: Add H atoms to Cp ring
    
    Returns:
        Piano-stool structure
    """
    if name not in PIANO_STOOL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_NAME", "message": f"Unknown complex '{name}'",
                      "available": list(PIANO_STOOL_DATABASE.keys())}
        }
    
    info = PIANO_STOOL_DATABASE[name]
    metal = info["metal"]
    n_co = info["carbonyls"]
    
    atoms = []
    
    # Metal center at origin
    atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    
    # Cp ring above metal
    z_cp = 1.7  # Typical M-Cp distance
    atoms.extend(generate_cyclopentadienyl_ring(1.21, z_cp, 0))
    
    # Carbonyl ligands below
    co_distance = 1.8  # M-CO distance
    co_bond = 1.14  # C=O bond
    
    for i in range(n_co):
        angle = 2 * np.pi * i / n_co
        x = co_distance * np.cos(angle) * 0.5
        y = co_distance * np.sin(angle) * 0.5
        z = -co_distance * 0.866
        
        # C of CO
        atoms.append({"element": "C", "cartesian": [x, y, z]})
        # O of CO
        atoms.append({"element": "O", "cartesian": [x * 1.6, y * 1.6, z - co_bond * 0.5]})
    
    # Add halide if present
    if info["halide"]:
        atoms.append({"element": info["halide"], "cartesian": [2.5, 0, -0.5]})
    
    # Add H to Cp
    if add_hydrogens:
        h_radius = 2.18
        for i in range(5):
            angle = 2 * np.pi * i / 5
            x = h_radius * np.cos(angle)
            y = h_radius * np.sin(angle)
            atoms.append({"element": "H", "cartesian": [x, y, z_cp]})
    
    return {
        "success": True,
        "name": name,
        "metal": metal,
        "description": info["description"],
        "n_carbonyls": n_co,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_nhc_carbene(
    metal: str = "Au",
    substituents: str = "Mes",
    oxidation: int = 1
) -> Dict[str, Any]:
    """
    Generate N-heterocyclic carbene complex.
    
    Args:
        metal: Metal center
        substituents: N-substituents (Mes, iPr, etc.)
        oxidation: Metal oxidation state
    
    Returns:
        NHC complex structure
    """
    atoms = []
    
    # Metal at origin
    atoms.append({"element": metal, "cartesian": [0, 0, 0]})
    
    # Imidazolylidene ring
    ring_dist = 2.0  # M-C distance
    
    # Carbene carbon
    atoms.append({"element": "C", "cartesian": [ring_dist, 0, 0]})
    
    # N atoms
    atoms.append({"element": "N", "cartesian": [ring_dist + 0.65, 0.95, 0]})
    atoms.append({"element": "N", "cartesian": [ring_dist + 0.65, -0.95, 0]})
    
    # Back C-C of imidazole
    atoms.append({"element": "C", "cartesian": [ring_dist + 1.7, 0.65, 0]})
    atoms.append({"element": "C", "cartesian": [ring_dist + 1.7, -0.65, 0]})
    
    # Add mesityl or iPr groups (simplified as single atoms for now)
    if substituents == "Mes":
        atoms.append({"element": "C", "cartesian": [ring_dist + 0.65, 2.3, 0]})
        atoms.append({"element": "C", "cartesian": [ring_dist + 0.65, -2.3, 0]})
    
    return {
        "success": True,
        "metal": metal,
        "ligand": "NHC",
        "substituents": substituents,
        "oxidation_state": oxidation,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
