"""
bulk/zeolites.py - Zeolite Structures

Comprehensive zeolite framework generation:
- Common zeolite framework types
- Guest molecule insertion
- Cation exchange
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Zeolite framework database
ZEOLITE_DATABASE = {
    # Common zeolites
    "MFI": {
        "name": "ZSM-5", "channels": "10-ring", "a": 20.09, "b": 19.74, "c": 13.14,
        "pore_diameter_A": 5.5, "Si_Al_ratio": 30, "catalytic": True
    },
    "FAU": {
        "name": "Faujasite", "channels": "12-ring", "a": 24.74, "cubic": True,
        "pore_diameter_A": 7.4, "supercage": True, "catalytic": True
    },
    "LTA": {
        "name": "Linde Type A", "channels": "8-ring", "a": 11.92, "cubic": True,
        "pore_diameter_A": 4.1, "molecular_sieve": True
    },
    "BEA": {
        "name": "Beta", "channels": "12-ring", "a": 12.63, "c": 26.41,
        "pore_diameter_A": 6.7, "disordered": True
    },
    "MOR": {
        "name": "Mordenite", "channels": "12-ring", "a": 18.09, "b": 20.52, "c": 7.52,
        "pore_diameter_A": 7.0
    },
    "CHA": {
        "name": "Chabazite", "channels": "8-ring", "a": 9.42, "alpha": 94.2,
        "pore_diameter_A": 3.8, "SCR_catalyst": True
    },
    "ZSM-22": {
        "framework": "TON", "channels": "10-ring", "a": 13.86, "b": 17.42, "c": 5.04,
        "pore_diameter_A": 4.6, "1D_channel": True
    },
    "ZSM-12": {
        "framework": "MTW", "channels": "12-ring", "a": 25.55, "b": 5.26, "c": 12.11,
        "pore_diameter_A": 6.0
    },
    
    # Specialty zeolites
    "SSZ-13": {
        "framework": "CHA", "channels": "8-ring", "high_Si": True,
        "NOx_reduction": True
    },
    "TS-1": {
        "framework": "MFI", "heteroatom": "Ti", "oxidation_catalyst": True
    },
    "AlPO-5": {
        "framework": "AFI", "aluminophosphate": True, "channels": "12-ring",
        "pore_diameter_A": 7.3
    },
    "SAPO-34": {
        "framework": "CHA", "silicoaluminophosphate": True,
        "MTO_catalyst": True
    },
}


# Extra-framework cation database
CATION_DATABASE = {
    "Na": {"charge": 1, "radius_A": 1.02, "hydration_enthalpy_kJ_mol": -405},
    "K": {"charge": 1, "radius_A": 1.38, "hydration_enthalpy_kJ_mol": -321},
    "Ca": {"charge": 2, "radius_A": 1.00, "hydration_enthalpy_kJ_mol": -1592},
    "Mg": {"charge": 2, "radius_A": 0.72, "hydration_enthalpy_kJ_mol": -1920},
    "Ba": {"charge": 2, "radius_A": 1.35},
    "NH4": {"charge": 1, "radius_A": 1.48, "molecular": True},
    "H": {"charge": 1, "radius_A": 0.00, "Bronsted_acid": True},
    "Cu": {"charge": 2, "radius_A": 0.73, "SCR_active": True},
    "Fe": {"charge": 3, "radius_A": 0.65, "redox_active": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_zeolite(
    framework: str = "MFI",
    cation: str = "Na",
    Si_Al_ratio: float = 30.0,
    supercell: List[int] = [1, 1, 1]
) -> Dict[str, Any]:
    """
    Generate zeolite structure.
    
    Args:
        framework: Zeolite framework type
        cation: Extra-framework cation
        Si_Al_ratio: Si/Al ratio
        supercell: Supercell dimensions
    
    Returns:
        Zeolite structure
    """
    if framework not in ZEOLITE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_FRAMEWORK", "message": f"Unknown framework",
                      "available": list(ZEOLITE_DATABASE.keys())}
        }
    
    info = ZEOLITE_DATABASE[framework]
    
    # Get lattice parameters
    a = info.get("a", 10.0)
    b = info.get("b", a)
    c = info.get("c", a)
    
    if info.get("cubic"):
        lattice = Lattice.cubic(a)
    else:
        alpha = info.get("alpha", 90)
        lattice = Lattice.from_parameters(a, b, c, alpha, 90, 90)
    
    species = []
    coords = []
    
    # Generate simplified framework (T-atoms at nodes)
    np.random.seed(42)
    
    n_T_atoms = int(a * b * c / 50)  # Approximate T-atom density
    n_Al = max(1, int(n_T_atoms / (Si_Al_ratio + 1)))
    n_Si = n_T_atoms - n_Al
    
    # Place T-atoms
    for i in range(n_Si):
        pos = np.random.random(3)
        species.append("Si")
        coords.append(list(pos))
    
    for i in range(n_Al):
        pos = np.random.random(3)
        species.append("Al")
        coords.append(list(pos))
    
    # Add bridging oxygens (2 per T-atom average)
    n_O = int(n_T_atoms * 2)
    for i in range(n_O):
        pos = np.random.random(3)
        species.append("O")
        coords.append(list(pos))
    
    # Add cations to balance Al charge
    if cation != "H" and n_Al > 0:
        cation_charge = CATION_DATABASE.get(cation, {}).get("charge", 1)
        n_cations = max(1, int(n_Al / cation_charge))
        
        for i in range(n_cations):
            pos = np.random.random(3)
            species.append(cation)
            coords.append(list(pos))
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "framework": framework,
        "name": info.get("name", framework),
        "channels": info.get("channels", "unknown"),
        "pore_diameter_A": info.get("pore_diameter_A", 0),
        "Si_Al_ratio": Si_Al_ratio,
        "cation": cation,
        "n_atoms": len(structure),
        "is_catalytic": info.get("catalytic", False),
        "structure": structure_to_dict(structure)
    }


def generate_zeolite_with_guest(
    framework: str = "FAU",
    guest_molecule: str = "H2O",
    loading: float = 0.5
) -> Dict[str, Any]:
    """
    Generate zeolite with guest molecules.
    
    Args:
        framework: Zeolite framework
        guest_molecule: Guest molecule (H2O, CO2, NH3, etc.)
        loading: Fractional loading (0-1)
    
    Returns:
        Zeolite with guests
    """
    # First generate empty framework
    result = generate_zeolite(framework)
    
    if not result["success"]:
        return result
    
    atoms = result["structure"]["atoms"]
    lattice = result["structure"]["lattice"]
    
    # Estimate max guest molecules
    pore_volume = ZEOLITE_DATABASE.get(framework, {}).get("pore_diameter_A", 5) ** 3 / 10
    n_guests = int(pore_volume * loading * 5)
    
    # Add guest molecules (simplified as single atoms)
    guest_elem = guest_molecule[0] if guest_molecule[0].isalpha() else guest_molecule[:2]
    
    np.random.seed(123)
    for _ in range(n_guests):
        pos = np.random.random(3).tolist()
        atoms.append({"element": guest_elem, "coords": pos, "is_guest": True})
    
    result["guest_molecule"] = guest_molecule
    result["loading_fraction"] = loading
    result["n_guest_molecules"] = n_guests
    result["n_atoms"] = len(atoms)
    result["structure"]["atoms"] = atoms
    
    return result


def ion_exchange_zeolite(
    framework: str = "FAU",
    original_cation: str = "Na",
    new_cation: str = "Ca",
    exchange_fraction: float = 0.5
) -> Dict[str, Any]:
    """
    Perform ion exchange on zeolite.
    
    Args:
        framework: Zeolite framework
        original_cation: Original cation
        new_cation: New cation after exchange
        exchange_fraction: Fraction exchanged
    
    Returns:
        Ion-exchanged zeolite
    """
    result = generate_zeolite(framework, original_cation)
    
    if not result["success"]:
        return result
    
    atoms = result["structure"]["atoms"]
    
    # Replace cations
    np.random.seed(42)
    for atom in atoms:
        if atom["element"] == original_cation:
            if np.random.random() < exchange_fraction:
                atom["element"] = new_cation
                atom["is_exchanged"] = True
    
    result["ion_exchange"] = {
        "original": original_cation,
        "new": new_cation,
        "fraction": exchange_fraction
    }
    result["structure"]["atoms"] = atoms
    
    return result
