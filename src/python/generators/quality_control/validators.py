"""
quality_control/validators.py - Structure Validation

Validates structure quality:
- Stoichiometry check
- Bond valence analysis
- Charge neutrality
- Duplicate removal
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Common oxidation states for bond valence
OXIDATION_STATES = {
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
    "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2,
    "Al": 3, "Ga": 3, "In": 3,
    "Si": 4, "Ge": 4, "Sn": 2,
    "Ti": 4, "Zr": 4, "Hf": 4,
    "V": 5, "Nb": 5, "Ta": 5,
    "Cr": 3, "Mo": 6, "W": 6,
    "Mn": 2, "Fe": 2, "Co": 2, "Ni": 2, "Cu": 2, "Zn": 2,
    "O": -2, "S": -2, "Se": -2, "Te": -2,
    "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "N": -3, "P": -3,
}


# Bond valence R0 parameters (Angstrom)
BV_R0 = {
    ("Li", "O"): 1.466, ("Na", "O"): 1.803, ("K", "O"): 2.132,
    ("Mg", "O"): 1.693, ("Ca", "O"): 1.967, ("Sr", "O"): 2.118,
    ("Al", "O"): 1.651, ("Si", "O"): 1.624, ("Ti", "O"): 1.815,
    ("Fe", "O"): 1.759, ("Mn", "O"): 1.790, ("Co", "O"): 1.692,
    ("Zn", "O"): 1.704, ("Cu", "O"): 1.679, ("Ni", "O"): 1.654,
}


def check_stoichiometry(
    structure: Dict[str, Any],
    expected_formula: Optional[str] = None
) -> Dict[str, Any]:
    """
    Check structure stoichiometry.
    
    Args:
        structure: Structure dictionary
        expected_formula: Expected formula (optional)
    
    Returns:
        Stoichiometry analysis
    """
    atoms = structure.get("structure", structure).get("atoms", [])
    
    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms"}}
    
    # Count elements
    composition = {}
    for atom in atoms:
        elem = atom.get("element", "X")
        composition[elem] = composition.get(elem, 0) + 1
    
    # Find GCD for reduced formula
    from math import gcd
    from functools import reduce
    
    counts = list(composition.values())
    common = reduce(gcd, counts)
    
    reduced = {k: v // common for k, v in composition.items()}
    
    # Build formula string
    formula = "".join(f"{k}{v if v > 1 else ''}" for k, v in sorted(composition.items()))
    reduced_formula = "".join(f"{k}{v if v > 1 else ''}" for k, v in sorted(reduced.items()))
    
    result = {
        "success": True,
        "composition": composition,
        "formula": formula,
        "reduced_formula": reduced_formula,
        "n_formula_units": common,
        "n_atoms": len(atoms),
    }
    
    if expected_formula:
        result["expected_formula"] = expected_formula
        result["matches_expected"] = reduced_formula == expected_formula or formula == expected_formula
    
    return result


def check_bond_valence(
    structure: Dict[str, Any],
    cutoff: float = 3.0
) -> Dict[str, Any]:
    """
    Check bond valence sums.
    
    Args:
        structure: Structure dictionary
        cutoff: Bond cutoff distance
    
    Returns:
        Bond valence analysis
    """
    atoms = structure.get("structure", structure).get("atoms", [])
    lattice_info = structure.get("structure", structure).get("lattice", {})
    
    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms"}}
    
    # Get lattice parameters for distance calculation
    a = lattice_info.get("a", 10.0)
    b = lattice_info.get("b", a)
    c = lattice_info.get("c", a)
    
    bv_sums = []
    bv_issues = []
    
    for i, atom_i in enumerate(atoms):
        elem_i = atom_i.get("element", "X")
        expected_valence = OXIDATION_STATES.get(elem_i, 0)
        
        bv_sum = 0
        
        for j, atom_j in enumerate(atoms):
            if i == j:
                continue
            
            elem_j = atom_j.get("element", "X")
            
            # Get R0 parameter
            pair = (elem_i, elem_j) if (elem_i, elem_j) in BV_R0 else (elem_j, elem_i)
            r0 = BV_R0.get(pair, None)
            
            if r0 is None:
                continue
            
            # Calculate distance (simplified)
            coord_i = atom_i.get("coords", [0, 0, 0])
            coord_j = atom_j.get("coords", [0, 0, 0])
            
            dx = (coord_i[0] - coord_j[0]) * a
            dy = (coord_i[1] - coord_j[1]) * b
            dz = (coord_i[2] - coord_j[2]) * c
            
            # Periodic boundary
            dx = dx - round(dx / a) * a
            dy = dy - round(dy / b) * b
            dz = dz - round(dz / c) * c
            
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            
            if r < cutoff and r > 0.5:
                bv = np.exp((r0 - r) / 0.37)
                bv_sum += bv
        
        bv_sums.append({"atom_index": i, "element": elem_i, "bv_sum": round(bv_sum, 2), "expected": abs(expected_valence)})
        
        if expected_valence != 0 and abs(bv_sum - abs(expected_valence)) > 0.3:
            bv_issues.append({"atom_index": i, "element": elem_i, "bv_sum": round(bv_sum, 2), "expected": abs(expected_valence)})
    
    return {
        "success": True,
        "bond_valence_sums": bv_sums[:10],  # First 10 for brevity
        "n_atoms_analyzed": len(bv_sums),
        "n_valence_issues": len(bv_issues),
        "issues": bv_issues[:5],
        "is_valid": len(bv_issues) == 0
    }


def check_charge_neutrality(
    structure: Dict[str, Any],
    tolerance: float = 0.1
) -> Dict[str, Any]:
    """
    Check charge neutrality.
    
    Args:
        structure: Structure dictionary
        tolerance: Charge tolerance per formula unit
    
    Returns:
        Charge neutrality analysis
    """
    atoms = structure.get("structure", structure).get("atoms", [])
    
    total_charge = 0
    element_charges = {}
    
    for atom in atoms:
        elem = atom.get("element", "X")
        charge = OXIDATION_STATES.get(elem, 0)
        total_charge += charge
        element_charges[elem] = element_charges.get(elem, 0) + charge
    
    is_neutral = abs(total_charge) <= tolerance * len(atoms) / 10
    
    return {
        "success": True,
        "total_charge": total_charge,
        "charge_per_atom": round(total_charge / len(atoms), 4) if atoms else 0,
        "element_charges": element_charges,
        "is_neutral": is_neutral,
        "tolerance": tolerance
    }


def remove_duplicates(
    structures: List[Dict[str, Any]],
    tolerance: float = 0.01,
    method: str = "composition"
) -> Dict[str, Any]:
    """
    Remove duplicate structures from list.
    
    Args:
        structures: List of structure dictionaries
        tolerance: Comparison tolerance
        method: 'composition', 'structure', 'graph'
    
    Returns:
        List of unique structures
    """
    if not structures:
        return {"success": True, "unique_structures": [], "n_duplicates": 0}
    
    unique = []
    duplicates = 0
    
    for struct in structures:
        is_duplicate = False
        
        for existing in unique:
            if method == "composition":
                # Compare compositions
                atoms1 = struct.get("structure", struct).get("atoms", [])
                atoms2 = existing.get("structure", existing).get("atoms", [])
                
                comp1 = {}
                for a in atoms1:
                    e = a.get("element", "X")
                    comp1[e] = comp1.get(e, 0) + 1
                
                comp2 = {}
                for a in atoms2:
                    e = a.get("element", "X")
                    comp2[e] = comp2.get(e, 0) + 1
                
                is_duplicate = comp1 == comp2
        
        if is_duplicate:
            duplicates += 1
        else:
            unique.append(struct)
    
    return {
        "success": True,
        "n_input": len(structures),
        "n_unique": len(unique),
        "n_duplicates": duplicates,
        "unique_structures": unique
    }
