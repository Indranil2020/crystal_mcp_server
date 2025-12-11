"""
bulk/disorder.py - Disordered Structure Generation

Generates chemically disordered structures:
- Random alloys
- Special Quasirandom Structures (SQS)
- Monte Carlo SQS (MCSQS)
- Short-range order parameters

Scientific basis:
- Cluster expansion theory
- Warren-Cowley short-range order parameters
- ATAT mcsqs algorithm principles
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice

from .base import structure_to_dict


def generate_random_alloy(
    structure_dict: Dict[str, Any],
    substitutions: Dict[str, Tuple[str, float]],
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate random alloy by random substitution.
    
    Args:
        structure_dict: Base structure dictionary
        substitutions: Dict mapping original element to (new element, fraction)
                      e.g., {"Si": ("Ge", 0.5)} for 50% Ge on Si sites
        seed: Random seed for reproducibility
    
    Returns:
        Structure with random substitutions
    
    Examples:
        >>> result = generate_random_alloy(si_struct, {"Si": ("Ge", 0.5)}, seed=42)
        >>> result["composition"]
        {'Si': 4, 'Ge': 4}
    """
    if seed is not None:
        np.random.seed(seed)
    
    atoms = structure_dict.get("atoms", [])
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY_STRUCTURE", "message": "No atoms"}}
    
    new_atoms = []
    substituted_sites = []
    
    for i, atom in enumerate(atoms):
        atom_copy = dict(atom)
        original_elem = atom["element"]
        
        if original_elem in substitutions:
            new_elem, fraction = substitutions[original_elem]
            if np.random.random() < fraction:
                atom_copy["element"] = new_elem
                substituted_sites.append(i)
        
        new_atoms.append(atom_copy)
    
    # Count composition
    from collections import Counter
    composition = Counter(a["element"] for a in new_atoms)
    
    return {
        "success": True,
        "seed": seed,
        "n_substituted": len(substituted_sites),
        "composition": dict(composition),
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": structure_dict.get("metadata", {})
        }
    }


def generate_sqs(
    structure_dict: Dict[str, Any],
    target_composition: Dict[str, float],
    n_iterations: int = 1000,
    correlation_cutoff: float = 6.0,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate Special Quasirandom Structure (SQS) using simulated annealing.
    
    SQS mimics the correlation functions of a random alloy while
    using a smaller supercell. This implementation uses a simplified
    pair correlation matching approach.
    
    Args:
        structure_dict: Base (ordered) structure
        target_composition: Target element fractions, e.g., {"Si": 0.5, "Ge": 0.5}
        n_iterations: Number of MC swap iterations
        correlation_cutoff: Distance cutoff for pair correlations (Å)
        seed: Random seed
    
    Returns:
        SQS structure with minimized correlation error
    """
    if seed is not None:
        np.random.seed(seed)
    
    atoms = structure_dict.get("atoms", [])
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY_STRUCTURE", "message": "No atoms"}}
    
    n_atoms = len(atoms)
    elements = list(target_composition.keys())
    
    # Initialize with target composition
    n_per_element = {el: int(round(frac * n_atoms)) for el, frac in target_composition.items()}
    
    # Adjust for rounding
    total = sum(n_per_element.values())
    if total != n_atoms:
        diff = n_atoms - total
        n_per_element[elements[0]] += diff
    
    # Create initial assignment
    assignment = []
    for el, count in n_per_element.items():
        assignment.extend([el] * count)
    np.random.shuffle(assignment)
    
    # Get positions
    lattice_matrix = np.array(structure_dict["lattice"]["matrix"])
    frac_coords = np.array([a["coords"] for a in atoms])
    cart_coords = frac_coords @ lattice_matrix
    
    def compute_pair_correlations(elem_list):
        """Compute pair correlation function."""
        correlations = {}
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = np.linalg.norm(cart_coords[i] - cart_coords[j])
                if dist < correlation_cutoff:
                    pair = tuple(sorted([elem_list[i], elem_list[j]]))
                    shell = int(dist / 1.0)  # 1 Å bins
                    key = (pair, shell)
                    correlations[key] = correlations.get(key, 0) + 1
        return correlations
    
    def correlation_objective(elem_list):
        """Objective: deviation from random alloy correlations."""
        corr = compute_pair_correlations(elem_list)
        # Target: uncorrelated (50/50 for each shell)
        obj = 0
        for key, count in corr.items():
            pair, shell = key
            if pair[0] == pair[1]:
                # Same element pairs
                expected = target_composition[pair[0]] ** 2
            else:
                expected = 2 * target_composition[pair[0]] * target_composition[pair[1]]
            obj += (count - expected * n_atoms) ** 2  # Simplified
        return obj
    
    # Simulated annealing
    current = list(assignment)
    current_obj = correlation_objective(current)
    best = current.copy()
    best_obj = current_obj
    
    temperature = 1.0
    cooling_rate = 0.995
    
    for iteration in range(n_iterations):
        # Propose swap
        i, j = np.random.choice(n_atoms, 2, replace=False)
        if current[i] != current[j]:
            # Swap
            current[i], current[j] = current[j], current[i]
            new_obj = correlation_objective(current)
            
            # Accept/reject
            delta = new_obj - current_obj
            if delta < 0 or np.random.random() < np.exp(-delta / temperature):
                current_obj = new_obj
                if new_obj < best_obj:
                    best = current.copy()
                    best_obj = new_obj
            else:
                # Reject - swap back
                current[i], current[j] = current[j], current[i]
        
        temperature *= cooling_rate
    
    # Build final structure
    new_atoms = []
    for i, atom in enumerate(atoms):
        atom_copy = dict(atom)
        atom_copy["element"] = best[i]
        new_atoms.append(atom_copy)
    
    from collections import Counter
    composition = Counter(best)
    
    return {
        "success": True,
        "method": "SQS (simulated annealing)",
        "n_iterations": n_iterations,
        "correlation_error": float(best_obj),
        "composition": dict(composition),
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": structure_dict.get("metadata", {})
        }
    }


def calculate_short_range_order(
    structure_dict: Dict[str, Any],
    element_a: str,
    element_b: str,
    cutoff_radii: List[float] = [3.0, 4.5, 6.0]
) -> Dict[str, Any]:
    """
    Calculate Warren-Cowley short-range order parameters.
    
    α = 1 - P(B|A) / c_B
    where P(B|A) is probability of B neighbor given A central atom
    and c_B is overall concentration of B.
    
    α = 0: random
    α > 0: A-A and B-B clustering (phase separation tendency)
    α < 0: A-B ordering (compound formation tendency)
    
    Args:
        structure_dict: Structure to analyze
        element_a: First element
        element_b: Second element
        cutoff_radii: List of coordination shell radii
    
    Returns:
        SRO parameters for each shell
    """
    atoms = structure_dict.get("atoms", [])
    lattice_matrix = np.array(structure_dict["lattice"]["matrix"])
    
    # Get positions and elements
    positions = []
    elements = []
    for atom in atoms:
        positions.append(np.array(atom["coords"]) @ lattice_matrix)
        elements.append(atom["element"])
    
    positions = np.array(positions)
    
    # Calculate concentrations
    c_a = sum(1 for e in elements if e == element_a) / len(elements)
    c_b = sum(1 for e in elements if e == element_b) / len(elements)
    
    if c_a == 0 or c_b == 0:
        return {"success": False, "error": {"code": "MISSING_ELEMENT", "message": f"Both {element_a} and {element_b} must be present"}}
    
    sro_params = []
    prev_cutoff = 0
    
    for cutoff in cutoff_radii:
        n_ab = 0  # A-B pairs
        n_a_neighbors = 0  # Total neighbors of A atoms
        
        for i, (pos_i, elem_i) in enumerate(zip(positions, elements)):
            if elem_i == element_a:
                for j, (pos_j, elem_j) in enumerate(zip(positions, elements)):
                    if i != j:
                        dist = np.linalg.norm(pos_i - pos_j)
                        if prev_cutoff < dist <= cutoff:
                            n_a_neighbors += 1
                            if elem_j == element_b:
                                n_ab += 1
        
        if n_a_neighbors > 0:
            p_b_given_a = n_ab / n_a_neighbors
            alpha = 1 - p_b_given_a / c_b
        else:
            alpha = 0
        
        sro_params.append({
            "shell_max": cutoff,
            "shell_min": prev_cutoff,
            "alpha": float(alpha),
            "interpretation": "ordered" if alpha < -0.1 else "clustered" if alpha > 0.1 else "random"
        })
        
        prev_cutoff = cutoff
    
    return {
        "success": True,
        "element_a": element_a,
        "element_b": element_b,
        "concentrations": {"a": c_a, "b": c_b},
        "shells": sro_params
    }
