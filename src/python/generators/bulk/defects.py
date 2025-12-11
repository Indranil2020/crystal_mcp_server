"""
bulk/defects.py - Point Defects in Bulk Crystals

Generates point defect structures:
- Vacancies (mono, di, tri, Schottky, Frenkel)
- Interstitials (octahedral, tetrahedral, split, crowdion)
- Antisites
- Dopant/impurity distributions

Scientific basis:
- Kröger-Vink notation for defect chemistry
- Defect formation energies
- Configurational entropy of defects
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice, Element

from .base import structure_to_dict


# Common interstitial sites in crystal structures
INTERSTITIAL_SITES = {
    "fcc": {
        "octahedral": [[0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
        "tetrahedral": [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
    },
    "bcc": {
        "octahedral": [[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        "tetrahedral": [[0.5, 0.25, 0], [0.5, 0, 0.25], [0, 0.5, 0.25], [0.25, 0.5, 0]]
    }
}


def generate_vacancy(
    structure_dict: Dict[str, Any],
    site_index: Optional[int] = None,
    element: Optional[str] = None,
    n_vacancies: int = 1,
    vacancy_type: str = "mono"
) -> Dict[str, Any]:
    """
    Generate vacancy defect structure.
    
    Args:
        structure_dict: Base structure
        site_index: Specific site to remove (if None, selects by element)
        element: Element to create vacancy for
        n_vacancies: Number of vacancies to create
        vacancy_type: 'mono', 'di' (divacancy), 'schottky', 'frenkel'
    
    Returns:
        Structure with vacancy/vacancies
    """
    atoms = structure_dict.get("atoms", [])
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY_STRUCTURE", "message": "No atoms"}}
    
    # Find sites to remove
    sites_to_remove = []
    
    if site_index is not None:
        if site_index >= len(atoms):
            return {"success": False, "error": {"code": "INVALID_INDEX", "message": f"Index {site_index} >= n_atoms"}}
        sites_to_remove.append(site_index)
    elif element:
        element_sites = [i for i, a in enumerate(atoms) if a["element"] == element]
        if not element_sites:
            return {"success": False, "error": {"code": "ELEMENT_NOT_FOUND", "message": f"No {element} found"}}
        sites_to_remove = element_sites[:n_vacancies]
    else:
        sites_to_remove = list(range(min(n_vacancies, len(atoms))))
    
    # Special vacancy types
    if vacancy_type == "schottky":
        # Remove stoichiometric vacancies (e.g., one cation + one anion)
        from collections import Counter
        elem_counts = Counter(a["element"] for a in atoms)
        if len(elem_counts) >= 2:
            elements = list(elem_counts.keys())
            for elem in elements[:2]:
                elem_sites = [i for i, a in enumerate(atoms) if a["element"] == elem]
                if elem_sites and elem_sites[0] not in sites_to_remove:
                    sites_to_remove.append(elem_sites[0])
    
    elif vacancy_type == "di":
        # Divacancy: nearest neighbor pair
        if len(sites_to_remove) == 1:
            # Find nearest neighbor
            base_idx = sites_to_remove[0]
            base_coords = np.array(atoms[base_idx]["coords"])
            min_dist = float('inf')
            nn_idx = None
            for i, a in enumerate(atoms):
                if i != base_idx:
                    dist = np.linalg.norm(np.array(a["coords"]) - base_coords)
                    if dist < min_dist:
                        min_dist = dist
                        nn_idx = i
            if nn_idx is not None:
                sites_to_remove.append(nn_idx)
    
    # Remove sites
    new_atoms = [a for i, a in enumerate(atoms) if i not in sites_to_remove]
    
    return {
        "success": True,
        "vacancy_type": vacancy_type,
        "removed_sites": sites_to_remove,
        "n_removed": len(sites_to_remove),
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": {"defect_type": "vacancy"}
        }
    }


def generate_interstitial(
    structure_dict: Dict[str, Any],
    element: str,
    site_type: str = "octahedral",
    crystal_type: str = "fcc",
    position: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Generate interstitial defect structure.
    
    Args:
        structure_dict: Base structure
        element: Element to insert
        site_type: 'octahedral', 'tetrahedral', 'split', 'crowdion'
        crystal_type: 'fcc', 'bcc' for automatic site selection
        position: Explicit fractional coordinates (overrides site_type)
    
    Returns:
        Structure with interstitial
    """
    atoms = list(structure_dict.get("atoms", []))
    lattice = structure_dict.get("lattice", {})
    
    if position is not None:
        interstitial_pos = position
    elif crystal_type in INTERSTITIAL_SITES and site_type in INTERSTITIAL_SITES[crystal_type]:
        interstitial_pos = INTERSTITIAL_SITES[crystal_type][site_type][0]
    else:
        # Default: center of cell
        interstitial_pos = [0.5, 0.5, 0.5]
    
    if site_type == "split":
        # Split interstitial: two atoms displaced from lattice site
        split_dist = 0.05
        atoms.append({
            "element": element,
            "coords": [interstitial_pos[0] + split_dist, interstitial_pos[1], interstitial_pos[2]],
            "cartesian": None
        })
        atoms.append({
            "element": element,
            "coords": [interstitial_pos[0] - split_dist, interstitial_pos[1], interstitial_pos[2]],
            "cartesian": None
        })
    elif site_type == "crowdion":
        # Crowdion: displaced along close-packed direction
        atoms.append({
            "element": element,
            "coords": [interstitial_pos[0] + 0.1, interstitial_pos[1] + 0.1, interstitial_pos[2] + 0.1],
            "cartesian": None
        })
    else:
        atoms.append({
            "element": element,
            "coords": interstitial_pos,
            "cartesian": None
        })
    
    return {
        "success": True,
        "interstitial_type": site_type,
        "element": element,
        "position": interstitial_pos,
        "structure": {
            "lattice": lattice,
            "atoms": atoms,
            "metadata": {"defect_type": "interstitial"}
        }
    }


def generate_antisite(
    structure_dict: Dict[str, Any],
    site_a: int,
    site_b: int
) -> Dict[str, Any]:
    """
    Generate antisite defect by swapping two atoms.
    
    Args:
        structure_dict: Base structure
        site_a: First site index
        site_b: Second site index
    
    Returns:
        Structure with antisite
    """
    atoms = structure_dict.get("atoms", [])
    if site_a >= len(atoms) or site_b >= len(atoms):
        return {"success": False, "error": {"code": "INVALID_INDEX", "message": "Site indices out of range"}}
    
    if atoms[site_a]["element"] == atoms[site_b]["element"]:
        return {"success": False, "error": {"code": "SAME_ELEMENT", "message": "Antisite requires different elements"}}
    
    new_atoms = [dict(a) for a in atoms]
    new_atoms[site_a]["element"], new_atoms[site_b]["element"] = \
        new_atoms[site_b]["element"], new_atoms[site_a]["element"]
    
    return {
        "success": True,
        "swapped": {
            "site_a": {"index": site_a, "original": atoms[site_a]["element"], "new": new_atoms[site_a]["element"]},
            "site_b": {"index": site_b, "original": atoms[site_b]["element"], "new": new_atoms[site_b]["element"]}
        },
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": {"defect_type": "antisite"}
        }
    }


def generate_dopant_distribution(
    structure_dict: Dict[str, Any],
    host_element: str,
    dopant_element: str,
    concentration: float,
    distribution: str = "random",
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate a doped structure with specified dopant distribution.
    
    Args:
        structure_dict: Base structure
        host_element: Element to substitute
        dopant_element: Dopant element
        concentration: Dopant concentration (0-1)
        distribution: 'random', 'gradient', 'clustered', 'periodic'
        seed: Random seed for reproducibility
    
    Returns:
        Doped structure
    """
    if seed is not None:
        np.random.seed(seed)
    
    atoms = structure_dict.get("atoms", [])
    host_sites = [(i, a) for i, a in enumerate(atoms) if a["element"] == host_element]
    
    if not host_sites:
        return {"success": False, "error": {"code": "NO_HOST", "message": f"No {host_element} atoms found"}}
    
    n_dopants = int(round(concentration * len(host_sites)))
    
    if distribution == "random":
        dopant_indices = np.random.choice([s[0] for s in host_sites], size=n_dopants, replace=False)
    elif distribution == "gradient":
        # Dopants concentrated at high-z positions
        sorted_sites = sorted(host_sites, key=lambda x: x[1]["coords"][2], reverse=True)
        dopant_indices = [s[0] for s in sorted_sites[:n_dopants]]
    elif distribution == "clustered":
        # Dopants clustered near center
        center = np.array([0.5, 0.5, 0.5])
        sorted_sites = sorted(host_sites, key=lambda x: np.linalg.norm(np.array(x[1]["coords"]) - center))
        dopant_indices = [s[0] for s in sorted_sites[:n_dopants]]
    else:
        dopant_indices = [s[0] for s in host_sites[:n_dopants]]
    
    new_atoms = [dict(a) for a in atoms]
    for idx in dopant_indices:
        new_atoms[idx]["element"] = dopant_element
    
    return {
        "success": True,
        "host": host_element,
        "dopant": dopant_element,
        "concentration": n_dopants / len(host_sites),
        "n_dopants": n_dopants,
        "distribution": distribution,
        "structure": {
            "lattice": structure_dict.get("lattice", {}),
            "atoms": new_atoms,
            "metadata": {"defect_type": "dopant"}
        }
    }
