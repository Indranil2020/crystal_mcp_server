"""
defect/point_defects.py - Point Defects

Comprehensive point defect generation per structure_catalogue.md Category 7:
(i) Point defects – vacancies (mono, di, Schottky, Frenkel), interstitials (split, crowdion),
    antisites, impurities (substitutional, interstitial), vacancy-impurity complexes, F-centers
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
import spglib
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# Point defect formation energies (approximate, eV)
DEFECT_FORMATION_ENERGIES = {
    "Si": {"vacancy": 3.6, "interstitial": 3.8, "divacancy": 5.0},
    "Ge": {"vacancy": 2.5, "interstitial": 3.0, "divacancy": 3.8},
    "GaAs": {"V_Ga": 3.5, "V_As": 3.0, "Ga_As": 2.0, "As_Ga": 2.5},
    "NaCl": {"V_Na": 1.5, "V_Cl": 1.8, "Schottky": 2.0, "Frenkel": 2.5},
    "MgO": {"V_Mg": 7.5, "V_O": 6.0, "Schottky": 6.5},
    "Al": {"vacancy": 0.67, "interstitial": 3.0, "divacancy": 0.9},
    "Cu": {"vacancy": 1.28, "interstitial": 3.0, "divacancy": 1.8},
    "Fe_bcc": {"vacancy": 2.0, "interstitial": 4.0},
}


# Common dopants/impurities
DOPANT_DATABASE = {
    "Si": {
        "n_type": ["P", "As", "Sb"],
        "p_type": ["B", "Al", "Ga", "In"],
        "deep": ["Au", "Fe", "Cu"]
    },
    "GaAs": {
        "n_type": ["Si", "Te", "S"],
        "p_type": ["Zn", "Be", "Mg"],
        "isoelectronic": ["N", "P"]
    },
    "SiC": {
        "n_type": ["N", "P"],
        "p_type": ["Al", "B"]
    },
    "ZnO": {
        "n_type": ["Al", "Ga", "In"],
        "p_type": ["N", "Li", "Na"]
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    
    # Get symmetry info via spglib
    cell = (
        structure.lattice.matrix,
        structure.frac_coords,
        structure.atomic_numbers
    )
    
    dataset = spglib.get_symmetry_dataset(cell, symprec=0.1)
    
    if dataset:
        sg_symbol = str(dataset["international"])
        sg_number = int(dataset["number"])
        hall = str(dataset["hall"])
        point_group = str(dataset["pointgroup"])
        
        if 1 <= sg_number <= 2: crystal_system = "triclinic"
        elif 3 <= sg_number <= 15: crystal_system = "monoclinic"
        elif 16 <= sg_number <= 74: crystal_system = "orthorhombic"
        elif 75 <= sg_number <= 142: crystal_system = "tetragonal"
        elif 143 <= sg_number <= 167: crystal_system = "trigonal"
        elif 168 <= sg_number <= 194: crystal_system = "hexagonal"
        elif 195 <= sg_number <= 230: crystal_system = "cubic"
        else: crystal_system = "unknown"
    else:
        # Fallback
        sg_symbol = "P1"
        sg_number = 1
        hall = "P 1"
        point_group = "1"
        crystal_system = "triclinic"

    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "space_group": {
            "number": sg_number,
            "symbol": sg_symbol,
            "hall_symbol": hall,
            "point_group": point_group,
            "crystal_system": crystal_system
        },
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_vacancy(
    host_structure: Union[Structure, Dict],
    vacancy_site: int = 0,
    vacancy_type: str = "single"
) -> Dict[str, Any]:
    """
    Generate vacancy defect.
    
    Args:
        host_structure: Host structure or dict
        vacancy_site: Site index for vacancy
        vacancy_type: 'single', 'divacancy', 'Schottky', 'Frenkel'
    
    Returns:
        Structure with vacancy
    """
    if isinstance(host_structure, dict):
        atoms = list(host_structure.get("atoms", []))
        lattice_info = host_structure.get("lattice", {})
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b, 
                        "c": host_structure.lattice.c}
    
    if vacancy_site >= len(atoms):
        return {
            "success": False,
            "error": {"code": "INVALID_SITE", "message": f"Site {vacancy_site} out of range"}
        }
    
    removed_atom = atoms[vacancy_site]["element"]
    removed_position = atoms[vacancy_site]["coords"]
    
    if vacancy_type == "single":
        del atoms[vacancy_site]
        n_vacancies = 1
        
    elif vacancy_type == "divacancy":
        # Remove two nearest neighbors
        del atoms[vacancy_site]
        if len(atoms) > 0:
            # Find nearest neighbor and remove it too
            nearest_idx = 0
            min_dist = float('inf')
            for i, atom in enumerate(atoms):
                dist = np.sqrt(sum((atom["coords"][j] - removed_position[j])**2 for j in range(3)))
                if dist < min_dist:
                    min_dist = dist
                    nearest_idx = i
            del atoms[nearest_idx]
        n_vacancies = 2
        
    elif vacancy_type == "Schottky":
        # Remove one cation and one anion
        del atoms[vacancy_site]
        # Find an atom of different type to remove
        other_type = None
        other_idx = None
        for i, atom in enumerate(atoms):
            if atom["element"] != removed_atom:
                other_type = atom["element"]
                other_idx = i
                break
        if other_idx is not None:
            del atoms[other_idx]
        n_vacancies = 2
        
    elif vacancy_type == "Frenkel":
        # Move atom to interstitial position (creates vacancy + interstitial)
        atoms[vacancy_site]["coords"] = [
            atoms[vacancy_site]["coords"][0] + 0.25,
            atoms[vacancy_site]["coords"][1] + 0.25,
            atoms[vacancy_site]["coords"][2] + 0.25
        ]
        atoms[vacancy_site]["is_interstitial"] = True
        n_vacancies = 1  # One vacancy created, atom moved
    
    else:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown vacancy type",
                      "available": ["single", "divacancy", "Schottky", "Frenkel"]}
        }
    
    return {
        "success": True,
        "defect_type": "vacancy",
        "vacancy_type": vacancy_type,
        "removed_atom": removed_atom,
        "vacancy_position": removed_position,
        "n_vacancies": n_vacancies,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "sites": atoms, "lattice": lattice_info, "space_group": {"number": 1, "symbol": "P1", "hall_symbol": "P 1", "point_group": "1", "crystal_system": "triclinic"}}
    }


def generate_interstitial(
    host_structure: Union[Structure, Dict],
    interstitial_species: str = "Si",
    interstitial_type: str = "tetrahedral",
    position: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Generate interstitial defect.
    
    Args:
        host_structure: Host structure
        interstitial_species: Element to add
        interstitial_type: 'tetrahedral', 'octahedral', 'split', 'crowdion', 'hex'
        position: Fractional coordinates (optional)
    
    Returns:
        Structure with interstitial
    """
    if isinstance(host_structure, dict):
        atoms = list(host_structure.get("atoms", []))
        lattice_info = host_structure.get("lattice", {})
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b, 
                        "c": host_structure.lattice.c}
    
    # Determine interstitial position
    interstitial_positions = {
        "tetrahedral": [0.25, 0.25, 0.25],
        "octahedral": [0.5, 0.5, 0.5],
        "split": [0.125, 0.125, 0.125],  # <110> split
        "crowdion": [0.5, 0.0, 0.0],      # Extended along <100>
        "hex": [0.5, 0.25, 0.125],        # Hexagonal interstitial
    }
    
    if interstitial_type not in interstitial_positions and position is None:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown interstitial type",
                      "available": list(interstitial_positions.keys())}
        }
    
    int_pos = position if position else interstitial_positions[interstitial_type]
    
    if interstitial_type == "split":
        # Split interstitial: add two atoms sharing one site
        atoms.append({
            "element": interstitial_species,
            "coords": [int_pos[0] - 0.05, int_pos[1] - 0.05, int_pos[2]],
            "is_interstitial": True
        })
        atoms.append({
            "element": interstitial_species,
            "coords": [int_pos[0] + 0.05, int_pos[1] + 0.05, int_pos[2]],
            "is_interstitial": True
        })
    elif interstitial_type == "crowdion":
        # Crowdion: extended defect
        for i in range(3):
            atoms.append({
                "element": interstitial_species,
                "coords": [int_pos[0] + i * 0.1, int_pos[1], int_pos[2]],
                "is_interstitial": True
            })
    else:
        atoms.append({
            "element": interstitial_species,
            "coords": int_pos,
            "is_interstitial": True
        })
    
    return {
        "success": True,
        "defect_type": "interstitial",
        "interstitial_type": interstitial_type,
        "interstitial_species": interstitial_species,
        "position": int_pos,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "sites": atoms, "lattice": lattice_info, "space_group": {"number": 1, "symbol": "P1", "hall_symbol": "P 1", "point_group": "1", "crystal_system": "triclinic"}}
    }


def generate_substitution(
    host_structure: Union[Structure, Dict],
    site: int,
    dopant: str,
    dopant_type: str = "substitutional"
) -> Dict[str, Any]:
    """
    Generate substitutional defect (dopant).
    
    Args:
        host_structure: Host structure
        site: Site index to substitute
        dopant: Dopant element
        dopant_type: 'substitutional' or 'antisite'
    
    Returns:
        Structure with substitution
    """
    if isinstance(host_structure, dict):
        atoms = list(host_structure.get("atoms", []))
        lattice_info = host_structure.get("lattice", {})
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b, 
                        "c": host_structure.lattice.c}
    
    if site >= len(atoms):
        return {
            "success": False,
            "error": {"code": "INVALID_SITE", "message": f"Site {site} out of range"}
        }
    
    original_element = atoms[site]["element"]
    atoms[site]["element"] = dopant
    atoms[site]["is_dopant"] = True
    atoms[site]["original_element"] = original_element
    
    return {
        "success": True,
        "defect_type": "substitution",
        "dopant_type": dopant_type,
        "dopant": dopant,
        "replaced_element": original_element,
        "site": site,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": lattice_info, "space_group": {"number": 1, "symbol": "P1", "hall_symbol": "P 1", "point_group": "1", "crystal_system": "triclinic"}}
    }


def generate_antisite(
    host_structure: Union[Structure, Dict],
    site1: int,
    site2: int
) -> Dict[str, Any]:
    """
    Generate antisite defect (swap two atoms).
    
    Args:
        host_structure: Host structure
        site1: First site index
        site2: Second site index
    
    Returns:
        Structure with antisite
    """
    if isinstance(host_structure, dict):
        atoms = list(host_structure.get("atoms", []))
        lattice_info = host_structure.get("lattice", {})
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b, 
                        "c": host_structure.lattice.c}
    
    if site1 >= len(atoms) or site2 >= len(atoms):
        return {
            "success": False,
            "error": {"code": "INVALID_SITE", "message": f"Site out of range"}
        }
    
    # Swap elements
    elem1 = atoms[site1]["element"]
    elem2 = atoms[site2]["element"]
    
    atoms[site1]["element"] = elem2
    atoms[site1]["is_antisite"] = True
    atoms[site2]["element"] = elem1
    atoms[site2]["is_antisite"] = True
    
    return {
        "success": True,
        "defect_type": "antisite",
        "swapped_elements": [elem1, elem2],
        "sites": [site1, site2],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": lattice_info, "space_group": {"number": 1, "symbol": "P1", "hall_symbol": "P 1", "point_group": "1", "crystal_system": "triclinic"}}
    }


def generate_f_center(
    host_structure: Union[Structure, Dict],
    anion_site: int
) -> Dict[str, Any]:
    """
    Generate F-center (color center) by removing anion.
    
    Args:
        host_structure: Host structure (ionic crystal)
        anion_site: Anion site index to remove
    
    Returns:
        Structure with F-center
    """
    result = generate_vacancy(host_structure, anion_site, "single")

    if result["success"]:
        result["defect_type"] = "F_center"
        result["trapped_electron"] = True
        result["optical_active"] = True

    return result


def create_defect(
    structure: Union[Structure, Dict],
    defect_type: str,
    site_index: int = 0,
    dopant: str = None,
    interstitial_species: str = None,
    interstitial_type: str = "tetrahedral"
) -> Dict[str, Any]:
    """
    Create a defect in a structure (general wrapper).

    Args:
        structure: Host structure (pymatgen Structure or dict)
        defect_type: Type of defect ('vacancy', 'interstitial', 'substitution', 'antisite')
        site_index: Site index for vacancy/substitution
        dopant: Dopant element for substitution
        interstitial_species: Element for interstitial
        interstitial_type: Type of interstitial site

    Returns:
        Structure with defect
    """
    # Convert to dict format if pymatgen Structure
    if isinstance(structure, Structure):
        host = {
            "lattice": {
                "a": structure.lattice.a, "b": structure.lattice.b, "c": structure.lattice.c,
                "matrix": structure.lattice.matrix.tolist()
            },
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    else:
        host = structure

    defect_type_lower = defect_type.lower()

    if defect_type_lower == "vacancy":
        return generate_vacancy(host, site_index, "single")
    elif defect_type_lower == "interstitial":
        species = interstitial_species or "Li"
        return generate_interstitial(host, species, interstitial_type)
    elif defect_type_lower == "substitution":
        if not dopant:
            return {"success": False, "error": {"code": "MISSING_DOPANT",
                    "message": "Must provide 'dopant' for substitution"}}
        return generate_substitution(host, site_index, dopant)
    elif defect_type_lower == "antisite":
        # Need second site for antisite
        site2 = (site_index + 1) % len(host.get("atoms", []))
        return generate_antisite(host, site_index, site2)
    elif defect_type_lower == "f_center":
        return generate_f_center(host, site_index)
    else:
        return {
            "success": False,
            "error": {
                "code": "INVALID_DEFECT_TYPE",
                "message": f"Unknown defect type: {defect_type}",
                "available": ["vacancy", "interstitial", "substitution", "antisite", "f_center"]
            }
        }


def generate_frenkel_pair(
    host_structure: Dict[str, Any],
    vacancy_site: int = 0,
    interstitial_position: Optional[List[float]] = None,
    element: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a Frenkel defect pair (vacancy + interstitial of same species).

    A Frenkel defect occurs when an atom leaves its lattice site and moves
    to an interstitial position, creating a vacancy-interstitial pair.
    Common in ionic crystals like AgCl, AgBr, CaF2.

    Args:
        host_structure: Host structure dictionary
        vacancy_site: Index of atom to move to interstitial (creates vacancy)
        interstitial_position: Fractional coords for interstitial (auto if None)
        element: Element type (auto-detected if None)

    Returns:
        Structure with Frenkel pair defect

    Examples:
        >>> result = generate_frenkel_pair(nacl_structure, vacancy_site=0)
        >>> result["defect_type"]
        'frenkel_pair'
    """
    atoms = host_structure.get("atoms", [])
    lattice_info = host_structure.get("lattice", {})

    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}

    if vacancy_site >= len(atoms):
        return {"success": False, "error": {"code": "INVALID_SITE", "message": f"Site {vacancy_site} out of range"}}

    # Get the displaced atom
    displaced_atom = atoms[vacancy_site]
    displaced_element = element or displaced_atom.get("element", "X")
    original_coords = displaced_atom.get("coords", [0, 0, 0])

    # Calculate interstitial position if not provided
    if interstitial_position is None:
        # Place at tetrahedral-like site offset from original
        interstitial_position = [
            (original_coords[0] + 0.25) % 1.0,
            (original_coords[1] + 0.25) % 1.0,
            (original_coords[2] + 0.25) % 1.0
        ]

    # Create new atom list: remove original, add at interstitial
    new_atoms = []
    for i, atom in enumerate(atoms):
        if i != vacancy_site:
            new_atoms.append(atom.copy())

    # Add the interstitial atom
    new_atoms.append({
        "element": displaced_element,
        "coords": interstitial_position
    })

    # Build pymatgen structure
    matrix = lattice_info.get("matrix", [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    lattice = Lattice(matrix)
    species = [a["element"] for a in new_atoms]
    coords = [a["coords"] for a in new_atoms]
    structure = Structure(lattice, species, coords)

    # Calculate separation
    vacancy_cart = np.dot(original_coords, lattice.matrix)
    interstitial_cart = np.dot(interstitial_position, lattice.matrix)
    separation = np.linalg.norm(interstitial_cart - vacancy_cart)

    return {
        "success": True,
        "defect_type": "frenkel_pair",
        "displaced_element": displaced_element,
        "vacancy_site": vacancy_site,
        "interstitial_position": interstitial_position,
        "separation_angstrom": round(float(separation), 3),
        "n_atoms": len(structure),
        "is_charge_neutral": True,
        "formation_energy_estimate_eV": DEFECT_FORMATION_ENERGIES.get("NaCl", {}).get("Frenkel", 2.5),
        "structure": structure_to_dict(structure)
    }


def generate_trivacancy(
    host_structure: Union[Structure, Dict],
    center_site: int = 0,
    trivacancy_type: str = "triangular"
) -> Dict[str, Any]:
    """
    Generate trivacancy defect (cluster of 3 adjacent vacancies).

    Trivacancies are important in radiation damage and diffusion studies.
    They can form different configurations: triangular (in-plane),
    linear, or mixed orientations.

    Args:
        host_structure: Host structure dictionary or pymatgen Structure
        center_site: Index of center vacancy site
        trivacancy_type: Configuration type:
            - 'triangular': Three vacancies in triangular arrangement
            - 'linear': Three vacancies in a line
            - 'bent': L-shaped arrangement

    Returns:
        Structure with trivacancy

    Examples:
        >>> result = generate_trivacancy(si_structure, center_site=0)
        >>> result["n_vacancies"]
        3
    """
    if isinstance(host_structure, dict):
        atoms = [a.copy() for a in host_structure.get("atoms", [])]
        lattice_info = host_structure.get("lattice", {})
        matrix = lattice_info.get("matrix", [[5, 0, 0], [0, 5, 0], [0, 0, 5]])
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b,
                        "c": host_structure.lattice.c, "matrix": host_structure.lattice.matrix.tolist()}
        matrix = host_structure.lattice.matrix

    if center_site >= len(atoms):
        return {"success": False, "error": {"code": "INVALID_SITE", "message": f"Site {center_site} out of range"}}

    if len(atoms) < 3:
        return {"success": False, "error": {"code": "TOO_FEW_ATOMS", "message": "Need at least 3 atoms for trivacancy"}}

    # Get center position
    center_coords = np.array(atoms[center_site]["coords"])
    center_elem = atoms[center_site]["element"]

    # Find two nearest neighbors to form trivacancy
    distances = []
    for i, atom in enumerate(atoms):
        if i != center_site:
            # Calculate distance in fractional coordinates (simplified)
            delta = np.array(atom["coords"]) - center_coords
            # Apply minimum image convention
            delta = delta - np.round(delta)
            dist = np.linalg.norm(np.dot(delta, matrix))
            distances.append((i, dist, delta))

    distances.sort(key=lambda x: x[1])

    # Select two neighbors based on trivacancy type
    if trivacancy_type == "triangular":
        # First two nearest neighbors
        neighbor1_idx = distances[0][0]
        neighbor2_idx = distances[1][0]
    elif trivacancy_type == "linear":
        # First neighbor and the one most opposite to it
        neighbor1_idx = distances[0][0]
        dir1 = distances[0][2] / np.linalg.norm(distances[0][2]) if np.linalg.norm(distances[0][2]) > 0 else np.array([1, 0, 0])
        best_idx = distances[1][0]
        best_dot = 1.0
        for i, dist, delta in distances[1:min(6, len(distances))]:
            dir2 = delta / np.linalg.norm(delta) if np.linalg.norm(delta) > 0 else np.array([1, 0, 0])
            dot = np.dot(dir1, dir2)
            if dot < best_dot:  # Most antiparallel
                best_dot = dot
                best_idx = i
        neighbor2_idx = best_idx
    elif trivacancy_type == "bent":
        # First neighbor and one perpendicular
        neighbor1_idx = distances[0][0]
        dir1 = distances[0][2] / np.linalg.norm(distances[0][2]) if np.linalg.norm(distances[0][2]) > 0 else np.array([1, 0, 0])
        best_idx = distances[1][0]
        best_perp = 1.0
        for i, dist, delta in distances[1:min(6, len(distances))]:
            dir2 = delta / np.linalg.norm(delta) if np.linalg.norm(delta) > 0 else np.array([1, 0, 0])
            perp = abs(np.dot(dir1, dir2))
            if perp < best_perp:  # Most perpendicular
                best_perp = perp
                best_idx = i
        neighbor2_idx = best_idx
    else:
        return {"success": False, "error": {"code": "INVALID_TYPE", "message": f"Unknown trivacancy type",
                "available": ["triangular", "linear", "bent"]}}

    # Remove vacancies (in reverse order of index to avoid shifting issues)
    vacancy_indices = sorted([center_site, neighbor1_idx, neighbor2_idx], reverse=True)
    vacancy_positions = [atoms[i]["coords"] for i in sorted([center_site, neighbor1_idx, neighbor2_idx])]
    vacancy_elements = [atoms[i]["element"] for i in sorted([center_site, neighbor1_idx, neighbor2_idx])]

    for idx in vacancy_indices:
        del atoms[idx]

    # Build result structure
    if atoms:
        lattice = Lattice(matrix)
        species = [a["element"] for a in atoms]
        coords = [a["coords"] for a in atoms]
        structure = Structure(lattice, species, coords)
    else:
        structure = None

    return {
        "success": True,
        "defect_type": "trivacancy",
        "trivacancy_type": trivacancy_type,
        "n_vacancies": 3,
        "vacancy_elements": vacancy_elements,
        "vacancy_positions": vacancy_positions,
        "formation_energy_estimate_eV": 3 * DEFECT_FORMATION_ENERGIES.get("Si", {}).get("vacancy", 3.6) * 0.85,  # ~15% binding
        "n_atoms": len(atoms),
        "structure": structure_to_dict(structure) if structure else {"atoms": atoms, "lattice": lattice_info}
    }


def generate_defect_gradient(
    host_structure: Union[Structure, Dict],
    defect_element: str = "N",
    gradient_direction: List[float] = None,
    min_concentration: float = 0.0,
    max_concentration: float = 0.1,
    gradient_profile: str = "linear",
    host_element: str = None,
    seed: int = None
) -> Dict[str, Any]:
    """
    Generate structure with spatial gradient of defect concentration.

    Defect gradients are important for:
    - Functionally graded materials
    - Diffusion profile studies
    - Interface engineering
    - Controlled property variation

    Args:
        host_structure: Host structure dictionary or pymatgen Structure
        defect_element: Element to substitute as defect
        gradient_direction: Direction of concentration gradient [x, y, z]
        min_concentration: Minimum defect fraction at gradient start (0-1)
        max_concentration: Maximum defect fraction at gradient end (0-1)
        gradient_profile: 'linear', 'exponential', 'gaussian', 'step', 'sigmoid'
        host_element: Element to replace (if None, replaces majority element)
        seed: Random seed for reproducibility

    Returns:
        Structure with defect gradient

    Examples:
        >>> result = generate_defect_gradient(gaas, defect_element='N', max_concentration=0.1)
        >>> result["gradient_profile"]
        'linear'
    """
    if seed is not None:
        np.random.seed(seed)

    if gradient_direction is None:
        gradient_direction = [0, 0, 1]  # Default: gradient along z

    gradient_dir = np.array(gradient_direction, dtype=float)
    gradient_dir = gradient_dir / np.linalg.norm(gradient_dir)

    if isinstance(host_structure, dict):
        atoms = [a.copy() for a in host_structure.get("atoms", [])]
        lattice_info = host_structure.get("lattice", {})
        matrix = lattice_info.get("matrix", [[5, 0, 0], [0, 5, 0], [0, 0, 5]])
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b,
                        "c": host_structure.lattice.c, "matrix": host_structure.lattice.matrix.tolist()}
        matrix = host_structure.lattice.matrix

    # Determine host element to replace
    if host_element is None:
        # Find majority element
        elem_counts = {}
        for atom in atoms:
            elem = atom["element"]
            elem_counts[elem] = elem_counts.get(elem, 0) + 1
        host_element = max(elem_counts.keys(), key=lambda k: elem_counts[k])

    # Get candidate sites (sites with host element)
    candidate_indices = [i for i, atom in enumerate(atoms) if atom["element"] == host_element]

    if not candidate_indices:
        return {"success": False, "error": {"code": "NO_HOST_SITES",
                "message": f"No sites with element {host_element} found"}}

    # Calculate position along gradient for each candidate
    positions = []
    for idx in candidate_indices:
        coords = np.array(atoms[idx]["coords"])
        cart_pos = np.dot(coords, matrix)
        # Project onto gradient direction
        pos_along_grad = np.dot(cart_pos, gradient_dir)
        positions.append((idx, pos_along_grad))

    # Normalize positions to [0, 1]
    pos_values = [p[1] for p in positions]
    pos_min, pos_max = min(pos_values), max(pos_values)
    if pos_max > pos_min:
        positions = [(idx, (p - pos_min) / (pos_max - pos_min)) for idx, p in positions]
    else:
        positions = [(idx, 0.5) for idx, p in positions]

    # Calculate local concentration based on profile
    def get_concentration(normalized_pos):
        if gradient_profile == "linear":
            return min_concentration + (max_concentration - min_concentration) * normalized_pos
        elif gradient_profile == "exponential":
            return min_concentration + (max_concentration - min_concentration) * (np.exp(normalized_pos) - 1) / (np.e - 1)
        elif gradient_profile == "gaussian":
            # Peak in middle
            return max_concentration * np.exp(-((normalized_pos - 0.5) ** 2) / 0.1) + min_concentration
        elif gradient_profile == "step":
            return max_concentration if normalized_pos > 0.5 else min_concentration
        elif gradient_profile == "sigmoid":
            x = 10 * (normalized_pos - 0.5)  # Scale for sharp transition
            sigmoid = 1 / (1 + np.exp(-x))
            return min_concentration + (max_concentration - min_concentration) * sigmoid
        else:
            return min_concentration + (max_concentration - min_concentration) * normalized_pos

    # Apply substitutions based on local concentration
    n_substituted = 0
    substitution_positions = []

    for idx, norm_pos in positions:
        local_conc = get_concentration(norm_pos)
        if np.random.random() < local_conc:
            atoms[idx]["element"] = defect_element
            atoms[idx]["is_defect"] = True
            n_substituted += 1
            substitution_positions.append(atoms[idx]["coords"])

    # Build result structure
    lattice = Lattice(matrix)
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)

    actual_concentration = n_substituted / len(candidate_indices) if candidate_indices else 0

    return {
        "success": True,
        "defect_type": "concentration_gradient",
        "gradient_profile": gradient_profile,
        "gradient_direction": gradient_direction,
        "defect_element": defect_element,
        "host_element": host_element,
        "min_concentration": min_concentration,
        "max_concentration": max_concentration,
        "actual_mean_concentration": round(actual_concentration, 4),
        "n_substituted": n_substituted,
        "n_total_host_sites": len(candidate_indices),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_vacancy_cluster(
    host_structure: Union[Structure, Dict],
    cluster_size: int = 4,
    cluster_center: List[float] = None,
    cluster_type: str = "compact",
    seed: int = None
) -> Dict[str, Any]:
    """
    Generate a cluster of vacancies.

    Vacancy clusters form during radiation damage and high-temperature
    processing. They can evolve into voids or dislocation loops.

    Args:
        host_structure: Host structure
        cluster_size: Number of vacancies (2-10)
        cluster_center: Center of cluster in fractional coords
        cluster_type: 'compact' (spherical), 'planar' (disk), 'linear'
        seed: Random seed

    Returns:
        Structure with vacancy cluster
    """
    if seed is not None:
        np.random.seed(seed)

    if cluster_size < 2 or cluster_size > 10:
        return {"success": False, "error": {"code": "INVALID_SIZE",
                "message": "Cluster size must be 2-10"}}

    if isinstance(host_structure, dict):
        atoms = [a.copy() for a in host_structure.get("atoms", [])]
        lattice_info = host_structure.get("lattice", {})
        matrix = np.array(lattice_info.get("matrix", [[5, 0, 0], [0, 5, 0], [0, 0, 5]]))
    else:
        atoms = [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in host_structure]
        lattice_info = {"a": host_structure.lattice.a, "b": host_structure.lattice.b,
                        "c": host_structure.lattice.c, "matrix": host_structure.lattice.matrix.tolist()}
        matrix = host_structure.lattice.matrix

    if len(atoms) < cluster_size:
        return {"success": False, "error": {"code": "TOO_FEW_ATOMS",
                "message": f"Need at least {cluster_size} atoms"}}

    # Determine cluster center
    if cluster_center is None:
        cluster_center = [0.5, 0.5, 0.5]
    center = np.array(cluster_center)

    # Calculate distances from center
    distances = []
    for i, atom in enumerate(atoms):
        delta = np.array(atom["coords"]) - center
        delta = delta - np.round(delta)  # Minimum image
        cart_delta = np.dot(delta, matrix)
        dist = np.linalg.norm(cart_delta)
        distances.append((i, dist, cart_delta))

    distances.sort(key=lambda x: x[1])

    # Select atoms for cluster based on type
    if cluster_type == "compact":
        # Simply take nearest atoms
        vacancy_indices = [d[0] for d in distances[:cluster_size]]
    elif cluster_type == "planar":
        # Select atoms close to a plane through center
        # Filter by distance to z=center[2] plane
        plane_distances = []
        for i, dist, delta in distances[:min(3*cluster_size, len(distances))]:
            plane_dist = abs(delta[2])  # Distance to xy plane through center
            plane_distances.append((i, dist, plane_dist))
        plane_distances.sort(key=lambda x: (x[2], x[1]))  # Sort by plane distance, then total distance
        vacancy_indices = [d[0] for d in plane_distances[:cluster_size]]
    elif cluster_type == "linear":
        # Select atoms along a line through center
        if distances:
            direction = distances[0][2] / np.linalg.norm(distances[0][2]) if np.linalg.norm(distances[0][2]) > 0 else np.array([1, 0, 0])
        else:
            direction = np.array([1, 0, 0])
        line_distances = []
        for i, dist, delta in distances[:min(3*cluster_size, len(distances))]:
            proj = np.dot(delta, direction)
            perp_dist = np.linalg.norm(delta - proj * direction)
            line_distances.append((i, dist, perp_dist))
        line_distances.sort(key=lambda x: (x[2], x[1]))
        vacancy_indices = [d[0] for d in line_distances[:cluster_size]]
    else:
        return {"success": False, "error": {"code": "INVALID_TYPE",
                "message": f"Unknown cluster type", "available": ["compact", "planar", "linear"]}}

    # Remove vacancies
    vacancy_positions = [atoms[i]["coords"] for i in sorted(vacancy_indices)]
    vacancy_elements = [atoms[i]["element"] for i in sorted(vacancy_indices)]

    for idx in sorted(vacancy_indices, reverse=True):
        del atoms[idx]

    # Build result
    if atoms:
        lattice = Lattice(matrix)
        species = [a["element"] for a in atoms]
        coords = [a["coords"] for a in atoms]
        structure = Structure(lattice, species, coords)
    else:
        structure = None

    return {
        "success": True,
        "defect_type": "vacancy_cluster",
        "cluster_type": cluster_type,
        "cluster_size": cluster_size,
        "cluster_center": cluster_center,
        "n_vacancies": len(vacancy_indices),
        "vacancy_elements": vacancy_elements,
        "vacancy_positions": vacancy_positions,
        "n_atoms": len(atoms),
        "structure": structure_to_dict(structure) if structure else {"atoms": atoms, "lattice": lattice_info}
    }
