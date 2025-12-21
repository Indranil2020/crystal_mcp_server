"""
defect/point_defects.py - Point Defects

Comprehensive point defect generation per structure_catalogue.md Category 7:
(i) Point defects – vacancies (mono, di, Schottky, Frenkel), interstitials (split, crowdion),
    antisites, impurities (substitutional, interstitial), vacancy-impurity complexes, F-centers
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
from pymatgen.core import Structure, Lattice


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
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
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
        "structure": {"atoms": atoms, "sites": atoms, "lattice": lattice_info}
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
        "structure": {"atoms": atoms, "sites": atoms, "lattice": lattice_info}
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
        "structure": {"atoms": atoms, "lattice": lattice_info}
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
        "structure": {"atoms": atoms, "lattice": lattice_info}
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
