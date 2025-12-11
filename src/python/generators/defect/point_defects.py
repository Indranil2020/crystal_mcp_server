"""
defect/point_defects.py - Point Defect Generation
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


def generate_vacancy(
    structure_dict: Dict[str, Any],
    site_index: int = 0,
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """Generate vacancy defect."""
    atoms = list(structure_dict.get("atoms", []))
    lattice_data = structure_dict.get("lattice", {})
    
    if not atoms:
        return {"success": False, "error": {"code": "EMPTY", "message": "Empty structure"}}
    
    # Build supercell first
    if "matrix" in lattice_data:
        lattice = Lattice(lattice_data["matrix"])
    else:
        lattice = Lattice.from_parameters(
            lattice_data.get("a", 5), lattice_data.get("b", 5), lattice_data.get("c", 5), 90, 90, 90
        )
    
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)
    structure.make_supercell(supercell)
    
    # Remove atom
    if site_index < len(structure):
        removed_element = str(structure[site_index].specie)
        structure.remove_sites([site_index])
    else:
        removed_element = "?"
    
    return {
        "success": True,
        "defect_type": "vacancy",
        "removed_element": removed_element,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }


def generate_interstitial(
    structure_dict: Dict[str, Any],
    element: str = "H",
    position: List[float] = [0.5, 0.5, 0.5],
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """Generate interstitial defect."""
    atoms = list(structure_dict.get("atoms", []))
    lattice_data = structure_dict.get("lattice", {})
    
    if "matrix" in lattice_data:
        lattice = Lattice(lattice_data["matrix"])
    else:
        lattice = Lattice.from_parameters(
            lattice_data.get("a", 5), lattice_data.get("b", 5), lattice_data.get("c", 5), 90, 90, 90
        )
    
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)
    structure.make_supercell(supercell)
    
    # Scale position for supercell
    scaled_pos = [p / s for p, s in zip(position, supercell)]
    structure.append(element, scaled_pos)
    
    return {
        "success": True,
        "defect_type": "interstitial",
        "added_element": element,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }


def generate_substitution(
    structure_dict: Dict[str, Any],
    site_index: int = 0,
    new_element: str = "Fe",
    supercell: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """Generate substitutional defect."""
    atoms = list(structure_dict.get("atoms", []))
    lattice_data = structure_dict.get("lattice", {})
    
    if "matrix" in lattice_data:
        lattice = Lattice(lattice_data["matrix"])
    else:
        lattice = Lattice.from_parameters(
            lattice_data.get("a", 5), lattice_data.get("b", 5), lattice_data.get("c", 5), 90, 90, 90
        )
    
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)
    structure.make_supercell(supercell)
    
    if site_index < len(structure):
        old_element = str(structure[site_index].specie)
        structure[site_index] = new_element
    else:
        old_element = "?"
    
    return {
        "success": True,
        "defect_type": "substitution",
        "old_element": old_element,
        "new_element": new_element,
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }
