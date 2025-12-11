"""
surface/interfaces.py - Polar and Oxide Interfaces

Generates polar interfaces and complex oxide heterostructures:
- LaAlO3/SrTiO3 (LAO/STO) with polar catastrophe
- NdGaO3/SrTiO3
- LaTiO3/SrTiO3
- Various terminations and reconstructions
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Perovskite oxide database for interfaces
PEROVSKITE_DATABASE = {
    "SrTiO3": {"a": 3.905, "A": "Sr", "B": "Ti", "polar": False, "bandgap": 3.2},
    "LaAlO3": {"a": 3.79, "A": "La", "B": "Al", "polar": True, "bandgap": 5.6},
    "NdGaO3": {"a": 3.86, "A": "Nd", "B": "Ga", "polar": True, "bandgap": 4.2},
    "LaTiO3": {"a": 3.96, "A": "La", "B": "Ti", "polar": True, "bandgap": 0.2},
    "LaVO3": {"a": 3.93, "A": "La", "B": "V", "polar": True, "bandgap": 1.1},
    "SrVO3": {"a": 3.84, "A": "Sr", "B": "V", "polar": False, "bandgap": 0},
    "KTaO3": {"a": 3.99, "A": "K", "B": "Ta", "polar": False, "bandgap": 3.6},
    "BaTiO3": {"a": 4.00, "A": "Ba", "B": "Ti", "polar": False, "bandgap": 3.2},
    "PbTiO3": {"a": 3.90, "A": "Pb", "B": "Ti", "polar": True, "bandgap": 3.4},
    "SrRuO3": {"a": 3.93, "A": "Sr", "B": "Ru", "polar": False, "bandgap": 0},
    "LaNiO3": {"a": 3.84, "A": "La", "B": "Ni", "polar": True, "bandgap": 0},
}


# Pre-defined polar interfaces
POLAR_INTERFACES = {
    "LAO/STO": {
        "top": "LaAlO3", "bottom": "SrTiO3",
        "description": "Prototypical 2DEG interface",
        "2deg_expected": True, "critical_thickness": 4
    },
    "NGO/STO": {
        "top": "NdGaO3", "bottom": "SrTiO3",
        "description": "Alternative 2DEG interface",
        "2deg_expected": True, "critical_thickness": 4
    },
    "LTO/STO": {
        "top": "LaTiO3", "bottom": "SrTiO3",
        "description": "Mott insulator interface",
        "2deg_expected": True, "critical_thickness": 2
    },
    "LAO/KTO": {
        "top": "LaAlO3", "bottom": "KTaO3",
        "description": "High-mobility alternative",
        "2deg_expected": True, "critical_thickness": 4
    },
    "LNO/STO": {
        "top": "LaNiO3", "bottom": "SrTiO3",
        "description": "Metallic oxide interface",
        "2deg_expected": False, "critical_thickness": None
    },
    "SRO/STO": {
        "top": "SrRuO3", "bottom": "SrTiO3",
        "description": "Ferromagnetic metal interface",
        "2deg_expected": False, "critical_thickness": None
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_polar_interface(
    interface: str = "LAO/STO",
    n_layers_top: int = 4,
    n_layers_bottom: int = 6,
    top_termination: str = "AO",
    bottom_termination: str = "TiO2",
    vacuum: float = 15.0,
    compensation: str = "none"
) -> Dict[str, Any]:
    """
    Generate polar oxide interface.
    
    Args:
        interface: Interface type (LAO/STO, NGO/STO, etc.)
        n_layers_top: Unit cells of top layer
        n_layers_bottom: Unit cells of bottom layer
        top_termination: Top surface termination
        bottom_termination: Bottom surface termination (interface)
        vacuum: Vacuum thickness
        compensation: Compensation mechanism ('none', 'oxygen_vacancy', 'La_vacancy', 'cation_intermixing')
    
    Returns:
        Polar interface structure
    """
    if interface not in POLAR_INTERFACES:
        return {
            "success": False,
            "error": {"code": "INVALID_INTERFACE", "message": f"Unknown interface '{interface}'",
                      "available": list(POLAR_INTERFACES.keys())}
        }
    
    int_info = POLAR_INTERFACES[interface]
    top_mat = int_info["top"]
    bottom_mat = int_info["bottom"]
    
    top_info = PEROVSKITE_DATABASE[top_mat]
    bottom_info = PEROVSKITE_DATABASE[bottom_mat]
    
    # Use average lattice constant
    a = (top_info["a"] + bottom_info["a"]) / 2
    
    # Total thickness
    c_perovskite = a  # Cubic perovskite
    total_layers = n_layers_bottom + n_layers_top
    c_total = total_layers * c_perovskite + vacuum
    
    lattice = Lattice.tetragonal(a, c_total)
    
    species = []
    coords = []
    
    # Build bottom (substrate)
    for layer in range(n_layers_bottom):
        z_base = (layer * c_perovskite + vacuum / 2) / c_total
        
        # AO layer (e.g., SrO)
        species.append(bottom_info["A"])
        coords.append([0, 0, z_base])
        species.append("O")
        coords.append([0.5, 0.5, z_base])
        
        # BO2 layer (e.g., TiO2)
        z_bo2 = z_base + 0.5 * c_perovskite / c_total
        species.append(bottom_info["B"])
        coords.append([0.5, 0.5, z_bo2])
        species.append("O")
        coords.append([0.5, 0, z_bo2])
        species.append("O")
        coords.append([0, 0.5, z_bo2])
    
    # Build top (polar film)
    for layer in range(n_layers_top):
        z_base = ((n_layers_bottom + layer) * c_perovskite + vacuum / 2) / c_total
        
        # AO layer (e.g., LaO)
        species.append(top_info["A"])
        coords.append([0, 0, z_base])
        species.append("O")
        coords.append([0.5, 0.5, z_base])
        
        # BO2 layer (e.g., AlO2)
        z_bo2 = z_base + 0.5 * c_perovskite / c_total
        species.append(top_info["B"])
        coords.append([0.5, 0.5, z_bo2])
        species.append("O")
        coords.append([0.5, 0, z_bo2])
        species.append("O")
        coords.append([0, 0.5, z_bo2])
    
    structure = Structure(lattice, species, coords)
    
    # Apply compensation mechanism
    if compensation == "oxygen_vacancy" and len(structure) > 5:
        # Remove one O from top surface
        o_indices = [i for i, s in enumerate(structure) if str(s.specie) == "O" and s.frac_coords[2] > 0.8]
        if o_indices:
            structure.remove_sites([o_indices[-1]])
    
    # Calculate polar discontinuity
    polar_discontinuity = 0.5 if int_info["2deg_expected"] else 0
    
    return {
        "success": True,
        "interface": interface,
        "top_material": top_mat,
        "bottom_material": bottom_mat,
        "n_layers_top": n_layers_top,
        "n_layers_bottom": n_layers_bottom,
        "top_termination": top_termination,
        "bottom_termination": bottom_termination,
        "lattice_mismatch_percent": round((top_info["a"] - bottom_info["a"]) / bottom_info["a"] * 100, 2),
        "polar_discontinuity_e": polar_discontinuity,
        "2deg_expected": int_info["2deg_expected"],
        "critical_thickness_uc": int_info["critical_thickness"],
        "above_critical": n_layers_top >= (int_info["critical_thickness"] or 0),
        "compensation": compensation,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_oxide_superlattice(
    materials: List[str] = ["LaAlO3", "SrTiO3"],
    layer_thicknesses: List[int] = [2, 2],
    n_periods: int = 5,
    vacuum: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate oxide superlattice.
    
    Args:
        materials: List of perovskite materials
        layer_thicknesses: Thickness of each layer in unit cells
        n_periods: Number of superlattice periods
        vacuum: Vacuum for slab (None for bulk)
    
    Returns:
        Oxide superlattice structure
    """
    if not all(m in PEROVSKITE_DATABASE for m in materials):
        return {
            "success": False,
            "error": {"code": "INVALID_MATERIAL", "message": "Unknown material",
                      "available": list(PEROVSKITE_DATABASE.keys())}
        }
    
    # Average lattice constant
    a = np.mean([PEROVSKITE_DATABASE[m]["a"] for m in materials])
    
    c_uc = a
    period_thickness = sum(layer_thicknesses)
    total_uc = period_thickness * n_periods
    
    c_total = total_uc * c_uc + (vacuum if vacuum else 0)
    
    lattice = Lattice.tetragonal(a, c_total)
    
    species = []
    coords = []
    
    uc_counter = 0
    for period in range(n_periods):
        for mat_idx, (mat, thickness) in enumerate(zip(materials, layer_thicknesses)):
            mat_info = PEROVSKITE_DATABASE[mat]
            
            for layer in range(thickness):
                z_base = (uc_counter * c_uc + (vacuum / 2 if vacuum else 0)) / c_total
                
                # AO layer
                species.append(mat_info["A"])
                coords.append([0, 0, z_base])
                species.append("O")
                coords.append([0.5, 0.5, z_base])
                
                # BO2 layer
                z_bo2 = z_base + 0.5 * c_uc / c_total
                species.append(mat_info["B"])
                coords.append([0.5, 0.5, z_bo2])
                species.append("O")
                coords.append([0.5, 0, z_bo2])
                species.append("O")
                coords.append([0, 0.5, z_bo2])
                
                uc_counter += 1
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "materials": materials,
        "layer_thicknesses": layer_thicknesses,
        "period_thickness_uc": period_thickness,
        "n_periods": n_periods,
        "total_thickness_uc": total_uc,
        "is_slab": vacuum is not None,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_interface_with_defects(
    interface: str = "LAO/STO",
    n_layers_top: int = 4,
    n_layers_bottom: int = 6,
    defect_type: str = "oxygen_vacancy",
    defect_concentration: float = 0.1,
    defect_location: str = "interface",
    vacuum: float = 15.0,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate interface with point defects.
    
    Args:
        interface: Interface type
        n_layers_top: Top layer thickness
        n_layers_bottom: Bottom layer thickness
        defect_type: 'oxygen_vacancy', 'cation_vacancy', 'antisite'
        defect_concentration: Fraction of sites with defects
        defect_location: 'interface', 'bulk', 'surface', 'all'
        vacuum: Vacuum thickness
        seed: Random seed
    
    Returns:
        Interface with defects
    """
    np.random.seed(seed)
    
    # First generate pristine interface
    result = generate_polar_interface(interface, n_layers_top, n_layers_bottom, vacuum=vacuum)
    
    if not result["success"]:
        return result
    
    # Get structure
    struct_dict = result["structure"]
    atoms = struct_dict["atoms"]
    lattice_data = struct_dict["lattice"]
    
    # Create structure
    lattice = Lattice.from_parameters(
        lattice_data["a"], lattice_data["b"], lattice_data["c"],
        lattice_data.get("alpha", 90), lattice_data.get("beta", 90), lattice_data.get("gamma", 90)
    )
    species = [a["element"] for a in atoms]
    coords = [a["coords"] for a in atoms]
    structure = Structure(lattice, species, coords)
    
    # Identify defect region
    total_layers = n_layers_top + n_layers_bottom
    interface_z_min = (n_layers_bottom - 1) / total_layers * 0.8
    interface_z_max = (n_layers_bottom + 1) / total_layers * 0.8
    
    # Find sites to modify
    defect_sites = []
    for i, site in enumerate(structure):
        z = site.frac_coords[2]
        
        in_region = False
        if defect_location == "interface":
            in_region = interface_z_min < z < interface_z_max
        elif defect_location == "surface":
            in_region = z > 0.7
        elif defect_location == "bulk":
            in_region = z < interface_z_min
        else:  # all
            in_region = True
        
        if in_region:
            if defect_type == "oxygen_vacancy" and str(site.specie) == "O":
                defect_sites.append(i)
            elif defect_type == "cation_vacancy" and str(site.specie) != "O":
                defect_sites.append(i)
    
    # Remove random subset
    n_defects = max(1, int(len(defect_sites) * defect_concentration))
    np.random.shuffle(defect_sites)
    sites_to_remove = sorted(defect_sites[:n_defects], reverse=True)
    
    for idx in sites_to_remove:
        structure.remove_sites([idx])
    
    result["defect_type"] = defect_type
    result["n_defects"] = n_defects
    result["defect_location"] = defect_location
    result["defect_concentration"] = defect_concentration
    result["n_atoms"] = len(structure)
    result["structure"] = structure_to_dict(structure)
    
    return result
