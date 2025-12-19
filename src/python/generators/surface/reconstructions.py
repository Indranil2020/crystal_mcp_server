"""
surface/reconstructions.py - Surface Reconstructions

Comprehensive surface reconstruction generation per structure_catalogue.md Category 4:
(ii) Classic reconstructions (Si(100)-2×1, Si(111)-7×7, Au(110)-2×1, Pt(110)-1×2)
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# Surface reconstruction database
RECONSTRUCTION_DATABASE = {
    # Silicon reconstructions
    "Si(100)-2x1": {
        "bulk": "Si", "miller": [1, 0, 0], "reconstruction": "dimer_row",
        "unit_cell": [2, 1], "description": "Asymmetric dimer rows"
    },
    "Si(100)-c4x2": {
        "bulk": "Si", "miller": [1, 0, 0], "reconstruction": "buckled_dimer",
        "unit_cell": [4, 2], "description": "Buckled dimer superstructure"
    },
    "Si(100)-4x2": {
        "bulk": "Si", "miller": [1, 0, 0], "reconstruction": "buckled_dimer_ordered",
        "unit_cell": [4, 2], "description": "Ordered buckled dimers"
    },
    "Si(111)-2x1": {
        "bulk": "Si", "miller": [1, 1, 1], "reconstruction": "Pandey_chain",
        "unit_cell": [2, 1], "description": "π-bonded chains"
    },
    "Si(111)-7x7": {
        "bulk": "Si", "miller": [1, 1, 1], "reconstruction": "DAS",
        "unit_cell": [7, 7], "description": "Dimer-adatom-stacking fault",
        "n_adatoms": 12, "n_rest_atoms": 6
    },
    "Si(111)-2x2": {
        "bulk": "Si", "miller": [1, 1, 1], "reconstruction": "adatom",
        "unit_cell": [2, 2], "description": "Adatom overlayer"
    },
    
    # Germanium
    "Ge(100)-2x1": {
        "bulk": "Ge", "miller": [1, 0, 0], "reconstruction": "dimer_row",
        "unit_cell": [2, 1], "description": "Similar to Si(100)"
    },
    "Ge(111)-c2x8": {
        "bulk": "Ge", "miller": [1, 1, 1], "reconstruction": "adatom",
        "unit_cell": [2, 8], "description": "Adatom superstructure"
    },
    
    # Gold
    "Au(110)-2x1": {
        "bulk": "Au", "miller": [1, 1, 0], "reconstruction": "missing_row",
        "unit_cell": [2, 1], "description": "Missing row reconstruction"
    },
    "Au(111)-22x√3": {
        "bulk": "Au", "miller": [1, 1, 1], "reconstruction": "herringbone",
        "unit_cell": [22, 1.73], "description": "Herringbone pattern",
        "stacking_fault": True
    },
    "Au(100)-hex": {
        "bulk": "Au", "miller": [1, 0, 0], "reconstruction": "hexagonal",
        "description": "Quasi-hexagonal overlayer on square substrate"
    },
    
    # Platinum
    "Pt(110)-2x1": {
        "bulk": "Pt", "miller": [1, 1, 0], "reconstruction": "missing_row",
        "unit_cell": [2, 1], "description": "Missing row"
    },
    "Pt(100)-hex": {
        "bulk": "Pt", "miller": [1, 0, 0], "reconstruction": "hexagonal",
        "description": "Hexagonal overlayer"
    },
    
    # III-V semiconductors
    "GaAs(100)-2x4": {
        "bulk": "GaAs", "miller": [1, 0, 0], "reconstruction": "beta2",
        "unit_cell": [2, 4], "description": "β2(2×4) As-rich"
    },
    "GaAs(100)-4x2": {
        "bulk": "GaAs", "miller": [1, 0, 0], "reconstruction": "Ga_rich",
        "unit_cell": [4, 2], "description": "Ga-rich surface"
    },
    "GaAs(110)-1x1": {
        "bulk": "GaAs", "miller": [1, 1, 0], "reconstruction": "relaxation",
        "unit_cell": [1, 1], "description": "Relaxed cleaved surface"
    },
    "InP(100)-2x4": {
        "bulk": "InP", "miller": [1, 0, 0], "reconstruction": "mixed_dimer",
        "unit_cell": [2, 4], "description": "Mixed dimer rows"
    },
    
    # Transition metals
    "W(100)-c2x2": {
        "bulk": "W", "miller": [1, 0, 0], "reconstruction": "lateral_shift",
        "unit_cell": [2, 2], "description": "Lateral displacement"
    },
    "Mo(100)-c2x2": {
        "bulk": "Mo", "miller": [1, 0, 0], "reconstruction": "lateral_shift",
        "unit_cell": [2, 2], "description": "Similar to W(100)"
    },
    "Ir(100)-5x1": {
        "bulk": "Ir", "miller": [1, 0, 0], "reconstruction": "hexagonal",
        "unit_cell": [5, 1], "description": "Quasi-hexagonal"
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


def generate_reconstruction(
    surface: str = "Si(100)-2x1",
    n_layers: int = 6,
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate reconstructed surface.
    
    Args:
        surface: Surface reconstruction from database
        n_layers: Number of bulk layers
        vacuum: Vacuum thickness in Å
    
    Returns:
        Reconstructed surface structure
    """
    if surface not in RECONSTRUCTION_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_SURFACE", "message": f"Unknown surface '{surface}'",
                      "available": list(RECONSTRUCTION_DATABASE.keys())}
        }
    
    info = RECONSTRUCTION_DATABASE[surface]
    element = info["bulk"]
    miller = info["miller"]
    recon_type = info["reconstruction"]
    
    # Get lattice constants
    lattice_params = {
        "Si": 5.43, "Ge": 5.66, "Au": 4.08, "Pt": 3.92, "W": 3.16, "Mo": 3.15, "Ir": 3.84,
        "GaAs": 5.65, "InP": 5.87
    }
    
    a = lattice_params.get(element.replace("As", "").replace("P", ""), 4.0)
    
    # Determine surface cell dimensions
    if miller == [1, 0, 0]:
        a_surf = a / np.sqrt(2)
        interlayer = a / 4
    elif miller == [1, 1, 1]:
        a_surf = a / np.sqrt(2)
        interlayer = a / np.sqrt(3)
    elif miller == [1, 1, 0]:
        a_surf = a
        interlayer = a / (2 * np.sqrt(2))
    else:
        a_surf = a
        interlayer = a / 4
    
    # Get reconstruction unit cell
    unit_cell = info.get("unit_cell", [1, 1])
    nx, ny = int(unit_cell[0]), int(unit_cell[1])
    
    c = n_layers * interlayer + vacuum
    
    lattice = Lattice.orthorhombic(a_surf * nx, a_surf * ny, c)
    
    species = []
    coords = []
    
    # Build bulk layers
    for layer in range(n_layers):
        z_layer = (layer + 0.5) * interlayer / c
        
        for i in range(nx):
            for j in range(ny):
                x = (i + 0.5 * ((layer + 1) % 2)) / nx
                y = (j + 0.5 * ((layer + 1) % 2)) / ny
                
                if "GaAs" in element:
                    species.append("Ga" if layer % 2 == 0 else "As")
                elif "InP" in element:
                    species.append("In" if layer % 2 == 0 else "P")
                else:
                    species.append(element)
                coords.append([x % 1, y % 1, z_layer])
    
    # Apply reconstruction
    z_surface = n_layers * interlayer / c
    
    if recon_type == "dimer_row":
        # Si(100)-2x1: form dimers
        dimer_length = 0.05
        for i in range(nx // 2):
            for j in range(ny):
                x = (2 * i + 0.5) / nx
                y = (j + 0.5) / ny
                
                species.append(element if element in ["Si", "Ge"] else element[0])
                coords.append([x - dimer_length, y, z_surface + 0.01])
                species.append(element if element in ["Si", "Ge"] else element[0])
                coords.append([x + dimer_length, y, z_surface + 0.02])  # Asymmetric
    
    elif recon_type == "missing_row":
        # Au/Pt(110)-2x1: remove every other row
        for i in range(nx):
            if i % 2 == 0:  # Keep every other row
                for j in range(ny):
                    x = (i + 0.5) / nx
                    y = (j + 0.5) / ny
                    species.append(element)
                    coords.append([x, y, z_surface])
    
    elif recon_type == "DAS":
        # Si(111)-7x7: complex DAS structure
        # Simplified: add adatoms and rest atoms
        n_adatoms = info.get("n_adatoms", 12)
        adatom_positions = []
        for i in range(n_adatoms):
            angle = 2 * np.pi * i / n_adatoms
            r = 0.3
            adatom_positions.append([0.5 + r * np.cos(angle), 0.5 + r * np.sin(angle)])
        
        for pos in adatom_positions:
            species.append("Si")
            coords.append([pos[0], pos[1], z_surface + 0.02])
    
    elif recon_type == "hexagonal":
        # Quasi-hexagonal overlayer
        for i in range(nx * 6 // 5):
            for j in range(ny):
                x = (i + (j % 2) * 0.5) / (nx * 6 / 5)
                y = (j * np.sqrt(3) / 2) / ny
                if x < 1 and y < 1:
                    species.append(element)
                    coords.append([x, y, z_surface])
    
    else:
        # Generic reconstruction - just add surface layer
        for i in range(nx):
            for j in range(ny):
                x = (i + 0.5) / nx
                y = (j + 0.5) / ny
                species.append(element if element in ["Si", "Ge", "Au", "Pt", "W", "Mo", "Ir"] else element[0])
                coords.append([x, y, z_surface])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "surface": surface,
        "bulk_material": element,
        "miller_index": miller,
        "reconstruction_type": recon_type,
        "unit_cell": unit_cell,
        "description": info["description"],
        "n_layers": n_layers,
        "vacuum_A": vacuum,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def get_reconstruction_database() -> Dict[str, Any]:
    """Get all available reconstructions."""
    return {
        "success": True,
        "n_reconstructions": len(RECONSTRUCTION_DATABASE),
        "reconstructions": list(RECONSTRUCTION_DATABASE.keys()),
        "by_material": {
            "Si": [k for k in RECONSTRUCTION_DATABASE if "Si" in k],
            "Ge": [k for k in RECONSTRUCTION_DATABASE if "Ge" in k],
            "Au": [k for k in RECONSTRUCTION_DATABASE if "Au" in k],
            "Pt": [k for k in RECONSTRUCTION_DATABASE if "Pt" in k],
            "III-V": [k for k in RECONSTRUCTION_DATABASE if "GaAs" in k or "InP" in k],
        }
    }
