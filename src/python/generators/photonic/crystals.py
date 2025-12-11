"""
photonic/crystals.py - Photonic Crystal Generation
"""

from typing import Dict, Any, List
import numpy as np
from pymatgen.core import Structure, Lattice


PHOTONIC_LATTICES = {
    "fcc": {"description": "Face-centered cubic opal"},
    "diamond": {"description": "Diamond lattice (inverse opal)"},
    "woodpile": {"description": "Layer-by-layer woodpile"},
    "gyroid": {"description": "Gyroid minimal surface"},
}


def generate_photonic_crystal(
    lattice_type: str = "fcc",
    period: float = 500.0,  # nm
    fill_fraction: float = 0.74,
    material: str = "Si",
    size: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate photonic crystal structure.
    
    Args:
        lattice_type: Crystal lattice type
        period: Lattice period (nm)
        fill_fraction: Volume fill fraction
        material: Material
        size: Supercell size
    
    Returns:
        Photonic crystal structure
    """
    if lattice_type not in PHOTONIC_LATTICES:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown lattice"}}
    
    a = period
    
    if lattice_type == "fcc":
        lattice = Lattice.cubic(a)
        species = [material] * 4
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        
    elif lattice_type == "diamond":
        lattice = Lattice.cubic(a)
        species = [material] * 8
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                  [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
        
    elif lattice_type == "woodpile":
        lattice = Lattice.tetragonal(a, a * np.sqrt(2))
        species = [material] * 4
        coords = [[0, 0, 0], [0.5, 0, 0.25], [0, 0.5, 0.5], [0.5, 0.5, 0.75]]
        
    else:
        lattice = Lattice.cubic(a)
        species = [material]
        coords = [[0.5, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Estimate bandgap wavelength
    n_eff = 1.5  # Effective refractive index
    bandgap_wavelength = 2 * a * n_eff / np.sqrt(3)
    
    return {
        "success": True,
        "lattice_type": lattice_type,
        "period_nm": period,
        "fill_fraction": fill_fraction,
        "material": material,
        "bandgap_wavelength_nm": round(bandgap_wavelength, 1),
        "n_atoms": len(structure),
        "structure": {
            "lattice": {"matrix": structure.lattice.matrix.tolist()},
            "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure]
        }
    }
