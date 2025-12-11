"""
quantum/dots.py - Semiconductor Quantum Dots

3D quantum dots and nanocrystals.
"""

from typing import Dict, Any, List, Optional
import numpy as np


QD_DATABASE = {
    "CdSe": {"a": 6.052, "structure": "zincblende", "gap_bulk_eV": 1.74, "bohr_nm": 5.6},
    "CdS": {"a": 5.832, "structure": "zincblende", "gap_bulk_eV": 2.42, "bohr_nm": 2.8},
    "CdTe": {"a": 6.482, "structure": "zincblende", "gap_bulk_eV": 1.50, "bohr_nm": 7.3},
    "PbS": {"a": 5.936, "structure": "rocksalt", "gap_bulk_eV": 0.41, "bohr_nm": 18.0},
    "PbSe": {"a": 6.128, "structure": "rocksalt", "gap_bulk_eV": 0.28, "bohr_nm": 46.0},
    "InP": {"a": 5.869, "structure": "zincblende", "gap_bulk_eV": 1.35, "bohr_nm": 9.6},
    "InAs": {"a": 6.058, "structure": "zincblende", "gap_bulk_eV": 0.36, "bohr_nm": 34.0},
    "GaAs": {"a": 5.653, "structure": "zincblende", "gap_bulk_eV": 1.42, "bohr_nm": 11.3},
    "ZnO": {"a": 4.58, "structure": "wurtzite", "gap_bulk_eV": 3.37, "bohr_nm": 1.8},
}


def generate_quantum_dot_3d(
    material: str = "CdSe",
    diameter_nm: float = 4.0,
    shape: str = "spherical",
    core_shell: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate 3D semiconductor quantum dot.
    
    Args:
        material: QD material
        diameter_nm: Dot diameter in nm
        shape: Dot shape ('spherical', 'cubic', 'tetrahedral')
        core_shell: Shell material for core-shell QDs
    
    Returns:
        Quantum dot structure with size-dependent gap
    """
    if material not in QD_DATABASE:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown QD: {material}",
                "available": list(QD_DATABASE.keys())}}
    
    info = QD_DATABASE[material]
    a = info["a"] / 10  # Convert to nm
    gap_bulk = info["gap_bulk_eV"]
    bohr = info["bohr_nm"]
    
    # Calculate size-dependent bandgap using effective mass approximation
    R = diameter_nm / 2  # radius in nm
    
    # Simplified Brus equation: E_gap = E_bulk + (h^2 / 8R^2) * (1/me* + 1/mh*) - 1.8e^2/(4πεR)
    # Using empirical fit: ΔE ≈ 1.8 / R^2 for most QDs
    confinement = 1.8 / (R ** 2)  # eV, simplified
    gap_qd = gap_bulk + confinement
    
    # Emission wavelength
    wavelength_nm = 1240 / gap_qd
    
    # Generate spherical nanocrystal
    n_atoms = int((diameter_nm / a) ** 3 * 8)  # Approximate
    radius_A = diameter_nm * 10 / 2
    
    positions = []
    elements = []
    
    # Simple FCC construction
    n_per_side = int(diameter_nm / a) + 1
    
    M, X = material[:2], material[2:]
    if len(material) == 3:
        M, X = material[:2], material[2]
    
    for i in range(-n_per_side, n_per_side + 1):
        for j in range(-n_per_side, n_per_side + 1):
            for k in range(-n_per_side, n_per_side + 1):
                pos = np.array([i, j, k]) * a * 10 / 2
                r = np.linalg.norm(pos)
                
                if shape == "spherical" and r < radius_A:
                    if (i + j + k) % 2 == 0:
                        elements.append(M)
                    else:
                        elements.append(X)
                    positions.append(list(pos + radius_A))
                elif shape == "cubic" and abs(i) < n_per_side//2 and abs(j) < n_per_side//2 and abs(k) < n_per_side//2:
                    elements.append(M if (i + j + k) % 2 == 0 else X)
                    positions.append(list(pos + radius_A))
    
    return {
        "success": True,
        "material": material,
        "diameter_nm": diameter_nm,
        "shape": shape,
        "bulk_gap_eV": gap_bulk,
        "qd_gap_eV": round(gap_qd, 2),
        "emission_wavelength_nm": round(wavelength_nm, 1),
        "bohr_radius_nm": bohr,
        "in_strong_confinement": diameter_nm < 2 * bohr,
        "n_atoms": len(positions),
        "structure": {
            "atoms": [{"element": e, "cartesian": p} for e, p in zip(elements, positions)]
        }
    }
