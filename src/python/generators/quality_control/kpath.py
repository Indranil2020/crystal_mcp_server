"""
quality_control/kpath.py - K-path Generation

Generates k-paths for band structure calculations:
- High-symmetry k-points
- Automatic k-path generation
- K-mesh for various dimensions
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


# High-symmetry points for different Bravais lattices
HIGH_SYMMETRY_POINTS = {
    "cubic": {
        "G": [0, 0, 0], "X": [0.5, 0, 0], "M": [0.5, 0.5, 0],
        "R": [0.5, 0.5, 0.5], "L": [0.5, 0.5, 0]
    },
    "fcc": {
        "G": [0, 0, 0], "X": [0.5, 0, 0.5], "L": [0.5, 0.5, 0.5],
        "W": [0.5, 0.25, 0.75], "U": [0.625, 0.25, 0.625], "K": [0.375, 0.375, 0.75]
    },
    "bcc": {
        "G": [0, 0, 0], "H": [0.5, -0.5, 0.5], "N": [0, 0, 0.5],
        "P": [0.25, 0.25, 0.25]
    },
    "hexagonal": {
        "G": [0, 0, 0], "M": [0.5, 0, 0], "K": [1/3, 1/3, 0],
        "A": [0, 0, 0.5], "L": [0.5, 0, 0.5], "H": [1/3, 1/3, 0.5]
    },
    "tetragonal": {
        "G": [0, 0, 0], "X": [0.5, 0, 0], "M": [0.5, 0.5, 0],
        "Z": [0, 0, 0.5], "R": [0.5, 0, 0.5], "A": [0.5, 0.5, 0.5]
    },
    "orthorhombic": {
        "G": [0, 0, 0], "X": [0.5, 0, 0], "Y": [0, 0.5, 0],
        "Z": [0, 0, 0.5], "S": [0.5, 0.5, 0], "T": [0, 0.5, 0.5],
        "U": [0.5, 0, 0.5], "R": [0.5, 0.5, 0.5]
    },
    "2D_hexagonal": {
        "G": [0, 0, 0], "M": [0.5, 0, 0], "K": [1/3, 1/3, 0]
    },
    "2D_square": {
        "G": [0, 0, 0], "X": [0.5, 0, 0], "M": [0.5, 0.5, 0]
    },
}


# Standard k-paths
STANDARD_PATHS = {
    "cubic": "G-X-M-G-R-X|M-R",
    "fcc": "G-X-W-L-G-K|U-X",
    "bcc": "G-H-N-G-P-H|P-N",
    "hexagonal": "G-M-K-G-A-L-H-A|L-M|K-H",
    "2D_hexagonal": "G-M-K-G",
    "2D_square": "G-X-M-G",
}


def get_high_symmetry_points(
    lattice_type: str = "cubic"
) -> Dict[str, Any]:
    """
    Get high-symmetry k-points for lattice type.
    
    Args:
        lattice_type: Bravais lattice type
    
    Returns:
        High-symmetry points
    """
    if lattice_type not in HIGH_SYMMETRY_POINTS:
        return {
            "success": False,
            "error": {"code": "INVALID_LATTICE", "message": f"Unknown lattice '{lattice_type}'",
                      "available": list(HIGH_SYMMETRY_POINTS.keys())}
        }
    
    points = HIGH_SYMMETRY_POINTS[lattice_type]
    
    return {
        "success": True,
        "lattice_type": lattice_type,
        "high_symmetry_points": points,
        "standard_path": STANDARD_PATHS.get(lattice_type, "G-X-M-G")
    }


def generate_kpath(
    lattice_type: str = "cubic",
    path: Optional[str] = None,
    n_points: int = 50
) -> Dict[str, Any]:
    """
    Generate k-path for band structure.
    
    Args:
        lattice_type: Bravais lattice type
        path: Custom path string (e.g., "G-X-M-G")
        n_points: Points per segment
    
    Returns:
        K-path coordinates
    """
    if lattice_type not in HIGH_SYMMETRY_POINTS:
        return {"success": False, "error": {"code": "INVALID_LATTICE", "message": f"Unknown lattice"}}
    
    hs_points = HIGH_SYMMETRY_POINTS[lattice_type]
    
    if path is None:
        path = STANDARD_PATHS.get(lattice_type, "G-X-M-G")
    
    # Parse path
    segments = []
    current_segment = []
    
    for part in path.split("|"):
        points_in_part = part.split("-")
        for i, point in enumerate(points_in_part):
            if point not in hs_points:
                return {"success": False, "error": {"code": "INVALID_POINT", "message": f"Unknown point '{point}'"}}
            
            if i > 0:
                segments.append((points_in_part[i-1], point))
    
    # Generate path
    kpoints = []
    labels = []
    distances = [0]
    
    for start, end in segments:
        k_start = np.array(hs_points[start])
        k_end = np.array(hs_points[end])
        
        for i in range(n_points):
            t = i / (n_points - 1) if n_points > 1 else 0
            k = k_start + t * (k_end - k_start)
            kpoints.append(list(k))
            
            if i == 0:
                labels.append({"index": len(kpoints) - 1, "label": start})
            elif i == n_points - 1:
                labels.append({"index": len(kpoints) - 1, "label": end})
        
        distances.append(distances[-1] + np.linalg.norm(k_end - k_start))
    
    return {
        "success": True,
        "lattice_type": lattice_type,
        "path": path,
        "n_segments": len(segments),
        "n_kpoints": len(kpoints),
        "n_points_per_segment": n_points,
        "kpoints": kpoints,
        "labels": labels,
        "total_path_length": round(distances[-1], 4)
    }


def calculate_kmesh(
    structure: Dict[str, Any] = None,
    density: float = 30.0,
    dimensions: str = "3D",
    method: str = "gamma_centered"
) -> Dict[str, Any]:
    """
    Calculate k-mesh for structure.
    
    Args:
        structure: Structure dictionary
        density: K-point density (points per Å⁻¹)
        dimensions: '3D', '2D', '1D', '0D'
        method: 'gamma_centered', 'monkhorst_pack'
    
    Returns:
        Recommended k-mesh
    """
    # Get lattice parameters
    if structure:
        lattice_info = structure.get("structure", structure).get("lattice", {})
        a = lattice_info.get("a", 10.0)
        b = lattice_info.get("b", a)
        c = lattice_info.get("c", a)
    else:
        a = b = c = 10.0
    
    # Calculate reciprocal lengths
    ka = 2 * np.pi / a
    kb = 2 * np.pi / b
    kc = 2 * np.pi / c
    
    # Calculate mesh
    if dimensions == "3D":
        nx = max(1, int(np.ceil(density * ka)))
        ny = max(1, int(np.ceil(density * kb)))
        nz = max(1, int(np.ceil(density * kc)))
    elif dimensions == "2D":
        nx = max(1, int(np.ceil(density * ka)))
        ny = max(1, int(np.ceil(density * kb)))
        nz = 1
    elif dimensions == "1D":
        nx = max(1, int(np.ceil(density * ka)))
        ny = 1
        nz = 1
    else:  # 0D
        nx = ny = nz = 1
    
    mesh = [nx, ny, nz]
    total_kpoints = nx * ny * nz
    
    # Adjust for gamma centering
    shift = [0.0, 0.0, 0.0] if method == "gamma_centered" else [0.5/nx, 0.5/ny, 0.5/nz]
    
    return {
        "success": True,
        "dimensions": dimensions,
        "method": method,
        "density_per_inv_angstrom": density,
        "mesh": mesh,
        "shift": shift,
        "total_kpoints": total_kpoints,
        "irreducible_kpoints_estimate": total_kpoints // 8 if dimensions == "3D" else total_kpoints // 4,
        "vasp_kpoints": f"{nx} {ny} {nz}",
        "qe_k_points": f"{nx} {ny} {nz} {int(shift[0])} {int(shift[1])} {int(shift[2])}"
    }
