"""
photonic/metamaterials.py - Metamaterial Generation
"""

from typing import Dict, Any, List
import numpy as np


def generate_metamaterial(
    meta_type: str = "split_ring",
    period: float = 1000.0,  # nm
    metal: str = "Au"
) -> Dict[str, Any]:
    """
    Generate metamaterial unit cell.
    
    Args:
        meta_type: Metamaterial type
        period: Unit cell period
        metal: Metal element
    
    Returns:
        Metamaterial structure
    """
    types = {
        "split_ring": {"description": "Split-ring resonator for negative μ"},
        "fishnet": {"description": "Fishnet for negative n"},
        "wire_array": {"description": "Wire array for negative ε"},
    }
    
    if meta_type not in types:
        return {"success": False, "error": {"code": "UNKNOWN", "message": f"Unknown type"}}
    
    a = period
    
    if meta_type == "split_ring":
        # Ring positions
        n_points = 20
        positions = []
        for i in range(n_points):
            if i < 17:  # Gap
                theta = 2 * np.pi * i / n_points
                x = a/2 + a/4 * np.cos(theta)
                y = a/2 + a/4 * np.sin(theta)
                positions.append([x, y, a/2])
    else:
        positions = [[a/2, a/2, a/2]]
    
    return {
        "success": True,
        "meta_type": meta_type,
        "period_nm": period,
        "metal": metal,
        "description": types[meta_type]["description"],
        "structure": {
            "atoms": [{"element": metal, "cartesian": list(p)} for p in positions]
        }
    }
