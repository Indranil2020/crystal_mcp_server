"""
catalyst/clusters.py - Supported Metal Clusters
"""

from typing import Dict, Any, List, Optional
import numpy as np


def generate_supported_cluster(
    cluster_element: str = "Pt",
    cluster_size: int = 13,
    support: str = "graphene",
    cluster_shape: str = "icosahedral"
) -> Dict[str, Any]:
    """
    Generate supported metal cluster.
    
    Args:
        cluster_element: Cluster metal
        cluster_size: Number of atoms in cluster
        support: Support material
        cluster_shape: Cluster geometry
    
    Returns:
        Supported cluster structure
    """
    # Lattice constants
    a_cluster = {"Pt": 3.92, "Au": 4.08, "Pd": 3.89, "Ag": 4.09}.get(cluster_element, 4.0)
    
    # Generate cluster positions
    positions = []
    
    if cluster_shape == "icosahedral" and cluster_size == 13:
        # Magic number 13-atom icosahedron
        phi = (1 + np.sqrt(5)) / 2
        r = a_cluster * 0.7
        
        # Central atom
        positions.append([0, 0, 0])
        
        # 12 vertices
        for sign1 in [-1, 1]:
            for sign2 in [-1, 1]:
                positions.append([0, sign1 * r, sign2 * r * phi])
                positions.append([sign1 * r, sign2 * r * phi, 0])
                positions.append([sign2 * r * phi, 0, sign1 * r])
    else:
        # Spherical cluster
        r = a_cluster * (cluster_size / 13) ** (1/3)
        positions.append([0, 0, 0])
        for i in range(1, cluster_size):
            theta = np.arccos(1 - 2 * i / cluster_size)
            phi = np.pi * (1 + 5**0.5) * i
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            positions.append([x, y, z])
    
    # Shift cluster above support
    cluster_height = 2.5
    positions = np.array(positions)
    positions[:, 2] += cluster_height - positions[:, 2].min()
    
    return {
        "success": True,
        "cluster_element": cluster_element,
        "cluster_size": cluster_size,
        "support": support,
        "cluster_shape": cluster_shape,
        "structure": {
            "atoms": [{"element": cluster_element, "cartesian": list(p)} for p in positions]
        }
    }
