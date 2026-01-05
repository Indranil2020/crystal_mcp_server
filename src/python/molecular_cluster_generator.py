#!/usr/bin/env python3
"""
molecular_cluster_generator.py - MCP Tool Backend for Molecular Cluster Generation

This module is the entry point for the build_molecular_cluster MCP tool.
It receives JSON input from the TypeScript layer and returns the cluster structure.
"""

import os
import sys

# Robustness: Ensure src/python is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import json
import numpy as np
import importlib.util
from typing import Dict, Any, Optional, List

# Check availability using importlib.util.find_spec
CLUSTER_AVAILABLE = importlib.util.find_spec("generators.molecule.molecular_cluster") is not None

# Import only if available
generate_molecular_cluster = None
if CLUSTER_AVAILABLE:
    from generators.molecule.molecular_cluster import generate_molecular_cluster


def generate_cluster(
    molecules: List[Dict[str, Any]],
    stacking: str = "auto",
    intermolecular_distance: Optional[float] = None,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    rotation_x: float = 0.0,
    rotation_y: float = 0.0,
    rotation_z: float = 0.0,
    rotation_per_molecule: float = 0.0,
    axis: str = "z",
    positions: Optional[List[Dict[str, float]]] = None,
    rotations: Optional[List[Dict[str, float]]] = None,
    optimize: bool = False,
    vacuum: float = 10.0,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate a molecular cluster from the MCP tool input.
    
    This function is called by the TypeScript handler via python-bridge.
    """
    if not CLUSTER_AVAILABLE or generate_molecular_cluster is None:
        return {
            "success": False,
            "error": {
                "code": "DEPENDENCY_ERROR",
                "message": "molecular_cluster module not available",
                "details": "Install required dependencies: rdkit, scipy"
            }
        }
    
    # Call the main generator
    result = generate_molecular_cluster(
        molecules=molecules,
        stacking=stacking,
        intermolecular_distance=intermolecular_distance,
        offset_x=offset_x,
        offset_y=offset_y,
        rotation_x=rotation_x,
        rotation_y=rotation_y,
        rotation_z=rotation_z,
        rotation_per_molecule=rotation_per_molecule,
        axis=axis,
        positions=positions,
        rotations=rotations,
        optimize=optimize,
        vacuum=vacuum,
    )
    
    if not result.get("success"):
        return result
    
    # Ensure structure is in correct MCP format
    if "structure" in result:
        structure = result["structure"]
        # Add required fields for MCP response
        if "lattice" in structure and "atoms" not in structure:
            # Create atoms list from species/coords
            atoms = []
            species = result.get("atoms", [])
            coords = result.get("coords", [])
            for i, (elem, coord) in enumerate(zip(species, coords)):
                atoms.append({
                    "element": elem,
                    "cartesian": coord,
                    "species": [{"element": elem, "occupation": 1.0}]
                })
            structure["atoms"] = atoms
            structure["sites"] = atoms
    
    return result


def main():
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False, 
            "error": {"message": "Usage: python molecular_cluster_generator.py <input.json>"}
        }))
        sys.exit(1)
        
    input_file = sys.argv[1]
    with open(input_file, 'r') as f:
        params = json.load(f)
        
    # Log input for debugging
    sys.stderr.write(f"DEBUG: Input params: {json.dumps(params)}\n")
    
    result = generate_cluster(**params)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
