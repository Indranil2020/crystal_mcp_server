#!/usr/bin/env python3
"""
molecular_cluster_generator.py - MCP Tool Backend for Molecular Cluster Generation

This module is the entry point for the build_molecular_cluster MCP tool.
It receives JSON input from the TypeScript layer and returns the cluster structure.

UPDATED: Now routes through arrangement_adapter to enable advanced features:
- Constraint-based arrangements (distance, angle, H-bond)
- Formula-based position generation
- Constraint solver optimization
"""

import os
import os
import sys
import tempfile
import time

# Robustness: Ensure src/python is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import json
import numpy as np
import importlib.util
from typing import Dict, Any, Optional, List

# Debug print helper - writes to stderr so stdout stays clean for JSON
def debug(msg: str) -> None:
    sys.stderr.write(f"[DEBUG molecular_cluster_generator] {msg}\n")
    sys.stderr.flush()

# Check availability - prefer advanced adapter, fallback to basic
ADAPTER_AVAILABLE = importlib.util.find_spec("generators.molecule.arrangement_adapter") is not None
CLUSTER_AVAILABLE = importlib.util.find_spec("generators.molecule.molecular_cluster") is not None

debug(f"Module availability check: adapter={ADAPTER_AVAILABLE}, cluster={CLUSTER_AVAILABLE}")

# Import the appropriate module
generate_molecular_cluster = None
if ADAPTER_AVAILABLE:
    debug("Using advanced arrangement_adapter")
    from generators.molecule.arrangement_adapter import generate_molecular_cluster
elif CLUSTER_AVAILABLE:
    debug("Fallback: Using basic molecular_cluster (advanced features unavailable)")
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
    output_file: Optional[str] = None,
    # NEW: Advanced engine parameters
    formulas: Optional[Dict[str, str]] = None,
    constraints: Optional[List[str]] = None,
    use_solver: bool = False,
    natural_language: Optional[str] = None,
    use_new_engine: bool = False,
    validate: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate a molecular cluster from the MCP tool input.
    
    This function is called by the TypeScript handler via python-bridge.
    
    Parameters:
        molecules: List of molecule specifications [{identifier, count}, ...]
        stacking: Stacking pattern (auto, pi_pi_parallel, t_shaped, etc.)
        intermolecular_distance: Distance between molecules in Angstroms
        
        # Legacy parameters
        offset_x, offset_y: Lateral offsets
        rotation_x, rotation_y, rotation_z: Global rotations
        rotation_per_molecule: Incremental rotation
        axis: Stacking axis
        positions, rotations: Custom positions/rotations
        optimize: Force field optimization
        vacuum: Vacuum padding
        
        # NEW: Advanced parameters
        formulas: Position formulas {"x": "...", "y": "...", "z": "..."}
        constraints: Chemical constraints ["distance(...)", ...]
        use_solver: Use constraint solver
        natural_language: Natural language arrangement request
        use_new_engine: Force use of new engine
        validate: Validate arrangement
    """
    debug(f"generate_cluster called")
    debug(f"  molecules: {molecules}")
    debug(f"  stacking: {stacking}")
    debug(f"  intermolecular_distance: {intermolecular_distance}")
    debug(f"  formulas: {formulas}")
    debug(f"  constraints: {constraints}")
    debug(f"  use_solver: {use_solver}")
    debug(f"  natural_language: {natural_language}")
    debug(f"  use_new_engine: {use_new_engine}")
    
    if generate_molecular_cluster is None:
        error_msg = "Neither arrangement_adapter nor molecular_cluster module available"
        debug(f"ERROR: {error_msg}")
        return {
            "success": False,
            "error": {
                "code": "DEPENDENCY_ERROR",
                "message": error_msg,
                "details": "Install required dependencies: rdkit, scipy"
            }
        }
    
    # Determine if we should force new engine based on advanced parameters
    has_advanced_params = bool(formulas or constraints or use_solver or natural_language)
    if has_advanced_params:
        debug(f"Advanced parameters detected, forcing use_new_engine=True")
        use_new_engine = True
    
    debug(f"Calling generate_molecular_cluster (use_new_engine={use_new_engine})")
    
    # Handle rotation arrays: if rotation_x/y/z is a list, convert to rotations format
    # This allows LLM to specify per-molecule rotations via rotation_x: [0, 0, 90] etc.
    per_mol_rotations = rotations  # Start with existing rotations if any
    
    # Check if any rotation param is a list (per-molecule)
    rx_is_list = isinstance(rotation_x, list)
    ry_is_list = isinstance(rotation_y, list)
    rz_is_list = isinstance(rotation_z, list)
    
    if rx_is_list or ry_is_list or rz_is_list:
        # Determine number of molecules from the lists
        n_mol = max(
            len(rotation_x) if rx_is_list else 0,
            len(rotation_y) if ry_is_list else 0,
            len(rotation_z) if rz_is_list else 0
        )
        if n_mol == 0:
            # Fallback to counting molecules from input
            n_mol = sum(m.get('count', 1) for m in molecules)
        
        # Build per-molecule rotations list
        per_mol_rotations = []
        for i in range(n_mol):
            rx = rotation_x[i] if rx_is_list and i < len(rotation_x) else (0.0 if rx_is_list else rotation_x)
            ry = rotation_y[i] if ry_is_list and i < len(rotation_y) else (0.0 if ry_is_list else rotation_y)
            rz = rotation_z[i] if rz_is_list and i < len(rotation_z) else (0.0 if rz_is_list else rotation_z)
            per_mol_rotations.append({'x': rx, 'y': ry, 'z': rz})
        
        debug(f"  Converted rotation arrays to per_mol_rotations: {per_mol_rotations}")
        # Reset the global rotation values since we're using per-molecule
        rotation_x = 0.0
        rotation_y = 0.0
        rotation_z = 0.0
    
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
        rotations=per_mol_rotations,  # Use converted per-molecule rotations
        optimize=optimize,
        vacuum=vacuum,
        # Advanced parameters
        formulas=formulas,
        constraints=constraints,
        use_solver=use_solver,
        natural_language=natural_language,
        use_new_engine=use_new_engine,
        validate=validate,
    )
    
    debug(f"generate_molecular_cluster returned: success={result.get('success')}")
    
    if not result.get("success"):
        debug(f"ERROR from generator: {result.get('error')}")
        return result
    
    debug(f"Result summary: n_molecules={result.get('n_molecules')}, n_atoms={result.get('n_atoms')}")
    debug(f"Stacking type used: {result.get('stacking_type')}")
    
    # Ensure structure is in correct MCP format
    if "structure" in result:
        structure = result["structure"]
        debug(f"Structure has keys: {list(structure.keys())}")
        # Add required fields for MCP response
        if "lattice" in structure and "atoms" not in structure:
            debug("Adding atoms array to structure for MCP format")
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
    
    # Add validation info if available
    if "validation" in result:
        validation = result["validation"]
        debug(f"Validation result: valid={validation.get('valid')}, warnings={validation.get('warnings')}, errors={validation.get('errors')}")
    
    debug("generate_cluster completed successfully")
    
    # Universal Auto-Save Logic
    # Always save to a file to ensure robust chaining, even if not requested
    saved_file_path = output_file
    
    if not saved_file_path and result.get("success"):
        # Auto-generate temp file
        timestamp = int(time.time())
        temp_dir = tempfile.gettempdir()
        temp_filename = f"molecular_cluster_{timestamp}_{os.getpid()}.json"
        saved_file_path = os.path.join(temp_dir, temp_filename)
        debug(f"Auto-offloading result to temp file: {saved_file_path}")

    # Save to file (either requested or auto-generated)
    if saved_file_path and result.get("success"):
        with open(saved_file_path, 'w') as f:
            json.dump(result, f, indent=2)
        debug(f"Saved result to {saved_file_path}")
        result["output_file"] = saved_file_path
        result["auto_saved_path"] = saved_file_path # Explicitly mark as auto-saved

    
    # Context Management / Frontend Compatibility
    # If the structure is extremely large (> 500 atoms), remove it from the direct JSON output
    # to prevent crashing the LLM context. Otherwise, keep it for the Frontend.
    n_atoms = result.get("n_atoms", 0)
    if n_atoms > 500:
        debug(f"Structure is large ({n_atoms} atoms). Truncating from response payload to save context.")
        if "structure" in result:
            del result["structure"]
        result["truncated"] = True
        result["message"] = f"Structure ({n_atoms} atoms) too large for chat. Saved to {saved_file_path}"
            
    return result


def main():
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False, 
            "error": {"message": "Usage: python molecular_cluster_generator.py <input.json>"}
        }))
        sys.exit(1)
        
    input_file = sys.argv[1]
    debug(f"Reading input from: {input_file}")
    
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    debug(f"Input params: {json.dumps(params, indent=2)}")
    
    result = generate_cluster(**params)
    
    debug(f"Final result success: {result.get('success')}")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
