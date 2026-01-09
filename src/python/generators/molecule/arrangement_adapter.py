"""
arrangement_adapter.py - Unified Arrangement API (Facade Pattern)

This module provides a unified interface for molecular arrangement that:
1. Maintains backward compatibility with existing molecular_cluster.py API
2. Routes to the new generic engines when advanced features are requested
3. Combines molecule generation (universal_molecule) with arrangement (new engines)

Architecture:
    User Request → arrangement_adapter → [molecular_cluster.py OR new engines]
                                              ↓
                                      universal_molecule.py (molecule generation)

The new engines (molecular_arrangement_engine.py and unified_molecular_arrangement.py)
ONLY handle positioning - they expect molecules to already be generated.
This adapter bridges that gap.

Author: Crystal MCP Server
License: MIT
"""

from typing import Dict, List, Optional, Any, Union
import numpy as np
import logging
import importlib.util

import sys

# New engines are in the same package - assume available, lazy import handles errors
_ENGINE_AVAILABLE = True
_UNIFIED_AVAILABLE = True

logger = logging.getLogger(__name__)

# Debug print helper - writes to stderr so stdout stays clean for JSON
def debug(msg: str) -> None:
    sys.stderr.write(f"[DEBUG arrangement_adapter] {msg}\n")
    sys.stderr.flush()

# Import legacy molecular cluster
from .molecular_cluster import (
    generate_molecular_cluster as _legacy_generate_cluster,
    StackingType,
    STACKING_DEFAULTS,
    combine_molecules,
    add_vacuum_box,
    classify_molecule,
    auto_select_stacking,
)

# Import molecule generation
from .universal_molecule import generate_molecule_universal

# Lazy import of new engines (only when needed)
_arrangement_engine = None
_unified_arrangement = None


def _get_arrangement_engine():
    """Lazy load the molecular_arrangement_engine module."""
    global _arrangement_engine
    if _arrangement_engine is None and _ENGINE_AVAILABLE:
        debug("Lazy loading molecular_arrangement_engine...")
        from . import molecular_arrangement_engine as mae
        _arrangement_engine = mae
        debug("molecular_arrangement_engine loaded successfully")
    return _arrangement_engine


def _get_unified_arrangement():
    """Lazy load the unified_molecular_arrangement module."""
    global _unified_arrangement
    if _unified_arrangement is None and _UNIFIED_AVAILABLE:
        debug("Lazy loading unified_molecular_arrangement...")
        from . import unified_molecular_arrangement as uma
        _unified_arrangement = uma
        debug("unified_molecular_arrangement loaded successfully")
    return _unified_arrangement


# =============================================================================
# PATTERN MAPPING (Legacy to New Engine)
# =============================================================================

LEGACY_TO_ENGINE_PATTERN = {
    # StackingType.value -> new engine pattern name
    "pi_pi_parallel": "pi_pi_parallel",
    "pi_pi_antiparallel": "pi_pi_antiparallel",
    "pi_pi_offset": "pi_pi_offset",
    "t_shaped": "t_shaped",
    "herringbone": "herringbone",
    "h_bonded": "h_bonded",
    "van_der_waals": "linear",  # Default to linear
    "linear": "linear",
    "circular": "circular",
    "spherical": "spherical",
    "custom": "custom",
    "swastika": "custom",  # Special handling
}


# =============================================================================
# FEATURE DETECTION
# =============================================================================

def _needs_new_engine(
    formulas: Optional[Dict] = None,
    constraints: Optional[List] = None,
    relative_poses: Optional[List] = None,
    use_solver: bool = False,
    natural_language: Optional[str] = None,
    **kwargs
) -> bool:
    """
    Determine if we should use the new arrangement engine.
    
    Returns True if any advanced feature is requested.
    """
    if formulas and any(formulas.values()):
        return True
    if constraints:
        return True
    if relative_poses:
        return True
    if use_solver:
        return True
    if natural_language:
        return True
    return False


# =============================================================================
# MOLECULE GENERATION
# =============================================================================

def _generate_molecules(
    molecule_specs: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """
    Generate actual molecule structures from specifications.
    
    Takes: [{"identifier": "benzene", "count": 2}, ...]
    Returns: [{"atoms": [...], "coords": [...], ...}, ...]
    """
    molecules = []
    
    for spec in molecule_specs:
        identifier = spec.get("identifier", spec.get("name", ""))
        count = spec.get("count", 1)
        input_type = spec.get("input_type", "auto")
        
        # Generate the molecule using universal_molecule
        result = generate_molecule_universal(
            identifier=identifier,
            input_type=input_type,
            optimize=spec.get("optimize", True),
            allow_external=spec.get("allow_external", True)
        )
        
        if not result.get("success", False):
            logger.warning(f"Failed to generate molecule '{identifier}': {result.get('error', 'unknown error')}")
            continue
        
        # Create molecule dict compatible with arrangement engines
        mol = {
            "atoms": result.get("atoms", []),
            "coords": result.get("coords", []),
            "identifier": identifier,
            "formula": result.get("formula", ""),
            "smiles": result.get("smiles", ""),
            "source": result.get("source", ""),
        }
        
        # Add count copies
        for _ in range(count):
            molecules.append(mol.copy())
    
    return molecules


# =============================================================================
# NEW ENGINE INTEGRATION
# =============================================================================

def _arrange_with_new_engine(
    molecules: List[Dict],
    arrangement: str = "auto",
    distance: Optional[float] = None,
    formulas: Optional[Dict[str, str]] = None,
    constraints: Optional[List[str]] = None,
    use_solver: bool = False,
    natural_language: Optional[str] = None,
    validate: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Arrange molecules using the new molecular_arrangement_engine.
    
    Supports:
    - Named patterns: "pi_pi_parallel", "circular", etc.
    - Position formulas: {"x": "5*cos(i)", "y": "5*sin(i)", "z": "0"}
    - Chemical constraints: ["distance(0:centroid(), 1:centroid(), 3.4)"]
    - Constraint solver optimization
    """
    engine = _get_arrangement_engine()
    if engine is None:
        raise ImportError("molecular_arrangement_engine not available")
    
    # Handle natural language parsing
    if natural_language:
        parsed = engine.parse_arrangement_request(natural_language)
        arrangement = parsed.get("pattern", arrangement)
        if parsed.get("distance"):
            distance = parsed["distance"]
        if parsed.get("n_molecules") and parsed["n_molecules"] > len(molecules):
            # Extend molecules if NL specifies more
            while len(molecules) < parsed["n_molecules"]:
                molecules.append(molecules[-1].copy())
    
    # Set default distance
    if distance is None:
        distance = 3.4
    
    # If formulas provided, use formula generator
    if formulas:
        positions = engine.generate_formula_positions(
            len(molecules),
            x_formula=formulas.get("x", "0"),
            y_formula=formulas.get("y", "0"),
            z_formula=formulas.get("z", f"{distance} * i")
        )
        orientations = engine.generate_fixed_orientations(len(molecules))
        poses = [
            engine.MoleculePose(
                position=pos,
                orientation=ori,
                molecule_idx=i
            )
            for i, (pos, ori) in enumerate(zip(positions, orientations))
        ]
    else:
        # Use generate_molecular_cluster from engine
        result = engine.generate_molecular_cluster(
            molecules=molecules,
            arrangement=arrangement,
            distance=distance,
            constraints=constraints or [],
            optimize=use_solver,
            validate=validate,
            **kwargs
        )
        return result
    
    # Build constraints if provided
    constraint_objects = []
    if constraints:
        for spec in constraints:
            if isinstance(spec, str):
                parsed = engine._parse_constraint_string(spec)
                if parsed:
                    constraint_objects.append(parsed)
            else:
                constraint_objects.append(spec)
    
    # Apply solver if requested
    if use_solver and constraint_objects:
        solver = engine.ConstraintSolver(molecules, poses, constraint_objects)
        poses = solver.solve()
    
    # Apply poses to molecules
    arranged = engine.apply_poses_to_molecules(molecules, poses)
    
    return {
        "success": True,
        "molecules": arranged,
        "poses": [
            {
                "position": [p.position.x, p.position.y, p.position.z],
                "orientation": [p.orientation.roll, p.orientation.pitch, p.orientation.yaw],
                "molecule_idx": p.molecule_idx
            }
            for p in poses
        ],
        "metadata": {
            "pattern": arrangement,
            "n_molecules": len(molecules),
            "distance": distance,
            "optimized": use_solver
        },
        "validation": {
            "valid": True  # Simplified for now
        }
    }


def _arrange_with_unified_engine(
    molecules: List[Dict],
    arrangement: str = "auto",
    distance: float = 3.4,
    **kwargs
) -> List[Dict]:
    """
    Arrange molecules using unified_molecular_arrangement.
    
    Returns: List of molecules with transformed coordinates.
    """
    unified = _get_unified_arrangement()
    if unified is None:
        raise ImportError("unified_molecular_arrangement not available")
    
    return unified.arrange_molecules_unified(
        molecules=molecules,
        arrangement=arrangement,
        distance=distance,
        **kwargs
    )


# =============================================================================
# MAIN UNIFIED API
# =============================================================================

def generate_molecular_cluster(
    molecules: List[Dict[str, Any]],
    stacking: str = "auto",
    intermolecular_distance: Optional[float] = None,
    # Legacy parameters (backward compatible)
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
    # New engine parameters
    formulas: Optional[Dict[str, str]] = None,
    constraints: Optional[List[str]] = None,
    use_solver: bool = False,
    natural_language: Optional[str] = None,
    use_new_engine: bool = False,
    validate: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Unified API for molecular cluster generation.
    
    This function provides full backward compatibility with the existing
    molecular_cluster.py API while enabling new advanced features.
    
    BASIC USAGE (Backward Compatible):
        generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 2}],
            stacking="pi_pi_parallel",
            intermolecular_distance=3.4
        )
    
    NATURAL LANGUAGE (New):
        generate_molecular_cluster(
            molecules=[...],
            natural_language="Stack 3 benzenes with π-stacking at 3.5 Å"
        )
    
    FORMULA-BASED (New):
        generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 10}],
            formulas={"x": "5*cos(2*pi*i/n)", "y": "5*sin(2*pi*i/n)", "z": "0"}
        )
    
    CONSTRAINT-BASED (New):
        generate_molecular_cluster(
            molecules=[...],
            constraints=["distance(0:ring_center(0), 1:ring_center(0), 3.4)"],
            use_solver=True
        )
    
    Args:
        molecules: List of molecule specifications
        stacking: Stacking type or pattern name
        intermolecular_distance: Distance between molecules
        
        # Legacy parameters
        offset_x, offset_y: Lateral offsets
        rotation_x, rotation_y, rotation_z: Global rotations
        rotation_per_molecule: Incremental rotation
        axis: Stacking axis
        positions, rotations: Custom positions/rotations
        optimize: Force field optimization
        vacuum: Vacuum padding
        
        # New engine parameters
        formulas: Position formulas {"x": ..., "y": ..., "z": ...}
        constraints: Chemical constraints ["distance(...)", ...]
        use_solver: Use constraint solver
        natural_language: Natural language arrangement request
        use_new_engine: Force use of new engine
        validate: Validate arrangement
    
    Returns:
        Dict with structure, atoms, coords, metadata
    """
    debug(f"generate_molecular_cluster called")
    debug(f"  molecules: {molecules}")
    debug(f"  stacking: {stacking}")
    debug(f"  use_new_engine: {use_new_engine}")
    debug(f"  formulas: {formulas}")
    debug(f"  constraints: {constraints}")
    debug(f"  use_solver: {use_solver}")
    debug(f"  natural_language: {natural_language}")
    
    # Check if we need the new engine
    needs_new = use_new_engine or _needs_new_engine(
        formulas=formulas,
        constraints=constraints,
        use_solver=use_solver,
        natural_language=natural_language
    )
    
    debug(f"needs_new_engine: {needs_new}")
    
    if needs_new:
        # Use new engine
        engine = _get_arrangement_engine()
        if engine is None:
            debug("ERROR: New engine requested but not available")
            return {
                "success": False,
                "error": {
                    "code": "ENGINE_UNAVAILABLE",
                    "message": "Advanced arrangement engine not available",
                    "details": "molecular_arrangement_engine could not be loaded"
                }
            }
        
        debug("Using new arrangement engine")
        
        # First, generate molecules from specs
        debug(f"Generating molecules from {len(molecules)} specifications...")
        mol_list = _generate_molecules(molecules)
        
        if not mol_list:
            debug("ERROR: No molecules could be generated")
            return {
                "success": False,
                "error": {"code": "NO_MOLECULES", "message": "No molecules could be generated"}
            }
        
        debug(f"Generated {len(mol_list)} molecules")
        
        # Arrange using new engine
        debug(f"Arranging molecules with pattern: {stacking}")
        result = _arrange_with_new_engine(
            molecules=mol_list,
            arrangement=stacking,
            distance=intermolecular_distance,
            formulas=formulas,
            constraints=constraints,
            use_solver=use_solver,
            natural_language=natural_language,
            validate=validate,
            **kwargs
        )
        
        debug(f"New engine result: success={result.get('success')}")
        
        # Convert to legacy output format for compatibility
        if result.get("success"):
            arranged = result.get("molecules", [])
            debug(f"Arranged molecules count: {len(arranged)}")
            
            # Combine molecules
            if arranged:
                combined = combine_molecules(arranged)
                combined = add_vacuum_box(combined, vacuum)
                
                debug(f"Combined structure: n_atoms={len(combined.get('atoms', []))}")
                
                # Build response in legacy format
                response = _build_legacy_response(
                    combined, 
                    result, 
                    vacuum,
                    intermolecular_distance or 3.4
                )
                debug(f"Built legacy response: success={response.get('success')}")
                return response
        
        # If new engine failed, propagate error (no silent fallback)
        if not result.get("success"):
            debug(f"ERROR from new engine: {result.get('error')}")
        return result
    
    # Use legacy molecular_cluster.py
    debug("Using legacy molecular_cluster.py")
    result = _legacy_generate_cluster(
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
    debug(f"Legacy result: success={result.get('success')}")
    return result


def _build_legacy_response(
    combined: Dict, 
    engine_result: Dict, 
    vacuum: float,
    distance: float
) -> Dict[str, Any]:
    """Build response in legacy format for backward compatibility."""
    cell = combined.get("cell", {"a": 20, "b": 20, "c": 20, "alpha": 90, "beta": 90, "gamma": 90})
    species = combined["atoms"]
    coords = combined["coords"]
    
    # Fractional coordinates
    frac_coords = []
    for c in coords:
        frac_coords.append([
            (c[0] + cell["a"]/2) / cell["a"],
            (c[1] + cell["b"]/2) / cell["b"],
            (c[2] + cell["c"]/2) / cell["c"],
        ])
    
    # Cartesian coordinates
    cart_coords = []
    for c in coords:
        cart_coords.append([
            c[0] + cell["a"]/2,
            c[1] + cell["b"]/2,
            c[2] + cell["c"]/2,
        ])
    
    # Build atoms array
    atoms = []
    for i, (elem, frac, cart) in enumerate(zip(species, frac_coords, cart_coords)):
        atoms.append({
            "element": elem,
            "coords": frac,
            "cartesian": cart,
            "species": [{"element": elem, "occupation": 1.0}]
        })
    
    metadata = engine_result.get("metadata", {})
    
    structure = {
        "atoms": atoms,
        "sites": atoms,
        "lattice": {
            "a": cell["a"],
            "b": cell["b"],
            "c": cell["c"],
            "alpha": cell["alpha"],
            "beta": cell["beta"],
            "gamma": cell["gamma"],
            "matrix": [
                [cell["a"], 0, 0],
                [0, cell["b"], 0],
                [0, 0, cell["c"]]
            ],
            "volume": cell["a"] * cell["b"] * cell["c"]
        },
        "space_group": {
            "number": 1,
            "symbol": "P1",
            "crystal_system": "triclinic"
        },
        "metadata": {
            "formula": "+".join(combined.get("formulas", [])),
            "natoms": len(species),
            "n_molecules": combined.get("n_molecules", 1),
            "stacking_type": metadata.get("pattern", "unknown"),
            "intermolecular_distance": distance,
            "engine": "new_arrangement_engine"
        }
    }
    
    return {
        "success": True,
        "n_molecules": combined.get("n_molecules", 1),
        "n_atoms": len(species),
        "atoms": species,
        "coords": coords,
        "formulas": combined.get("formulas", []),
        "stacking_type": metadata.get("pattern", "unknown"),
        "intermolecular_distance": distance,
        "cell": cell,
        "structure": structure,
        "source": "arrangement_engine",
        "validation": engine_result.get("validation", {}),
        "poses": engine_result.get("poses", []),
    }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def create_dimer(
    molecule: str,
    stacking: str = "pi_pi_parallel",
    distance: float = 3.4,
    use_new_engine: bool = False
) -> Dict[str, Any]:
    """Create a homo-dimer of the specified molecule."""
    return generate_molecular_cluster(
        molecules=[{"identifier": molecule, "count": 2}],
        stacking=stacking,
        intermolecular_distance=distance,
        use_new_engine=use_new_engine
    )


def create_hetero_dimer(
    molecule1: str,
    molecule2: str,
    stacking: str = "auto",
    distance: Optional[float] = None,
    use_new_engine: bool = False
) -> Dict[str, Any]:
    """Create a hetero-dimer of two different molecules."""
    return generate_molecular_cluster(
        molecules=[
            {"identifier": molecule1, "count": 1},
            {"identifier": molecule2, "count": 1}
        ],
        stacking=stacking,
        intermolecular_distance=distance,
        use_new_engine=use_new_engine
    )


def create_stack(
    molecules: List[str],
    distance: float = 3.4,
    stacking: str = "pi_pi_parallel",
    use_new_engine: bool = False
) -> Dict[str, Any]:
    """Create a stack of multiple molecules."""
    mol_list = [{"identifier": m, "count": 1} for m in molecules]
    return generate_molecular_cluster(
        molecules=mol_list,
        stacking=stacking,
        intermolecular_distance=distance,
        use_new_engine=use_new_engine
    )


def create_arrangement_from_formula(
    molecules: List[Dict[str, Any]],
    x_formula: str = "0",
    y_formula: str = "0",
    z_formula: str = "3.4 * i",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Create a molecular arrangement using position formulas.
    
    Variables available in formulas:
        i: molecule index (0, 1, 2, ...)
        n: total number of molecules
        t: normalized index = i/(n-1), range [0, 1]
        pi, e: mathematical constants
    
    Functions available:
        sin, cos, tan, sqrt, exp, log, abs, atan2
    
    Example (circular):
        create_arrangement_from_formula(
            molecules=[{"identifier": "benzene", "count": 6}],
            x_formula="5 * cos(2 * pi * i / n)",
            y_formula="5 * sin(2 * pi * i / n)",
            z_formula="0"
        )
    """
    return generate_molecular_cluster(
        molecules=molecules,
        formulas={"x": x_formula, "y": y_formula, "z": z_formula},
        vacuum=vacuum,
        use_new_engine=True
    )


def create_arrangement_with_constraints(
    molecules: List[Dict[str, Any]],
    constraints: List[str],
    initial_arrangement: str = "linear",
    vacuum: float = 10.0
) -> Dict[str, Any]:
    """
    Create a molecular arrangement satisfying chemical constraints.
    
    Constraint syntax examples:
        "distance(0:centroid(), 1:centroid(), 3.4)"
        "distance(sel1=0:ring_center(0), sel2=1:ring_center(0), target=3.5)"
        "h_bond(donor=0, acceptor=1)"
        "plane_parallel(mol1=0, mol2=1)"
    
    AtomSelector targets:
        centroid()      - Molecular center
        atom(O)         - First oxygen atom
        ring_center(0)  - Center of first aromatic ring
        plane_normal()  - Normal to molecular plane
        donor_hs()      - H-bond donor hydrogens
        acceptor_sites() - H-bond acceptor atoms
    """
    return generate_molecular_cluster(
        molecules=molecules,
        stacking=initial_arrangement,
        constraints=constraints,
        use_solver=True,
        vacuum=vacuum,
        use_new_engine=True
    )


def parse_natural_language(text: str) -> Dict[str, Any]:
    """
    Parse natural language arrangement request.
    
    Examples:
        "Stack 3 benzenes with π-stacking at 3.5 Å"
        "Arrange 6 water molecules in a circular H-bonded ring"
        "Create a T-shaped dimer at 5 angstroms"
    
    Returns:
        Dict with parsed arrangement specification
    """
    engine = _get_arrangement_engine()
    if engine is None:
        return {
            "pattern": "linear",
            "n_molecules": 2,
            "distance": 3.4,
            "error": "Arrangement engine not available"
        }
    return engine.parse_arrangement_request(text)


# =============================================================================
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Main API
    "generate_molecular_cluster",
    # Convenience functions
    "create_dimer",
    "create_hetero_dimer",
    "create_stack",
    "create_arrangement_from_formula",
    "create_arrangement_with_constraints",
    "parse_natural_language",
    # Re-exports from legacy
    "StackingType",
    "STACKING_DEFAULTS",
]
