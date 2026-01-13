"""
arrangement_adapter.py - Unified Arrangement API

This module provides a clean interface for molecular cluster generation:
1. Converts molecule specifications to actual molecule structures
2. Delegates all arrangement logic to molecular_arrangement.py

Architecture:
    User Request → arrangement_adapter → molecular_arrangement.py
                                              ↓
                                      universal_molecule.py (molecule generation)

The molecular_arrangement module handles all patterns, formulas, constraints,
and validation. This adapter is a thin translation layer.

Mantra: Accuracy, Simplicity, Generic, and Modular

Author: Crystal MCP Server
License: MIT
"""

from typing import Dict, List, Optional, Any
import logging
import sys

logger = logging.getLogger(__name__)

# Debug helper - writes to stderr so stdout stays clean for JSON
def debug(msg: str) -> None:
    sys.stderr.write(f"[DEBUG arrangement_adapter] {msg}\n")
    sys.stderr.flush()

# Import molecule generation
from .universal_molecule import generate_molecule_universal

# Import the unified arrangement engine
from . import molecular_arrangement as engine

# Hunter Tracing Integration
import os
import importlib.util

if os.environ.get("ENABLE_TRACE") == "1":
    if importlib.util.find_spec("hunter"):
        from hunter import trace, Q, CallPrinter
        # Detailed Trace: Filter to reduce noise
        # 1. generated.molecule packages only
        # 2. Ignore private methods starting with __ (except __init__)? No, keep them for now.
        # 3. Use CallPrinter to show function entry associated with args.
        print("Enable Hunter Tracing (Optimized)...", file=sys.stderr)
        trace(
            # Filter: Only trace 'call' events (entering functions) to see flow,
            # or 'line' events if debugging specific logic.
            # User wants to analyze result -> 'call' with inputs/outputs is best.
            Q(module_startswith="generators.molecule", kind="call"),
            # Also trace 'return' to see what came out
            Q(module_startswith="generators.molecule", kind="return"),
            
            action=CallPrinter(
                repr_limit=500,
                stream=sys.stderr
            )
        )
    else:
        logger.warning("Hunter not installed. Tracing disabled.")


# =============================================================================
# MOLECULE GENERATION
# =============================================================================

    # ... (skipping unchanged code) ...



def _generate_molecules(
    molecule_specs: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """
    Generate actual molecule structures from specifications.

    Takes: [{"identifier": "benzene", "count": 2}, ...]
    Returns: [{"atoms": [...], "coords": [...], ...}, ...]

    This is the bridge between user-friendly specs and the arrangement engine.
    """
    molecules = []

    for spec in molecule_specs:
        identifier = spec.get("identifier", spec.get("name", ""))
        count = spec.get("count", 1)
        input_type = spec.get("input_type", "auto")

        debug(f"Resolving '{identifier}' (type={input_type}, count={count})")

        # Generate the molecule using universal_molecule
        result = generate_molecule_universal(
            identifier=identifier,
            input_type=input_type,
            optimize=spec.get("optimize", True),
            allow_external=spec.get("allow_external", True)
        )

        debug(f"  Result: success={result.get('success')}, source={result.get('source', 'N/A')}")

        if not result.get("success", False):
            error_info = result.get('error', {})
            error_msg = error_info.get('message', str(error_info)) if isinstance(error_info, dict) else str(error_info)
            debug(f"  FAILED: {error_msg}")
            logger.warning(f"Failed to generate molecule '{identifier}': {error_msg}")
            continue

        # Create molecule dict for arrangement engine
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
# MAIN API
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
    use_new_engine: bool = False,  # Kept for API compatibility (ignored)
    validate: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Generate a molecular cluster with specified arrangement.

    This function provides a unified API for molecular cluster generation.

    BASIC USAGE:
        generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 2}],
            stacking="pi_pi_parallel",
            intermolecular_distance=3.4
        )

    FORMULA-BASED:
        generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 10}],
            formulas={"x": "5*cos(2*pi*i/n)", "y": "5*sin(2*pi*i/n)", "z": "0"}
        )

    CONSTRAINT-BASED:
        generate_molecular_cluster(
            molecules=[...],
            constraints=["distance(0:centroid(), 1:centroid(), 3.4)"],
            use_solver=True
        )

    Args:
        molecules: List of molecule specifications [{identifier, count}, ...]
        stacking: Arrangement pattern (auto, pi_pi_parallel, circular, etc.)
        intermolecular_distance: Distance between molecules in Angstroms
        formulas: Position formulas {"x": ..., "y": ..., "z": ...}
        constraints: Chemical constraints ["distance(...)", ...]
        use_solver: Use constraint solver
        validate: Validate arrangement for clashes
        vacuum: Vacuum padding around the cluster

    Returns:
        Dict with structure, atoms, coords, validation info
    """
    debug(f"generate_molecular_cluster called")
    debug(f"  molecules: {molecules}")
    debug(f"  stacking: {stacking}")
    debug(f"  intermolecular_distance: {intermolecular_distance}")
    debug(f"  formulas: {formulas}")
    debug(f"  constraints: {constraints}")
    debug(f"  use_solver: {use_solver}")

    # Step 1: Generate molecules from specifications
    debug(f"Generating molecules from {len(molecules)} specifications...")
    mol_list = _generate_molecules(molecules)

    if not mol_list:
        debug("ERROR: No molecules could be generated")
        return {
            "success": False,
            "error": {
                "code": "NO_MOLECULES",
                "message": "No molecules could be generated from specifications"
            }
        }

    debug(f"Generated {len(mol_list)} molecules")

    # Step 2: Arrange molecules using the unified engine
    debug(f"Arranging molecules with pattern: {stacking}")

    # Handle parameter aliasing (distance vs intermolecular_distance)
    dist = intermolecular_distance
    if dist is None and 'distance' in kwargs:
        dist = kwargs.pop('distance')

    result = engine.arrange_molecules(
        molecules=mol_list,
        pattern=stacking,
        distance=dist,
        constraints=constraints,
        optimize=use_solver or optimize,
        validate=validate,
        formulas=formulas,
        vacuum=vacuum,
        axis=axis,  # Explicitly pass axis
        **kwargs
    )

    debug(f"arrange_molecules result: success={result.get('success')}")

    if not result.get("success"):
        debug(f"ERROR: {result.get('error')}")
        return result

    # Step 3: Build final response (compatible format)
    response = _build_response(result, intermolecular_distance or 3.4, vacuum, stacking)

    debug(f"Response built: success={response.get('success')}, n_atoms={response.get('n_atoms')}")
    return response


def _build_response(
    engine_result: Dict[str, Any],
    distance: float,
    vacuum: float,
    stacking: str
) -> Dict[str, Any]:
    """Build response in compatible format."""
    combined = engine_result.get("combined", {})
    atoms = combined.get("atoms", [])
    coords = combined.get("coords", [])
    lattice = combined.get("lattice", {"a": 20, "b": 20, "c": 20, "alpha": 90, "beta": 90, "gamma": 90})

    # Build fractional coordinates
    frac_coords = []
    cart_coords = []
    for c in coords:
        frac_coords.append([
            (c[0] + lattice["a"]/2) / lattice["a"],
            (c[1] + lattice["b"]/2) / lattice["b"],
            (c[2] + lattice["c"]/2) / lattice["c"],
        ])
        cart_coords.append([
            c[0] + lattice["a"]/2,
            c[1] + lattice["b"]/2,
            c[2] + lattice["c"]/2,
        ])

    # Build atoms array
    sites = []
    for i, (elem, frac, cart) in enumerate(zip(atoms, frac_coords, cart_coords)):
        sites.append({
            "element": elem,
            "coords": frac,
            "cartesian": cart,
            "species": [{"element": elem, "occupation": 1.0}]
        })

    metadata = engine_result.get("metadata", {})
    pattern = metadata.get("pattern", stacking)

    structure = {
        "atoms": sites,
        "sites": sites,
        "lattice": {
            "a": lattice["a"],
            "b": lattice["b"],
            "c": lattice["c"],
            "alpha": lattice.get("alpha", 90),
            "beta": lattice.get("beta", 90),
            "gamma": lattice.get("gamma", 90),
            "matrix": [
                [lattice["a"], 0, 0],
                [0, lattice["b"], 0],
                [0, 0, lattice["c"]]
            ],
            "volume": lattice["a"] * lattice["b"] * lattice["c"]
        },
        "space_group": {
            "number": 1,
            "symbol": "P1",
            "crystal_system": "triclinic"
        },
        "metadata": {
            "formula": "+".join([m.get("formula", "") for m in engine_result.get("molecules", [])]),
            "natoms": len(atoms),
            "n_molecules": engine_result.get("n_molecules", 1),
            "stacking_type": pattern,
            "intermolecular_distance": distance,
            "engine": "molecular_arrangement"
        }
    }

    return {
        "success": True,
        "n_molecules": engine_result.get("n_molecules", 1),
        "n_atoms": len(atoms),
        "atoms": atoms,
        "coords": coords,
        "formulas": [m.get("formula", "") for m in engine_result.get("molecules", [])],
        "stacking_type": pattern,
        "intermolecular_distance": distance,
        "cell": lattice,
        "structure": structure,
        "source": "molecular_arrangement",
        "validation": engine_result.get("validation", {}),
        "poses": engine_result.get("poses", []),
    }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def create_dimer(
    molecule: str,
    stacking: str = "pi_pi_parallel",
    distance: float = 3.4
) -> Dict[str, Any]:
    """Create a homo-dimer of the specified molecule."""
    return generate_molecular_cluster(
        molecules=[{"identifier": molecule, "count": 2}],
        stacking=stacking,
        intermolecular_distance=distance
    )


def create_hetero_dimer(
    molecule1: str,
    molecule2: str,
    stacking: str = "auto",
    distance: Optional[float] = None
) -> Dict[str, Any]:
    """Create a hetero-dimer of two different molecules."""
    return generate_molecular_cluster(
        molecules=[
            {"identifier": molecule1, "count": 1},
            {"identifier": molecule2, "count": 1}
        ],
        stacking=stacking,
        intermolecular_distance=distance
    )


def create_stack(
    molecules: List[str],
    distance: float = 3.4,
    stacking: str = "pi_pi_parallel"
) -> Dict[str, Any]:
    """Create a stack of multiple molecules."""
    mol_list = [{"identifier": m, "count": 1} for m in molecules]
    return generate_molecular_cluster(
        molecules=mol_list,
        stacking=stacking,
        intermolecular_distance=distance
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
        vacuum=vacuum
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
        "h_bond(donor=0, acceptor=1)"
        "plane_parallel(mol1=0, mol2=1)"

    AtomSelector targets:
        centroid()      - Molecular center
        atom(O)         - First oxygen atom
        ring_center(0)  - Center of first aromatic ring
        plane_normal()  - Normal to molecular plane
        donor_h(0)      - H-bond donor hydrogen
        acceptor(0)     - H-bond acceptor atom
    """
    return generate_molecular_cluster(
        molecules=molecules,
        stacking=initial_arrangement,
        constraints=constraints,
        use_solver=True,
        vacuum=vacuum
    )


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
]
