"""
molecular_cluster.py - Molecular Cluster Generation

Generate molecular clusters for quantum chemistry calculations:
- Homo-dimers, hetero-dimers, trimers, n-mers
- π-π stacking (parallel, antiparallel, offset)
- T-shaped (edge-to-face) arrangements
- H-bonded clusters (water, alcohols, acids)
- Van der Waals complexes
- Custom arrangements with full rotation control

Designed for materials scientists, chemists, and quantum chemists.
Fully generic and modular - handles any combination of molecules.
"""

from typing import Dict, List, Optional, Any, Tuple, Union
import numpy as np
from dataclasses import dataclass
from enum import Enum
import importlib.util
import logging

# Check dependencies
SCIPY_AVAILABLE = importlib.util.find_spec("scipy") is not None
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None
PYMATGEN_AVAILABLE = importlib.util.find_spec("pymatgen") is not None

if SCIPY_AVAILABLE:
    from scipy.spatial.transform import Rotation

if PYMATGEN_AVAILABLE:
    from pymatgen.core import Molecule, Structure, Lattice

# Import existing universal molecule generator
from .universal_molecule import generate_molecule_universal

logger = logging.getLogger(__name__)


# =============================================================================
# Standard Stacking Parameters (Quantum Chemistry Standards)
# =============================================================================

class StackingType(Enum):
    """Standard stacking types in quantum chemistry."""
    PI_PI_PARALLEL = "pi_pi_parallel"         # Parallel π-stacking (face-to-face)
    PI_PI_ANTIPARALLEL = "pi_pi_antiparallel" # Antiparallel (180° rotated)
    PI_PI_OFFSET = "pi_pi_offset"             # Offset/slip-stacked
    T_SHAPED = "t_shaped"                      # Edge-to-face (perpendicular)
    HERRINGBONE = "herringbone"                # Alternating tilted arrangement
    H_BONDED = "h_bonded"                      # Hydrogen bonded
    VAN_DER_WAALS = "van_der_waals"           # General vdW contact
    LINEAR = "linear"                          # In a line
    CIRCULAR = "circular"                      # Ring arrangement
    SPHERICAL = "spherical"                    # 3D sphere distribution
    HELICAL = "helical"                        # Helical/spiral arrangement
    CUSTOM = "custom"                          # User-defined positions/rotations
    SWASTIKA = "swastika"                      # 4-molecule cross pattern


@dataclass
class StackingParameters:
    """Standard parameters for each stacking type."""
    distance: float          # Å
    offset_x: float = 0.0    # Å
    offset_y: float = 0.0    # Å
    rotation_z: float = 0.0  # degrees
    tilt: float = 0.0        # degrees (for herringbone)
    description: str = ""


# Quantum chemistry standard values
STACKING_DEFAULTS = {
    StackingType.PI_PI_PARALLEL: StackingParameters(
        distance=3.4, 
        description="Face-to-face π-stacking, typical for aromatics"
    ),
    StackingType.PI_PI_ANTIPARALLEL: StackingParameters(
        distance=3.4, 
        rotation_z=180.0,
        description="Antiparallel π-stacking with 180° rotation"
    ),
    StackingType.PI_PI_OFFSET: StackingParameters(
        distance=3.4, 
        offset_x=1.5,  # Typical slip distance
        description="Offset (slip-stacked) π-stacking"
    ),
    StackingType.T_SHAPED: StackingParameters(
        distance=5.0,  # Center-to-center for T-shaped
        tilt=90.0,
        description="Edge-to-face T-shaped arrangement"
    ),
    StackingType.HERRINGBONE: StackingParameters(
        distance=5.5,
        tilt=60.0,
        description="Herringbone pattern (common in organic crystals)"
    ),
    StackingType.H_BONDED: StackingParameters(
        distance=2.8,  # O-O distance in water
        description="Hydrogen bonded arrangement"
    ),
    StackingType.VAN_DER_WAALS: StackingParameters(
        distance=3.5,
        description="Van der Waals contact distance"
    ),
    StackingType.LINEAR: StackingParameters(
        distance=5.0,
        description="Linear arrangement along z-axis"
    ),
    StackingType.CIRCULAR: StackingParameters(
        distance=4.0,
        description="Circular ring arrangement in xy-plane"
    ),
    StackingType.SPHERICAL: StackingParameters(
        distance=5.0,
        description="Spherical distribution"
    ),
    StackingType.HELICAL: StackingParameters(
        distance=3.4,
        rotation_z=30.0,  # Default 30° per step
        description="Helical/spiral arrangement with rotation per step"
    ),
    StackingType.SWASTIKA: StackingParameters(
        distance=3.4,  # Center-to-center
        description="Swastika pattern (4-fold rotation)"
    ),
}


# =============================================================================
# Molecule Type Detection (for auto-selecting stacking)
# =============================================================================

# SMARTS patterns for molecule classification
MOLECULE_PATTERNS = {
    "aromatic": ["c", "a"],  # Any aromatic
    "h_bond_donor": ["[OH]", "[NH]", "[NH2]"],
    "h_bond_acceptor": ["[O]", "[N]", "[F]"],
    "planar": ["c1ccccc1", "c1ccncc1"],  # Benzene, pyridine
    "linear": ["C#C", "C#N", "N#N"],
}


def classify_molecule(atoms: List[str], smiles: str = "") -> Dict[str, bool]:
    """
    Classify molecule properties for auto-stacking.
    
    Returns dict with:
    - is_aromatic: Contains aromatic rings
    - is_planar: Likely planar structure
    - has_h_donors: Has H-bond donors
    - has_h_acceptors: Has H-bond acceptors
    """
    # Simple heuristics based on atoms
    has_c = "C" in atoms
    has_n = "N" in atoms
    has_o = "O" in atoms
    has_h = "H" in atoms
    has_aromatic_c = smiles.count("c") > 0 if smiles else False
    
    # Check for common aromatic patterns
    is_aromatic = has_aromatic_c or "c1" in smiles.lower() if smiles else False
    
    # Planar if aromatic and not too many sp3 carbons
    is_planar = is_aromatic and smiles.count("C") < smiles.count("c") if smiles else False
    
    # H-bond capability
    has_h_donors = has_h and (has_o or has_n)
    has_h_acceptors = has_o or has_n
    
    return {
        "is_aromatic": is_aromatic,
        "is_planar": is_planar,
        "has_h_donors": has_h_donors,
        "has_h_acceptors": has_h_acceptors,
    }


def auto_select_stacking(mol_properties: List[Dict[str, bool]]) -> StackingType:
    """
    Automatically select optimal stacking based on molecule properties.
    
    Rules (quantum chemistry best practices):
    1. All aromatic/planar → π-π parallel stacking (3.4Å)
    2. Mixed aromatic + H-bond capable → T-shaped or H-bonded
    3. All H-bond capable (water, alcohols) → H-bonded (2.8Å)
    4. Otherwise → Van der Waals contact
    """
    n_mols = len(mol_properties)
    if n_mols == 0:
        return StackingType.VAN_DER_WAALS
    
    all_aromatic = all(p["is_aromatic"] for p in mol_properties)
    any_aromatic = any(p["is_aromatic"] for p in mol_properties)
    all_h_capable = all(p["has_h_donors"] and p["has_h_acceptors"] for p in mol_properties)
    
    if all_aromatic:
        return StackingType.PI_PI_PARALLEL
    elif all_h_capable:
        return StackingType.H_BONDED
    elif any_aromatic:
        return StackingType.T_SHAPED
    else:
        return StackingType.VAN_DER_WAALS


# =============================================================================
# Coordinate Transformation Utilities
# =============================================================================

def center_molecule(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Center molecule at origin, return centered coords and original center."""
    center = np.mean(coords, axis=0)
    return coords - center, center


def rotate_molecule(
    coords: np.ndarray, 
    rotation_x: float = 0.0,
    rotation_y: float = 0.0, 
    rotation_z: float = 0.0
) -> np.ndarray:
    """
    Rotate molecule around x, y, z axes.
    
    Args:
        coords: Nx3 array of coordinates
        rotation_x: Rotation around x-axis in degrees
        rotation_y: Rotation around y-axis in degrees
        rotation_z: Rotation around z-axis in degrees
    
    Returns:
        Rotated coordinates
    """
    if not SCIPY_AVAILABLE:
        # Manual rotation matrices
        coords = _rotate_axis(coords, rotation_x, axis=0)
        coords = _rotate_axis(coords, rotation_y, axis=1)
        coords = _rotate_axis(coords, rotation_z, axis=2)
        return coords
    
    # Use scipy for rotation
    angles = [np.radians(rotation_x), np.radians(rotation_y), np.radians(rotation_z)]
    r = Rotation.from_euler('xyz', angles)
    return r.apply(coords)


def _rotate_axis(coords: np.ndarray, angle_deg: float, axis: int) -> np.ndarray:
    """Rotate coordinates around a single axis."""
    if abs(angle_deg) < 1e-6:
        return coords
    
    angle = np.radians(angle_deg)
    c, s = np.cos(angle), np.sin(angle)
    
    if axis == 0:  # X-axis
        R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    elif axis == 1:  # Y-axis
        R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    else:  # Z-axis
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    
    return coords @ R.T


def translate_molecule(coords: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Translate molecule by a vector."""
    return coords + translation


def align_to_z_axis(coords: np.ndarray) -> np.ndarray:
    """
    Align molecule's principal axis to z-axis.
    Useful for consistent stacking orientation.
    """
    # Compute principal axes using SVD
    centered, _ = center_molecule(coords)
    u, s, vh = np.linalg.svd(centered)
    
    # Align smallest principal component (normal to plane) to z
    return centered @ vh.T


# =============================================================================
# Arrangement Functions
# =============================================================================

def arrange_stacked(
    molecules: List[Dict[str, Any]],
    stacking_type: StackingType,
    distance: float,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    rotation_per_molecule: float = 0.0
) -> List[Dict[str, Any]]:
    """
    Arrange molecules in a stacked configuration along z-axis.
    
    Args:
        molecules: List of molecule dicts with 'atoms' and 'coords'
        stacking_type: Type of stacking
        distance: Inter-plane distance
        offset_x: Lateral offset in x
        offset_y: Lateral offset in y
        rotation_per_molecule: Rotation increment per molecule (degrees)
    """
    params = STACKING_DEFAULTS.get(stacking_type, STACKING_DEFAULTS[StackingType.PI_PI_PARALLEL])
    
    # Use provided or default values
    d = distance if distance else params.distance
    ox = offset_x if offset_x else params.offset_x
    oy = offset_y if offset_y else params.offset_y
    rot_z = rotation_per_molecule if rotation_per_molecule else params.rotation_z
    
    # For tilted arrangements (herringbone, T-shaped), enforce minimum distance
    # based on molecule extent to prevent atomic clashes
    if stacking_type == StackingType.HERRINGBONE and molecules:
        first_coords = np.array(molecules[0]["coords"])
        mol_extent = np.max(np.abs(first_coords)) * 2  # Approximate molecule span
        tilt_rad = np.radians(params.tilt)
        # When tilted, molecule projects further in z direction
        min_distance = mol_extent * np.abs(np.sin(tilt_rad)) + 1.0  # +1Å buffer
        if d < min_distance:
            logger.info(f"Herringbone: adjusting distance from {d:.2f} to {min_distance:.2f} Å to prevent clashes")
            d = max(d, min_distance)
    
    arranged = []
    z_position = 0.0
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        
        # Center molecule
        coords, _ = center_molecule(coords)
        
        # Align to z-axis for consistent stacking
        coords = align_to_z_axis(coords)
        
        # Apply rotation (cumulative or per-molecule)
        if stacking_type == StackingType.PI_PI_ANTIPARALLEL and i % 2 == 1:
            coords = rotate_molecule(coords, rotation_z=180.0)
        elif rot_z:
            coords = rotate_molecule(coords, rotation_z=rot_z * i)
        
        # Apply offset (alternating for some stacking types)
        if stacking_type == StackingType.PI_PI_OFFSET and i % 2 == 1:
            coords = translate_molecule(coords, np.array([ox, oy, 0]))
        elif stacking_type == StackingType.HERRINGBONE:
            tilt = params.tilt if i % 2 == 0 else -params.tilt
            coords = rotate_molecule(coords, rotation_y=tilt)
            # Add x-offset for alternating herringbone pattern
            x_offset = (mol_extent * 0.3) if i % 2 == 1 else 0.0
            coords = translate_molecule(coords, np.array([x_offset, 0, 0]))
        
        # Stack along z
        coords = translate_molecule(coords, np.array([0, 0, z_position]))
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
        
        z_position += d
    
    return arranged


def arrange_t_shaped(
    molecules: List[Dict[str, Any]],
    distance: float = 5.0
) -> List[Dict[str, Any]]:
    """
    Arrange molecules in T-shaped (edge-to-face) configuration.
    First molecule horizontal, others perpendicular.
    """
    arranged = []
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        coords = align_to_z_axis(coords)
        
        if i == 0:
            # First molecule: horizontal (in xy-plane)
            pass
        else:
            # Perpendicular molecules
            coords = rotate_molecule(coords, rotation_x=90.0)
            # Position around the first molecule
            angle = 2 * np.pi * (i - 1) / max(1, len(molecules) - 1)
            x_pos = distance * np.cos(angle)
            y_pos = distance * np.sin(angle)
            coords = translate_molecule(coords, np.array([x_pos, y_pos, 0]))
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
    
    return arranged


def arrange_circular(
    molecules: List[Dict[str, Any]],
    radius: float = 4.0,
    facing: str = "inward"  # "inward", "outward", "tangent"
) -> List[Dict[str, Any]]:
    """
    Arrange molecules in a circular pattern (good for H-bonded clusters).
    """
    n = len(molecules)
    arranged = []
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        
        # Angular position
        angle = 2 * np.pi * i / n
        
        # Position on circle
        x_pos = radius * np.cos(angle)
        y_pos = radius * np.sin(angle)
        
        # Orient molecule
        if facing == "inward":
            rot_angle = np.degrees(angle) + 180
        elif facing == "outward":
            rot_angle = np.degrees(angle)
        else:  # tangent
            rot_angle = np.degrees(angle) + 90
        
        coords = rotate_molecule(coords, rotation_z=rot_angle)
        coords = translate_molecule(coords, np.array([x_pos, y_pos, 0]))
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
    
    return arranged


def arrange_helical(
    molecules: List[Dict[str, Any]],
    distance: float = 3.4,
    radius: float = 0.0,
    rotation_per_step: float = 30.0,
    pitch: float = 3.4
) -> List[Dict[str, Any]]:
    """
    Arrange molecules in a helical (spiral) pattern.
    
    Args:
        molecules: List of molecules to arrange
        distance: Not used directly if pitch is set
        radius: Helix radius in xy-plane (0 = linear helix along z)
        rotation_per_step: Rotation in degrees per molecule
        pitch: Z-distance per molecule
    
    Returns:
        List of arranged molecules with transformed coordinates
    """
    arranged = []
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        
        # Rotation around z-axis
        rot_angle = rotation_per_step * i
        coords = rotate_molecule(coords, rotation_z=rot_angle)
        
        # Calculate position on helix
        angle_rad = np.radians(rot_angle)
        x_pos = radius * np.cos(angle_rad) if radius > 0 else 0.0
        y_pos = radius * np.sin(angle_rad) if radius > 0 else 0.0
        z_pos = pitch * i
        
        coords = translate_molecule(coords, np.array([x_pos, y_pos, z_pos]))
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
    
    return arranged


def arrange_linear(
    molecules: List[Dict[str, Any]],
    distance: float = 5.0,
    axis: str = "z",  # "x", "y", "z"
    rotations: Optional[List[Dict[str, float]]] = None  # Per-molecule rotations
) -> List[Dict[str, Any]]:
    """Arrange molecules in a line along specified axis with optional per-molecule rotations."""
    axis_idx = {"x": 0, "y": 1, "z": 2}.get(axis.lower(), 2)
    arranged = []
    
    position = 0.0
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        
        # Apply per-molecule rotation if specified
        if rotations and i < len(rotations):
            rot = rotations[i]
            coords = rotate_molecule(
                coords,
                rotation_x=rot.get("x", 0),
                rotation_y=rot.get("y", 0),
                rotation_z=rot.get("z", 0)
            )
        
        translation = np.zeros(3)
        translation[axis_idx] = position
        coords = translate_molecule(coords, translation)
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
        
        position += distance
    
    return arranged


def arrange_custom(
    molecules: List[Dict[str, Any]],
    positions: List[Dict[str, float]],
    rotations: List[Dict[str, float]]
) -> List[Dict[str, Any]]:
    """
    Arrange molecules with custom positions and rotations.
    
    Args:
        molecules: List of molecules
        positions: List of {"x": float, "y": float, "z": float} for each molecule
        rotations: List of {"x": float, "y": float, "z": float} (degrees) for each molecule
    """
    arranged = []
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        
        # Apply rotation if specified
        if i < len(rotations):
            rot = rotations[i]
            coords = rotate_molecule(
                coords,
                rotation_x=rot.get("x", 0),
                rotation_y=rot.get("y", 0),
                rotation_z=rot.get("z", 0)
            )
        
        # Apply translation if specified
        if i < len(positions):
            pos = positions[i]
            translation = np.array([
                pos.get("x", 0),
                pos.get("y", 0),
                pos.get("z", 0)
            ])
            coords = translate_molecule(coords, translation)
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
    
    return arranged


def arrange_swastika(
    molecules: List[Dict[str, Any]],
    arm_length: float = 10.0
) -> List[Dict[str, Any]]:
    """
    Arrange 4 planar molecules in a swastika pattern.
    Each molecule perpendicular to its neighbors.
    """
    if len(molecules) != 4:
        # Pad or truncate to 4
        while len(molecules) < 4:
            molecules.append(molecules[-1].copy())
        molecules = molecules[:4]
    
    arranged = []
    angles = [0, 90, 180, 270]  # Rotation for each arm
    offsets = [
        [arm_length/2, 0, 0],   # Right
        [0, arm_length/2, 0],   # Top
        [-arm_length/2, 0, 0],  # Left
        [0, -arm_length/2, 0],  # Bottom
    ]
    
    for i, mol in enumerate(molecules):
        coords = np.array(mol["coords"])
        coords, _ = center_molecule(coords)
        coords = align_to_z_axis(coords)
        
        # Rotate to point outward
        coords = rotate_molecule(coords, rotation_z=angles[i])
        
        # Translate to arm position
        coords = translate_molecule(coords, np.array(offsets[i]))
        
        arranged.append({
            "atoms": mol["atoms"],
            "coords": coords.tolist(),
            "identifier": mol.get("identifier", f"molecule_{i}"),
            "formula": mol.get("formula", ""),
        })
    
    return arranged


# =============================================================================
# Combine Molecules into Single Structure
# =============================================================================

def combine_molecules(molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Combine arranged molecules into a single structure."""
    all_atoms = []
    all_coords = []
    
    for mol in molecules:
        all_atoms.extend(mol["atoms"])
        all_coords.extend(mol["coords"])
    
    return {
        "atoms": all_atoms,
        "coords": all_coords,
        "n_molecules": len(molecules),
        "formulas": [mol.get("formula", "") for mol in molecules],
    }


def add_vacuum_box(combined: Dict[str, Any], vacuum: float = 10.0) -> Dict[str, Any]:
    """Add vacuum box around the cluster."""
    coords = np.array(combined["coords"])
    
    # Find bounding box
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    
    # Center the structure
    center = (min_coords + max_coords) / 2
    coords = coords - center
    
    # Calculate box size
    extent = max_coords - min_coords
    box_size = extent + 2 * vacuum
    
    combined["coords"] = coords.tolist()
    combined["cell"] = {
        "a": float(box_size[0]),
        "b": float(box_size[1]),
        "c": float(box_size[2]),
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    
    return combined


# =============================================================================
# Main Cluster Generation Function
# =============================================================================

def generate_molecular_cluster(
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
    allow_external: bool = True
) -> Dict[str, Any]:
    """
    Generate a molecular cluster from multiple molecules.
    
    This function:
    1. Generates each molecule using universal_molecule.py
    2. Classifies molecules for auto-stacking selection
    3. Arranges molecules according to specified/auto-detected stacking
    4. Combines into a single structure with vacuum box
    
    Args:
        molecules: List of molecule specifications:
            [{"identifier": "benzene", "count": 2}, {"identifier": "water"}]
        
        stacking: Stacking type or "auto" for automatic selection:
            - "auto": Auto-detect based on molecule types
            - "pi_pi_parallel": Parallel π-stacking
            - "pi_pi_antiparallel": Antiparallel π-stacking
            - "pi_pi_offset": Offset/slip-stacked
            - "t_shaped": Edge-to-face
            - "herringbone": Alternating tilted
            - "h_bonded": Hydrogen bonded
            - "van_der_waals": General vdW contact
            - "linear": In a line
            - "circular": Ring arrangement
            - "swastika": 4-molecule cross pattern
            - "custom": Use provided positions/rotations
        
        intermolecular_distance: Distance between molecules (Å).
            If None, uses standard value for stacking type.
        
        offset_x, offset_y: Lateral offsets (Å)
        
        rotation_x, rotation_y, rotation_z: Global rotation (degrees)
            Applied to entire cluster after arrangement.
        
        rotation_per_molecule: Incremental rotation per molecule (degrees)
        
        axis: Axis for linear arrangement ("x", "y", "z")
        
        positions: Custom positions for each molecule
            [{"x": 0, "y": 0, "z": 0}, {"x": 5, "y": 0, "z": 3.5}, ...]
        
        rotations: Custom rotations for each molecule (degrees)
            [{"x": 0, "y": 0, "z": 0}, {"x": 90, "y": 0, "z": 45}, ...]
        
        optimize: Optimize cluster geometry with force field
        
        vacuum: Vacuum padding (Å)
        
        allow_external: Allow online molecule lookup
    
    Returns:
        Dictionary with cluster structure:
        {
            "success": True,
            "n_molecules": 3,
            "n_atoms": 42,
            "atoms": [...],
            "coords": [...],
            "stacking_type": "pi_pi_parallel",
            "intermolecular_distance": 3.4,
            "cell": {...},
            "structure": {...}  # For MCP response
        }
    """
    # Validate input
    if not molecules:
        return {
            "success": False,
            "error": {
                "code": "NO_MOLECULES",
                "message": "No molecules specified",
                "suggestion": "Provide at least one molecule in the 'molecules' array"
            }
        }
    
    # Step 1: Generate each molecule
    mol_structures = []
    mol_properties = []
    
    for mol_spec in molecules:
        identifier = mol_spec.get("identifier", mol_spec.get("name", ""))
        count = mol_spec.get("count", 1)
        input_type = mol_spec.get("input_type", "auto")
        
        if not identifier:
            continue
        
        # Use existing universal molecule generation
        result = generate_molecule_universal(
            identifier,
            input_type=input_type,
            optimize=True,
            allow_external=allow_external
        )
        
        if not result.get("success"):
            return {
                "success": False,
                "error": {
                    "code": "MOLECULE_GENERATION_FAILED",
                    "message": f"Failed to generate molecule: {identifier}",
                    "details": result.get("error"),
                    "suggestion": "Check molecule name or provide SMILES directly"
                }
            }
        
        # Classify molecule
        props = classify_molecule(
            result["atoms"],
            result.get("smiles", "")
        )
        
        # Add copies
        for _ in range(count):
            mol_structures.append({
                "atoms": result["atoms"],
                "coords": result["coords"],
                "identifier": identifier,
                "formula": result.get("formula", identifier),
                "smiles": result.get("smiles", ""),
            })
            mol_properties.append(props)
    
    if not mol_structures:
        return {
            "success": False,
            "error": {
                "code": "NO_VALID_MOLECULES",
                "message": "No valid molecules could be generated",
            }
        }
    
    # Step 2: Determine stacking type
    if not stacking:
        stacking = "auto"
        
    stacking_lower = stacking.lower().replace("-", "_").replace(" ", "_")
    
    if stacking_lower == "auto":
        # Heuristic: If axis is X or Y, user likely implies LINEAR arrangement
        # (Standard stacking is usually Z/face-to-face)
        if axis.lower() in ["x", "y"]:
            stacking_type = StackingType.LINEAR
            logger.info(f"Auto-selected stacking: {stacking_type.value} (inferred from axis={axis})")
        else:
            stacking_type = auto_select_stacking(mol_properties)
            logger.info(f"Auto-selected stacking: {stacking_type.value}")
    else:
        # Map string to enum
        stacking_map = {
            "pi_pi_parallel": StackingType.PI_PI_PARALLEL,
            "parallel": StackingType.PI_PI_PARALLEL,
            "stacked": StackingType.PI_PI_PARALLEL,
            "pi_pi_antiparallel": StackingType.PI_PI_ANTIPARALLEL,
            "antiparallel": StackingType.PI_PI_ANTIPARALLEL,
            "pi_pi_offset": StackingType.PI_PI_OFFSET,
            "offset": StackingType.PI_PI_OFFSET,
            "slip_stacked": StackingType.PI_PI_OFFSET,
            "t_shaped": StackingType.T_SHAPED,
            "t-shaped": StackingType.T_SHAPED,
            "edge_to_face": StackingType.T_SHAPED,
            "herringbone": StackingType.HERRINGBONE,
            "h_bonded": StackingType.H_BONDED,
            "hydrogen_bonded": StackingType.H_BONDED,
            "van_der_waals": StackingType.VAN_DER_WAALS,
            "vdw": StackingType.VAN_DER_WAALS,
            "linear": StackingType.LINEAR,
            "circular": StackingType.CIRCULAR,
            "ring": StackingType.CIRCULAR,
            "spherical": StackingType.SPHERICAL,
            "custom": StackingType.CUSTOM,
            "swastika": StackingType.SWASTIKA,
            "swastic": StackingType.SWASTIKA,  # Typo tolerance
        }
        stacking_type = stacking_map.get(stacking_lower, StackingType.PI_PI_PARALLEL)
    
    # Track warnings for user feedback
    warnings = []
    
    # Get default distance if not specified
    if intermolecular_distance is None:
        params = STACKING_DEFAULTS.get(stacking_type, STACKING_DEFAULTS[StackingType.PI_PI_PARALLEL])
        intermolecular_distance = params.distance
        
        # For LINEAR stacking, calculate smart distance based on molecule size
        if stacking_type == StackingType.LINEAR and mol_structures:
            # Get the largest molecule extent
            max_extent = 0
            for mol in mol_structures:
                coords = np.array(mol["coords"])
                extent = np.max(coords, axis=0) - np.min(coords, axis=0)
                max_extent = max(max_extent, np.max(extent))
            
            # Smart distance = molecule size + 2Å buffer
            smart_distance = max_extent + 2.0
            if smart_distance > intermolecular_distance:
                logger.info(f"Linear stacking: auto-calculated distance {smart_distance:.2f}Å based on molecule size")
                warnings.append(f"Distance auto-calculated to {smart_distance:.1f}Å based on molecule size (~{max_extent:.1f}Å)")
                intermolecular_distance = smart_distance
    
    # Step 3: Arrange molecules
    if stacking_type == StackingType.SWASTIKA:
        arranged = arrange_swastika(mol_structures, arm_length=intermolecular_distance * 2)
    elif stacking_type == StackingType.CUSTOM and positions:
        arranged = arrange_custom(mol_structures, positions, rotations or [])
    elif stacking_type == StackingType.T_SHAPED:
        arranged = arrange_t_shaped(mol_structures, distance=intermolecular_distance)
    elif stacking_type == StackingType.CIRCULAR:
        arranged = arrange_circular(mol_structures, radius=intermolecular_distance)
    elif stacking_type == StackingType.HELICAL:
        arranged = arrange_helical(
            mol_structures, 
            distance=intermolecular_distance,
            rotation_per_step=rotation_per_molecule if rotation_per_molecule else 30.0,
            pitch=intermolecular_distance
        )
    elif stacking_type == StackingType.LINEAR:
        arranged = arrange_linear(mol_structures, distance=intermolecular_distance, axis=axis, rotations=rotations)
    else:
        # Stacked arrangements
        arranged = arrange_stacked(
            mol_structures,
            stacking_type=stacking_type,
            distance=intermolecular_distance,
            offset_x=offset_x,
            offset_y=offset_y,
            rotation_per_molecule=rotation_per_molecule
        )
    
    # Step 4: Combine into single structure
    combined = combine_molecules(arranged)
    
    # Step 5: Apply global rotation if specified
    if rotation_x or rotation_y or rotation_z:
        coords = np.array(combined["coords"])
        coords = rotate_molecule(coords, rotation_x, rotation_y, rotation_z)
        combined["coords"] = coords.tolist()
    
    # Step 6: Add vacuum box
    combined = add_vacuum_box(combined, vacuum)
    
    # Step 7: Optimize if requested (using RDKit or ASE)
    if optimize and RDKIT_AVAILABLE:
        combined = _optimize_cluster_rdkit(combined)
    
    # Build response structure
    cell = combined.get("cell", {"a": 20, "b": 20, "c": 20, "alpha": 90, "beta": 90, "gamma": 90})
    
    # Create structure dict for MCP response
    species = combined["atoms"]
    coords = combined["coords"]
    
    # Fractional coordinates (centered in box)
    frac_coords = []
    for c in coords:
        frac_coords.append([
            (c[0] + cell["a"]/2) / cell["a"],
            (c[1] + cell["b"]/2) / cell["b"],
            (c[2] + cell["c"]/2) / cell["c"],
        ])
    
    # Cartesian coordinates (also centered in box)
    cart_coords = []
    for c in coords:
        cart_coords.append([
            c[0] + cell["a"]/2,
            c[1] + cell["b"]/2,
            c[2] + cell["c"]/2,
        ])
    
    # Build atoms array for GUI viewer (requires element, coords, cartesian)
    atoms = []
    for i, (elem, frac, cart) in enumerate(zip(species, frac_coords, cart_coords)):
        atoms.append({
            "element": elem,
            "coords": frac,
            "cartesian": cart,
            "species": [{"element": elem, "occupation": 1.0}]
        })
    
    structure = {
        "atoms": atoms,
        "sites": atoms,  # Alias for compatibility
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
            "formula": "+".join(combined["formulas"]),
            "natoms": len(species),
            "n_molecules": combined["n_molecules"],
            "stacking_type": stacking_type.value,
            "intermolecular_distance": intermolecular_distance,
            "warnings": warnings,
        }
    }
    
    return {
        "success": True,
        "n_molecules": combined["n_molecules"],
        "n_atoms": len(species),
        "atoms": species,
        "coords": coords,
        "formulas": combined["formulas"],
        "stacking_type": stacking_type.value,
        "intermolecular_distance": intermolecular_distance,
        "cell": cell,
        "structure": structure,
        "source": "molecular_cluster",
    }


def _optimize_cluster_rdkit(combined: Dict[str, Any]) -> Dict[str, Any]:
    """Optimize cluster geometry using RDKit (if available)."""
    if not RDKIT_AVAILABLE:
        return combined
    
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # Convert to RDKit molecule (simplified - just use UFF on combined)
    # Note: This is a simplified optimization; for production, consider using
    # a proper force field package like ASE or OpenMM
    
    return combined  # For now, return unchanged


# =============================================================================
# Convenience Functions
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
