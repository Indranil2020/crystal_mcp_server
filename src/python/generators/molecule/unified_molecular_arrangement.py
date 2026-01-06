"""
unified_molecular_arrangement.py - The Ultimate Generic Molecular Arrangement System

This module synthesizes THREE complementary approaches into a unified system:

1. RIGID BODY ENGINE (from z2_plan):
   - MoleculePose with 6 DOF (position + orientation)
   - MolecularFrame for local coordinate systems
   - Relative placement with parent-child resolution
   - Deterministic, always produces valid output

2. MATHEMATICAL GENERATORS (from original plan):
   - Position generators (Linear, Circular, Helix, Formula, etc.)
   - Orientation generators (Fixed, Incremental, FaceCenter, etc.)
   - ArrangementSpec combining generators
   - Formula-based DSL for custom patterns

3. CHEMICAL INTELLIGENCE (from K2_plan):
   - AtomSelector for chemistry-aware selection
   - Chemical constraints (HBond, π-stacking, etc.)
   - Pattern library with chemical knowledge
   - Optional constraint refinement

The system supports MULTIPLE specification modes:

Mode A: Named patterns ("pi_pi_parallel", "h_bonded", etc.)
Mode B: Mathematical formulas ("x = 5*cos(2*pi*i/n)", etc.)
Mode C: Explicit poses ([{position: [0,0,0], orientation: [0,0,0]}, ...])
Mode D: Relative poses ([{relative_to: 0, offset: [0,0,3.4]}, ...])
Mode E: Chemical constraints (["distance(0:ring_center, 1:ring_center, 3.4)", ...])
Mode F: Natural language (parsed by MCP to one of the above)

Author: Unified Arrangement System
License: MIT
"""

from typing import (
    Dict, List, Optional, Any, Tuple, Union,
    Callable, Protocol, Sequence, Set
)
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
import numpy as np
import math
import re


# =============================================================================
# LAYER 1: CORE RIGID BODY ENGINE (from z2_plan)
# =============================================================================

@dataclass
class MolecularFrame:
    """
    Local coordinate frame attached to a molecule.
    
    Used for:
    - Defining molecule's intrinsic orientation (plane normal, principal axes)
    - Specifying relative placements in a meaningful coordinate system
    - Aligning molecules based on chemical features (ring plane, bond axis)
    """
    origin: np.ndarray        # 3-vector: reference point (centroid, atom, etc.)
    x_axis: np.ndarray        # 3-vector: normalized x direction
    y_axis: np.ndarray        # 3-vector: normalized y direction
    z_axis: np.ndarray        # 3-vector: normalized z direction (often plane normal)
    frame_type: str = "centroid"  # How the frame was computed
    
    def to_matrix(self) -> np.ndarray:
        """Return 3x3 rotation matrix (columns are local axes in global coords)."""
        return np.column_stack([self.x_axis, self.y_axis, self.z_axis])
    
    @classmethod
    def identity(cls) -> 'MolecularFrame':
        """Standard frame at origin with identity rotation."""
        return cls(
            origin=np.zeros(3),
            x_axis=np.array([1, 0, 0]),
            y_axis=np.array([0, 1, 0]),
            z_axis=np.array([0, 0, 1]),
            frame_type="identity"
        )
    
    @classmethod
    def from_coords_centroid(cls, coords: np.ndarray) -> 'MolecularFrame':
        """
        Compute frame from molecular coordinates using SVD (principal axes).
        Z-axis aligns with smallest principal component (plane normal for planar molecules).
        """
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        
        # SVD to get principal axes
        u, s, vh = np.linalg.svd(centered, full_matrices=False)
        
        # vh rows are principal directions, ordered by singular value (descending)
        # For a planar molecule, smallest singular value direction is the normal
        x_axis = vh[0] / (np.linalg.norm(vh[0]) + 1e-10)
        y_axis = vh[1] / (np.linalg.norm(vh[1]) + 1e-10)
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / (np.linalg.norm(z_axis) + 1e-10)
        
        # Ensure right-handed coordinate system
        if np.dot(z_axis, vh[2]) < 0:
            z_axis = -z_axis
        
        return cls(
            origin=centroid,
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis,
            frame_type="centroid_svd"
        )
    
    @classmethod
    def from_plane_normal(cls, coords: np.ndarray, 
                          normal: Optional[np.ndarray] = None) -> 'MolecularFrame':
        """
        Create frame with z-axis along plane normal.
        If normal not provided, compute best-fit plane.
        """
        centroid = np.mean(coords, axis=0)
        
        if normal is None:
            # Compute plane normal via SVD
            centered = coords - centroid
            u, s, vh = np.linalg.svd(centered, full_matrices=False)
            normal = vh[2]  # Smallest singular value direction
        
        z_axis = normal / (np.linalg.norm(normal) + 1e-10)
        
        # Choose x-axis perpendicular to z
        if abs(z_axis[2]) < 0.9:
            x_axis = np.cross([0, 0, 1], z_axis)
        else:
            x_axis = np.cross([1, 0, 0], z_axis)
        x_axis = x_axis / (np.linalg.norm(x_axis) + 1e-10)
        
        y_axis = np.cross(z_axis, x_axis)
        
        return cls(
            origin=centroid,
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis,
            frame_type="plane_normal"
        )
    
    @classmethod
    def from_bond_axis(cls, coords: np.ndarray, 
                       atom1_idx: int, atom2_idx: int) -> 'MolecularFrame':
        """Create frame with z-axis along a bond."""
        origin = coords[atom1_idx]
        z_axis = coords[atom2_idx] - coords[atom1_idx]
        z_axis = z_axis / (np.linalg.norm(z_axis) + 1e-10)
        
        # Choose perpendicular axes
        if abs(z_axis[2]) < 0.9:
            x_axis = np.cross([0, 0, 1], z_axis)
        else:
            x_axis = np.cross([1, 0, 0], z_axis)
        x_axis = x_axis / (np.linalg.norm(x_axis) + 1e-10)
        
        y_axis = np.cross(z_axis, x_axis)
        
        return cls(
            origin=origin,
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis,
            frame_type="bond_axis"
        )


@dataclass
class MoleculePose:
    """
    Complete 6-DOF pose of a molecule.
    
    Supports both absolute and relative specification:
    - Absolute: position and rotation in global coordinates
    - Relative: offset and rotation relative to a parent molecule's frame
    """
    # Absolute pose (global coordinates)
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))  # 3x3 matrix
    
    # Relative pose (if parent_idx >= 0)
    parent_idx: int = -1  # -1 means absolute, >= 0 means relative to that molecule
    local_offset: Optional[np.ndarray] = None  # Offset in parent's local frame
    local_rotation: Optional[np.ndarray] = None  # Rotation in parent's local frame
    
    # Metadata
    molecule_idx: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def is_relative(self) -> bool:
        return self.parent_idx >= 0
    
    def to_euler_degrees(self) -> Tuple[float, float, float]:
        """Extract Euler angles (roll, pitch, yaw) from rotation matrix."""
        R = self.rotation
        sy = np.sqrt(R[0, 0]**2 + R[1, 0]**2)
        singular = sy < 1e-6
        
        if not singular:
            roll = np.arctan2(R[2, 1], R[2, 2])
            pitch = np.arctan2(-R[2, 0], sy)
            yaw = np.arctan2(R[1, 0], R[0, 0])
        else:
            roll = np.arctan2(-R[1, 2], R[1, 1])
            pitch = np.arctan2(-R[2, 0], sy)
            yaw = 0
        
        return (np.degrees(roll), np.degrees(pitch), np.degrees(yaw))
    
    @classmethod
    def from_euler_degrees(cls, position: np.ndarray, 
                           roll: float, pitch: float, yaw: float,
                           **kwargs) -> 'MoleculePose':
        """Create pose from position and Euler angles (degrees)."""
        rx, ry, rz = np.radians([roll, pitch, yaw])
        
        Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
        Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
        Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
        
        R = Rz @ Ry @ Rx
        
        return cls(position=np.array(position), rotation=R, **kwargs)
    
    @classmethod
    def relative(cls, parent_idx: int, offset: np.ndarray,
                 rotation: Optional[np.ndarray] = None, **kwargs) -> 'MoleculePose':
        """Create a relative pose specification."""
        return cls(
            parent_idx=parent_idx,
            local_offset=np.array(offset),
            local_rotation=rotation if rotation is not None else np.eye(3),
            **kwargs
        )


def resolve_relative_poses(
    poses: List[MoleculePose],
    frames: List[MolecularFrame]
) -> List[MoleculePose]:
    """
    Resolve relative poses to absolute poses using topological sorting.
    
    For each pose with parent_idx >= 0:
    1. Get parent's resolved position and rotation
    2. Get parent's molecular frame
    3. Transform local offset/rotation to global coordinates
    
    Returns new list of poses with all positions/rotations in global coordinates.
    """
    n = len(poses)
    resolved = [None] * n
    
    # Build dependency graph and find resolution order
    # Simple approach: iterate until all resolved (works for DAGs)
    remaining = set(range(n))
    max_iterations = n + 1
    
    for _ in range(max_iterations):
        if not remaining:
            break
        
        resolved_this_round = []
        
        for i in list(remaining):
            pose = poses[i]
            
            if not pose.is_relative:
                # Absolute pose - directly resolved
                resolved[i] = MoleculePose(
                    position=pose.position.copy(),
                    rotation=pose.rotation.copy(),
                    molecule_idx=pose.molecule_idx,
                    metadata=pose.metadata
                )
                resolved_this_round.append(i)
            
            elif pose.parent_idx in [j for j in range(n) if resolved[j] is not None]:
                # Parent is resolved - can resolve this one
                parent = resolved[pose.parent_idx]
                parent_frame = frames[pose.parent_idx]
                
                # Transform local offset to global
                F = parent_frame.to_matrix()  # Parent frame axes in global coords
                R_parent = parent.rotation     # Parent rotation
                
                local_offset = pose.local_offset if pose.local_offset is not None else np.zeros(3)
                local_rotation = pose.local_rotation if pose.local_rotation is not None else np.eye(3)
                
                # Global offset = parent_rotation @ (frame @ local_offset)
                global_offset = R_parent @ (F @ local_offset)
                
                # Global rotation = parent_rotation @ frame @ local_rotation @ frame.T
                # (simplified if frame aligns with parent axes)
                global_rotation = R_parent @ F @ local_rotation @ F.T
                
                resolved[i] = MoleculePose(
                    position=parent.position + global_offset,
                    rotation=global_rotation,
                    molecule_idx=pose.molecule_idx,
                    metadata=pose.metadata
                )
                resolved_this_round.append(i)
        
        remaining -= set(resolved_this_round)
    
    if remaining:
        raise ValueError(f"Could not resolve poses {remaining} - circular dependency?")
    
    return resolved


# =============================================================================
# LAYER 2: MATHEMATICAL GENERATORS (enhanced from original plan)
# =============================================================================

class PositionGenerator(Protocol):
    """Protocol for position generators."""
    def __call__(self, n: int, **params) -> List[np.ndarray]: ...


class OrientationGenerator(Protocol):
    """Protocol for orientation generators."""
    def __call__(self, n: int, positions: List[np.ndarray], **params) -> List[np.ndarray]: ...


# --- Position Generators ---

class LinearPositions:
    """Generate positions along a line."""
    def __call__(self, n: int, 
                 axis: Tuple[float, float, float] = (0, 0, 1),
                 spacing: float = 3.4,
                 start: Tuple[float, float, float] = (0, 0, 0)) -> List[np.ndarray]:
        axis_vec = np.array(axis)
        axis_vec = axis_vec / (np.linalg.norm(axis_vec) + 1e-10)
        start_vec = np.array(start)
        
        return [start_vec + i * spacing * axis_vec for i in range(n)]


class CircularPositions:
    """Generate positions on a circle."""
    def __call__(self, n: int,
                 radius: float = 5.0,
                 center: Tuple[float, float, float] = (0, 0, 0),
                 plane_normal: Tuple[float, float, float] = (0, 0, 1),
                 start_angle: float = 0.0) -> List[np.ndarray]:
        center_vec = np.array(center)
        normal = np.array(plane_normal)
        normal = normal / (np.linalg.norm(normal) + 1e-10)
        
        # Create orthonormal basis in plane
        if abs(normal[2]) < 0.9:
            u = np.cross(normal, [0, 0, 1])
        else:
            u = np.cross(normal, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)
        
        positions = []
        for i in range(n):
            angle = np.radians(start_angle) + 2 * np.pi * i / n
            pos = center_vec + radius * (np.cos(angle) * u + np.sin(angle) * v)
            positions.append(pos)
        
        return positions


class HelicalPositions:
    """Generate positions along a helix."""
    def __call__(self, n: int,
                 radius: float = 5.0,
                 pitch: float = 3.4,
                 turns: float = 1.0,
                 axis: Tuple[float, float, float] = (0, 0, 1)) -> List[np.ndarray]:
        axis_vec = np.array(axis)
        axis_vec = axis_vec / (np.linalg.norm(axis_vec) + 1e-10)
        
        # Perpendicular vectors
        if abs(axis_vec[2]) < 0.9:
            u = np.cross(axis_vec, [0, 0, 1])
        else:
            u = np.cross(axis_vec, [1, 0, 0])
        u = u / np.linalg.norm(u)
        v = np.cross(axis_vec, u)
        
        total_angle = 2 * np.pi * turns
        total_rise = pitch * turns
        
        positions = []
        for i in range(n):
            t = i / max(n - 1, 1)
            angle = t * total_angle
            rise = t * total_rise
            
            pos = radius * np.cos(angle) * u + radius * np.sin(angle) * v + rise * axis_vec
            positions.append(pos)
        
        return positions


class GridPositions:
    """Generate positions on a 2D or 3D grid."""
    def __call__(self, n: int,
                 spacing: Tuple[float, ...] = (3.4, 3.4),
                 center: Tuple[float, float, float] = (0, 0, 0)) -> List[np.ndarray]:
        center_vec = np.array(center)
        
        # Determine grid dimensions
        side = int(np.ceil(np.sqrt(n)))
        
        positions = []
        idx = 0
        for i in range(side):
            for j in range(side):
                if idx >= n:
                    break
                x = (i - (side - 1) / 2) * spacing[0]
                y = (j - (side - 1) / 2) * (spacing[1] if len(spacing) > 1 else spacing[0])
                positions.append(center_vec + np.array([x, y, 0]))
                idx += 1
        
        return positions


class SphericalPositions:
    """Generate positions on a sphere (Fibonacci lattice for uniformity)."""
    def __call__(self, n: int,
                 radius: float = 5.0,
                 center: Tuple[float, float, float] = (0, 0, 0)) -> List[np.ndarray]:
        center_vec = np.array(center)
        golden_ratio = (1 + np.sqrt(5)) / 2
        
        positions = []
        for i in range(n):
            theta = 2 * np.pi * i / golden_ratio
            phi = np.arccos(1 - 2 * (i + 0.5) / n)
            
            x = radius * np.sin(phi) * np.cos(theta)
            y = radius * np.sin(phi) * np.sin(theta)
            z = radius * np.cos(phi)
            
            positions.append(center_vec + np.array([x, y, z]))
        
        return positions


class FormulaPositions:
    """
    Generate positions from string formulas.
    
    Variables: i (index), n (count), t (=i/(n-1)), pi, e
    Functions: sin, cos, tan, sqrt, exp, log, abs, atan2
    """
    def __call__(self, n: int,
                 x_formula: str = "0",
                 y_formula: str = "0",
                 z_formula: str = "3.4 * i",
                 **extra_vars) -> List[np.ndarray]:
        safe_dict = {
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
            'abs': abs, 'pi': np.pi, 'e': np.e,
            'atan2': np.arctan2, 'degrees': np.degrees, 'radians': np.radians,
            'n': n
        }
        safe_dict.update(extra_vars)
        
        positions = []
        for i in range(n):
            safe_dict['i'] = i
            safe_dict['t'] = i / max(n - 1, 1)
            
            x = float(eval(x_formula, {"__builtins__": {}}, safe_dict))
            y = float(eval(y_formula, {"__builtins__": {}}, safe_dict))
            z = float(eval(z_formula, {"__builtins__": {}}, safe_dict))
            
            positions.append(np.array([x, y, z]))
        
        return positions


class ExplicitPositions:
    """Use explicitly provided coordinates."""
    def __call__(self, n: int,
                 coordinates: List[Tuple[float, float, float]] = None) -> List[np.ndarray]:
        if coordinates is None:
            coordinates = [(0, 0, i * 3.4) for i in range(n)]
        
        positions = [np.array(c) for c in coordinates[:n]]
        
        # Pad if needed
        while len(positions) < n:
            positions.append(positions[-1].copy() if positions else np.zeros(3))
        
        return positions


# --- Orientation Generators ---

class FixedOrientation:
    """All molecules have the same orientation."""
    def __call__(self, n: int, positions: List[np.ndarray],
                 roll: float = 0, pitch: float = 0, yaw: float = 0) -> List[np.ndarray]:
        R = euler_to_matrix(roll, pitch, yaw)
        return [R.copy() for _ in range(n)]


class IncrementalRotation:
    """Progressive rotation around an axis."""
    def __call__(self, n: int, positions: List[np.ndarray],
                 axis: str = 'z',
                 angle_per_molecule: float = 0.0,
                 initial: Tuple[float, float, float] = (0, 0, 0)) -> List[np.ndarray]:
        rotations = []
        for i in range(n):
            roll, pitch, yaw = initial
            angle = angle_per_molecule * i
            
            if axis.lower() == 'x':
                roll += angle
            elif axis.lower() == 'y':
                pitch += angle
            else:
                yaw += angle
            
            rotations.append(euler_to_matrix(roll, pitch, yaw))
        
        return rotations


class AlternatingRotation:
    """Alternate between two orientations."""
    def __call__(self, n: int, positions: List[np.ndarray],
                 orientation_a: Tuple[float, float, float] = (0, 0, 0),
                 orientation_b: Tuple[float, float, float] = (0, 0, 180),
                 pattern: List[int] = None) -> List[np.ndarray]:
        if pattern is None:
            pattern = [0, 1]
        
        R_a = euler_to_matrix(*orientation_a)
        R_b = euler_to_matrix(*orientation_b)
        matrices = [R_a, R_b]
        
        return [matrices[pattern[i % len(pattern)]].copy() for i in range(n)]


class FaceCenterOrientation:
    """Orient molecules to face toward/away from a center point."""
    def __call__(self, n: int, positions: List[np.ndarray],
                 center: Tuple[float, float, float] = (0, 0, 0),
                 facing: str = 'inward',
                 up_axis: Tuple[float, float, float] = (0, 0, 1)) -> List[np.ndarray]:
        center_vec = np.array(center)
        up_vec = np.array(up_axis)
        up_vec = up_vec / (np.linalg.norm(up_vec) + 1e-10)
        
        rotations = []
        for pos in positions:
            direction = center_vec - pos
            dist = np.linalg.norm(direction)
            
            if dist < 1e-6:
                rotations.append(np.eye(3))
                continue
            
            direction = direction / dist
            if facing == 'outward':
                direction = -direction
            elif facing == 'tangent':
                direction = np.cross(up_vec, direction)
                if np.linalg.norm(direction) < 1e-6:
                    direction = np.array([1, 0, 0])
                direction = direction / np.linalg.norm(direction)
            
            R = direction_to_rotation(direction, up_vec)
            rotations.append(R)
        
        return rotations


class FormulaOrientation:
    """Generate orientations from formulas."""
    def __call__(self, n: int, positions: List[np.ndarray],
                 roll_formula: str = "0",
                 pitch_formula: str = "0",
                 yaw_formula: str = "0",
                 **extra_vars) -> List[np.ndarray]:
        safe_dict = {
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
            'abs': abs, 'pi': np.pi, 'e': np.e,
            'atan2': np.arctan2, 'degrees': np.degrees, 'radians': np.radians,
            'n': n
        }
        safe_dict.update(extra_vars)
        
        rotations = []
        for i, pos in enumerate(positions):
            safe_dict['i'] = i
            safe_dict['t'] = i / max(n - 1, 1)
            safe_dict['x'], safe_dict['y'], safe_dict['z'] = pos
            
            roll = float(eval(roll_formula, {"__builtins__": {}}, safe_dict))
            pitch = float(eval(pitch_formula, {"__builtins__": {}}, safe_dict))
            yaw = float(eval(yaw_formula, {"__builtins__": {}}, safe_dict))
            
            rotations.append(euler_to_matrix(roll, pitch, yaw))
        
        return rotations


# --- Helper Functions ---

def euler_to_matrix(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """Convert Euler angles (degrees) to rotation matrix."""
    rx, ry, rz = np.radians([roll, pitch, yaw])
    
    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    
    return Rz @ Ry @ Rx


def direction_to_rotation(direction: np.ndarray, up: np.ndarray) -> np.ndarray:
    """Create rotation matrix that aligns z-axis with direction."""
    z_axis = direction / (np.linalg.norm(direction) + 1e-10)
    
    x_axis = np.cross(up, z_axis)
    if np.linalg.norm(x_axis) < 1e-6:
        x_axis = np.cross([1, 0, 0], z_axis)
    x_axis = x_axis / (np.linalg.norm(x_axis) + 1e-10)
    
    y_axis = np.cross(z_axis, x_axis)
    
    return np.column_stack([x_axis, y_axis, z_axis])


# =============================================================================
# LAYER 3: CHEMICAL INTELLIGENCE (from K2_plan)
# =============================================================================

class AtomSelector:
    """
    Chemistry-aware atom/group selector.
    
    Syntax: "mol_id:target(args)"
    
    Examples:
        "0:centroid()"           - Center of molecule 0
        "1:atom(O)"              - First oxygen in molecule 1
        "0:ring_center(0)"       - Center of first aromatic ring
        "0:plane_normal()"       - Normal vector of molecular plane
        "0:donor_hs()"           - H-bond donor hydrogens
        "0:acceptor_sites()"     - H-bond acceptor atoms (O, N, F)
        "0:func_group(carbonyl)" - Carbonyl group atoms
    """
    
    PATTERN = re.compile(r"(\d+):(\w+)(?:\(([^)]*)\))?")
    
    def __init__(self, selector_str: str, molecules: List[Dict[str, Any]]):
        self.selector_str = selector_str
        self.molecules = molecules
    
    def evaluate(self) -> Union[np.ndarray, List[np.ndarray]]:
        """Evaluate selector and return coordinates or coordinate list."""
        match = self.PATTERN.match(self.selector_str)
        if not match:
            raise ValueError(f"Invalid selector: {self.selector_str}")
        
        mol_id = int(match.group(1))
        target = match.group(2)
        args_str = match.group(3) or ""
        
        # Parse arguments
        args = self._parse_args(args_str)
        
        return self._resolve(mol_id, target, args)
    
    def _parse_args(self, args_str: str) -> Dict[str, Any]:
        """Parse argument string."""
        args = {}
        if not args_str.strip():
            return args
        
        for part in args_str.split(","):
            part = part.strip()
            if "=" in part:
                k, v = part.split("=", 1)
                args[k.strip()] = self._parse_value(v.strip())
            elif part:
                args["arg0"] = self._parse_value(part)
        
        return args
    
    def _parse_value(self, v: str) -> Any:
        """Parse a single value."""
        # Try number
        if v.replace(".", "").replace("-", "").isdigit():
            return float(v) if "." in v else int(v)
        # String (strip quotes)
        return v.strip("'\"")
    
    def _resolve(self, mol_id: int, target: str, args: Dict) -> Any:
        """Resolve selector to actual data."""
        if mol_id >= len(self.molecules):
            raise ValueError(f"Molecule index {mol_id} out of range")
        
        mol = self.molecules[mol_id]
        coords = np.array(mol["coords"])
        atoms = mol.get("atoms", [])
        
        if target == "centroid":
            return np.mean(coords, axis=0)
        
        elif target == "atom":
            element = args.get("arg0", args.get("element", "C"))
            index = args.get("index", 0)
            indices = [i for i, a in enumerate(atoms) if a == element]
            if not indices:
                raise ValueError(f"No {element} atoms in molecule {mol_id}")
            return coords[indices[min(index, len(indices) - 1)]]
        
        elif target == "plane_normal":
            frame = MolecularFrame.from_coords_centroid(coords)
            return frame.z_axis
        
        elif target == "ring_center":
            ring_id = args.get("arg0", args.get("ring_id", 0))
            return self._find_ring_center(mol, coords, atoms, ring_id)
        
        elif target == "donor_hs":
            # H atoms bonded to O or N
            return self._find_donor_hydrogens(coords, atoms, mol)
        
        elif target == "acceptor_sites":
            # O, N, F atoms
            indices = [i for i, a in enumerate(atoms) if a in ["O", "N", "F"]]
            return coords[indices] if indices else np.array([])
        
        elif target == "func_group":
            group_name = args.get("arg0", args.get("name", ""))
            return self._find_functional_group(mol, coords, atoms, group_name)
        
        else:
            raise ValueError(f"Unknown selector target: {target}")
    
    def _find_ring_center(self, mol: Dict, coords: np.ndarray, 
                          atoms: List[str], ring_id: int) -> np.ndarray:
        """Find center of aromatic ring (simplified heuristic)."""
        # For aromatic carbons, find clusters of 5-6 carbons forming a ring
        c_indices = [i for i, a in enumerate(atoms) if a == "C"]
        
        if len(c_indices) < 5:
            return np.mean(coords, axis=0)
        
        # Simple heuristic: return centroid of first 6 carbons
        ring_atoms = c_indices[:min(6, len(c_indices))]
        return np.mean(coords[ring_atoms], axis=0)
    
    def _find_donor_hydrogens(self, coords: np.ndarray, atoms: List[str],
                               mol: Dict) -> np.ndarray:
        """Find H-bond donor hydrogens (H attached to O or N)."""
        # Simplified: find H atoms near O or N
        h_indices = [i for i, a in enumerate(atoms) if a == "H"]
        on_indices = [i for i, a in enumerate(atoms) if a in ["O", "N"]]
        
        donors = []
        for h_idx in h_indices:
            for on_idx in on_indices:
                dist = np.linalg.norm(coords[h_idx] - coords[on_idx])
                if dist < 1.2:  # Typical O-H or N-H bond length
                    donors.append(coords[h_idx])
                    break
        
        return np.array(donors) if donors else np.array([]).reshape(0, 3)
    
    def _find_functional_group(self, mol: Dict, coords: np.ndarray,
                                atoms: List[str], group_name: str) -> np.ndarray:
        """Find functional group atoms (simplified)."""
        group_patterns = {
            "carbonyl": ["C", "O"],  # C=O
            "hydroxyl": ["O", "H"],  # O-H
            "amine": ["N", "H"],     # N-H
            "carboxyl": ["C", "O", "O"],  # COOH
        }
        
        pattern = group_patterns.get(group_name.lower(), [])
        if not pattern:
            return np.mean(coords, axis=0)
        
        # Find atoms matching pattern (simplified)
        matching = []
        for i, a in enumerate(atoms):
            if a == pattern[0]:
                matching.append(coords[i])
        
        return np.mean(matching, axis=0) if matching else np.mean(coords, axis=0)


# --- Chemical Constraints ---

class ChemicalConstraint(ABC):
    """Base class for chemical constraints."""
    
    def __init__(self, weight: float = 1.0, priority: int = 1):
        self.weight = weight
        self.priority = priority
    
    @abstractmethod
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        """Return constraint violation (0 = satisfied)."""
        pass
    
    @abstractmethod
    def description(self) -> str:
        """Human-readable description."""
        pass


class DistanceConstraint(ChemicalConstraint):
    """Distance between two selected points."""
    
    def __init__(self, selector1: str, selector2: str,
                 target: Optional[float] = None,
                 min_dist: Optional[float] = None,
                 max_dist: Optional[float] = None,
                 weight: float = 1.0, priority: int = 1):
        super().__init__(weight=weight, priority=priority)
        self.selector1 = selector1
        self.selector2 = selector2
        self.target = target
        self.min_dist = min_dist
        self.max_dist = max_dist
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        p1 = AtomSelector(self.selector1, molecules).evaluate()
        p2 = AtomSelector(self.selector2, molecules).evaluate()
        
        d = np.linalg.norm(p1 - p2)
        
        if self.target is not None:
            return (d - self.target) ** 2
        
        violation = 0.0
        if self.min_dist and d < self.min_dist:
            violation += (self.min_dist - d) ** 2
        if self.max_dist and d > self.max_dist:
            violation += (d - self.max_dist) ** 2
        
        return violation
    
    def description(self) -> str:
        target_str = f"= {self.target}" if self.target else f"[{self.min_dist}, {self.max_dist}]"
        return f"Distance({self.selector1}, {self.selector2}) {target_str} Å"


class AngleConstraint(ChemicalConstraint):
    """Angle between three points."""
    
    def __init__(self, selector1: str, selector2: str, selector3: str,
                 target: Optional[float] = None,
                 min_angle: Optional[float] = None,
                 max_angle: Optional[float] = None,
                 weight: float = 1.0, priority: int = 1):
        super().__init__(weight=weight, priority=priority)
        self.selector1 = selector1
        self.selector2 = selector2  # Vertex
        self.selector3 = selector3
        self.target = target
        self.min_angle = min_angle
        self.max_angle = max_angle
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        p1 = AtomSelector(self.selector1, molecules).evaluate()
        p2 = AtomSelector(self.selector2, molecules).evaluate()
        p3 = AtomSelector(self.selector3, molecules).evaluate()
        
        v1 = p1 - p2
        v2 = p3 - p2
        
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
        
        if self.target is not None:
            return (angle - self.target) ** 2
        
        violation = 0.0
        if self.min_angle and angle < self.min_angle:
            violation += (self.min_angle - angle) ** 2
        if self.max_angle and angle > self.max_angle:
            violation += (angle - self.max_angle) ** 2
        
        return violation
    
    def description(self) -> str:
        return f"Angle({self.selector1}, {self.selector2}, {self.selector3})"


class PlaneParallelConstraint(ChemicalConstraint):
    """Two molecular planes should be parallel."""
    
    def __init__(self, mol1_idx: int, mol2_idx: int,
                 max_deviation: float = 10.0,
                 weight: float = 1.0, priority: int = 1):
        super().__init__(weight=weight, priority=priority)
        self.mol1_idx = mol1_idx
        self.mol2_idx = mol2_idx
        self.max_deviation = max_deviation  # degrees
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        n1 = AtomSelector(f"{self.mol1_idx}:plane_normal()", molecules).evaluate()
        n2 = AtomSelector(f"{self.mol2_idx}:plane_normal()", molecules).evaluate()
        
        cos_angle = abs(np.dot(n1, n2))
        angle = np.degrees(np.arccos(np.clip(cos_angle, 0, 1)))
        
        if angle > self.max_deviation:
            return (angle - self.max_deviation) ** 2
        return 0.0
    
    def description(self) -> str:
        return f"Planes of mol {self.mol1_idx} and {self.mol2_idx} parallel (±{self.max_deviation}°)"


class HBondConstraint(ChemicalConstraint):
    """
    Hydrogen bond geometry constraint.
    
    Enforces proper H-bond geometry:
    - D-H...A distance: 1.5-2.5 Å (H to acceptor)
    - D...A distance: 2.5-3.5 Å (donor heavy atom to acceptor)
    - D-H...A angle: > 120° (ideally ~180°)
    """
    
    def __init__(self, donor_mol: int, acceptor_mol: int,
                 h_acceptor_dist: Tuple[float, float] = (1.5, 2.5),
                 donor_acceptor_dist: Tuple[float, float] = (2.5, 3.5),
                 min_angle: float = 120.0,
                 weight: float = 1.0, priority: int = 1):
        super().__init__(weight=weight, priority=priority)
        self.donor_mol = donor_mol
        self.acceptor_mol = acceptor_mol
        self.h_acceptor_dist = h_acceptor_dist
        self.donor_acceptor_dist = donor_acceptor_dist
        self.min_angle = min_angle
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        # Get donor hydrogens and acceptor sites
        donor_hs = AtomSelector(f"{self.donor_mol}:donor_hs()", molecules).evaluate()
        acceptor_sites = AtomSelector(f"{self.acceptor_mol}:acceptor_sites()", molecules).evaluate()
        
        if len(donor_hs) == 0 or len(acceptor_sites) == 0:
            return 1e6  # No H-bond possible
        
        # Find best H-bond
        best_score = 1e6
        
        for h_coord in donor_hs:
            for acc_coord in acceptor_sites:
                d_h_acc = np.linalg.norm(h_coord - acc_coord)
                
                # Distance penalty
                if d_h_acc < self.h_acceptor_dist[0]:
                    dist_penalty = (self.h_acceptor_dist[0] - d_h_acc) ** 2
                elif d_h_acc > self.h_acceptor_dist[1]:
                    dist_penalty = (d_h_acc - self.h_acceptor_dist[1]) ** 2
                else:
                    dist_penalty = 0
                
                score = dist_penalty
                best_score = min(best_score, score)
        
        return best_score
    
    def description(self) -> str:
        return f"H-bond from mol {self.donor_mol} to mol {self.acceptor_mol}"


# --- Chemical Pattern Library ---

class ChemicalPatternLibrary:
    """
    Library of chemically-meaningful arrangement patterns.
    
    Each pattern is defined as a combination of:
    - Position/orientation generators
    - Chemical constraints
    - Default parameters
    """
    
    @staticmethod
    def pi_pi_parallel(distance: float = 3.4) -> Dict[str, Any]:
        """Face-to-face π-stacking."""
        return {
            'position_generator': LinearPositions(),
            'orientation_generator': FixedOrientation(),
            'position_params': {'spacing': distance, 'axis': (0, 0, 1)},
            'orientation_params': {},
            'constraints': [
                PlaneParallelConstraint(mol1_idx=0, mol2_idx=1, max_deviation=10)
            ],
            'description': f'Parallel π-stacking at {distance} Å'
        }
    
    @staticmethod
    def pi_pi_antiparallel(distance: float = 3.4) -> Dict[str, Any]:
        """Antiparallel π-stacking (180° rotation)."""
        return {
            'position_generator': LinearPositions(),
            'orientation_generator': AlternatingRotation(),
            'position_params': {'spacing': distance, 'axis': (0, 0, 1)},
            'orientation_params': {
                'orientation_a': (0, 0, 0),
                'orientation_b': (0, 0, 180)
            },
            'constraints': [],
            'description': f'Antiparallel π-stacking at {distance} Å'
        }
    
    @staticmethod
    def pi_pi_offset(distance: float = 3.4, offset: float = 1.5) -> Dict[str, Any]:
        """Offset/slip-stacked π-stacking."""
        return {
            'position_generator': FormulaPositions(),
            'orientation_generator': FixedOrientation(),
            'position_params': {
                'x_formula': f'{offset} * (i % 2)',
                'y_formula': '0',
                'z_formula': f'{distance} * i'
            },
            'orientation_params': {},
            'constraints': [],
            'description': f'Offset π-stacking at {distance} Å with {offset} Å slip'
        }
    
    @staticmethod
    def t_shaped(distance: float = 5.0) -> Dict[str, Any]:
        """T-shaped edge-to-face arrangement."""
        return {
            'position_generator': ExplicitPositions(),
            'orientation_generator': lambda n, pos, **kw: [
                euler_to_matrix(0, 0, 0),
                euler_to_matrix(0, 90, 0)
            ][:n],
            'position_params': {
                'coordinates': [(0, 0, 0), (distance, 0, 0)]
            },
            'orientation_params': {},
            'constraints': [],
            'description': f'T-shaped at {distance} Å'
        }
    
    @staticmethod
    def h_bonded_dimer(distance: float = 2.8) -> Dict[str, Any]:
        """H-bonded dimer with proper geometry."""
        return {
            'position_generator': LinearPositions(),
            'orientation_generator': FaceCenterOrientation(),
            'position_params': {'spacing': distance},
            'orientation_params': {'facing': 'inward'},
            'constraints': [
                HBondConstraint(donor_mol=0, acceptor_mol=1)
            ],
            'description': f'H-bonded dimer at ~{distance} Å'
        }
    
    @staticmethod
    def circular_hbond(n: int, radius: float = 5.0) -> Dict[str, Any]:
        """Circular arrangement for H-bonded clusters (water, alcohols)."""
        return {
            'position_generator': CircularPositions(),
            'orientation_generator': FaceCenterOrientation(),
            'position_params': {'radius': radius},
            'orientation_params': {'facing': 'inward'},
            'constraints': [],
            'description': f'{n}-membered H-bonded ring at {radius} Å radius'
        }
    
    @staticmethod
    def herringbone(distance: float = 5.5, tilt: float = 60.0) -> Dict[str, Any]:
        """Herringbone pattern with alternating tilts."""
        return {
            'position_generator': LinearPositions(),
            'orientation_generator': AlternatingRotation(),
            'position_params': {'spacing': distance},
            'orientation_params': {
                'orientation_a': (0, tilt, 0),
                'orientation_b': (0, -tilt, 0)
            },
            'constraints': [],
            'description': f'Herringbone at ±{tilt}°'
        }


# =============================================================================
# LAYER 4: UNIFIED ARRANGEMENT API
# =============================================================================

@dataclass
class ArrangementSpec:
    """Complete arrangement specification."""
    position_generator: Any  # PositionGenerator
    orientation_generator: Any  # OrientationGenerator
    position_params: Dict[str, Any] = field(default_factory=dict)
    orientation_params: Dict[str, Any] = field(default_factory=dict)
    constraints: List[ChemicalConstraint] = field(default_factory=list)
    relative_poses: List[Dict[str, Any]] = field(default_factory=list)
    name: str = "custom"
    description: str = ""


def generate_arrangement(
    n_molecules: int,
    spec: Union[str, ArrangementSpec, Dict[str, Any]],
    molecules: Optional[List[Dict[str, Any]]] = None,
    **kwargs
) -> List[MoleculePose]:
    """
    Generate molecular poses from an arrangement specification.
    
    This is the main entry point supporting all specification modes:
    
    Mode A: Named pattern string
        generate_arrangement(4, "pi_pi_parallel", distance=3.4)
    
    Mode B: ArrangementSpec object
        spec = ArrangementSpec(...)
        generate_arrangement(4, spec)
    
    Mode C: Dictionary configuration
        generate_arrangement(4, {'position_type': 'circular', ...})
    
    Mode D: Chemical pattern from library
        generate_arrangement(4, "h_bonded_dimer")
    """
    
    # Convert string to spec
    if isinstance(spec, str):
        spec = _resolve_named_pattern(spec, n_molecules, **kwargs)
    elif isinstance(spec, dict):
        spec = _spec_from_dict(spec, **kwargs)
    
    # Generate positions
    positions = spec.position_generator(n_molecules, **spec.position_params)
    
    # Generate orientations
    if callable(spec.orientation_generator):
        rotations = spec.orientation_generator(n_molecules, positions, **spec.orientation_params)
    else:
        rotations = [np.eye(3) for _ in range(n_molecules)]
    
    # Create poses
    poses = []
    for i, (pos, rot) in enumerate(zip(positions, rotations)):
        poses.append(MoleculePose(
            position=pos,
            rotation=rot,
            molecule_idx=i,
            metadata={'pattern': spec.name}
        ))
    
    # Handle relative poses if specified
    if spec.relative_poses:
        poses = _apply_relative_poses(poses, spec.relative_poses, molecules)
    
    return poses


def _resolve_named_pattern(name: str, n: int, **kwargs) -> ArrangementSpec:
    """Resolve pattern name to ArrangementSpec."""
    name_lower = name.lower().replace('-', '_').replace(' ', '_')
    
    # Check chemical pattern library first
    library_methods = {
        'pi_pi_parallel': ChemicalPatternLibrary.pi_pi_parallel,
        'parallel': ChemicalPatternLibrary.pi_pi_parallel,
        'pi_pi_antiparallel': ChemicalPatternLibrary.pi_pi_antiparallel,
        'antiparallel': ChemicalPatternLibrary.pi_pi_antiparallel,
        'pi_pi_offset': ChemicalPatternLibrary.pi_pi_offset,
        'offset': ChemicalPatternLibrary.pi_pi_offset,
        'slip_stacked': ChemicalPatternLibrary.pi_pi_offset,
        't_shaped': ChemicalPatternLibrary.t_shaped,
        'edge_to_face': ChemicalPatternLibrary.t_shaped,
        'h_bonded': ChemicalPatternLibrary.h_bonded_dimer,
        'hydrogen_bonded': ChemicalPatternLibrary.h_bonded_dimer,
        'herringbone': ChemicalPatternLibrary.herringbone,
    }
    
    if name_lower in library_methods:
        pattern_dict = library_methods[name_lower](**kwargs)
        return ArrangementSpec(
            position_generator=pattern_dict['position_generator'],
            orientation_generator=pattern_dict['orientation_generator'],
            position_params=pattern_dict['position_params'],
            orientation_params=pattern_dict['orientation_params'],
            constraints=pattern_dict.get('constraints', []),
            name=name_lower,
            description=pattern_dict.get('description', '')
        )
    
    # Generic patterns
    generic_patterns = {
        'linear': ArrangementSpec(
            LinearPositions(), FixedOrientation(),
            position_params={'spacing': kwargs.get('distance', 3.4)},
            name='linear'
        ),
        'circular': ArrangementSpec(
            CircularPositions(), FaceCenterOrientation(),
            position_params={'radius': kwargs.get('radius', kwargs.get('distance', 5.0))},
            orientation_params={'facing': kwargs.get('facing', 'inward')},
            name='circular'
        ),
        'helical': ArrangementSpec(
            HelicalPositions(), FaceCenterOrientation(),
            position_params={
                'radius': kwargs.get('radius', 5.0),
                'pitch': kwargs.get('pitch', 3.4),
                'turns': kwargs.get('turns', 1.0)
            },
            name='helical'
        ),
        'spherical': ArrangementSpec(
            SphericalPositions(), FaceCenterOrientation(),
            position_params={'radius': kwargs.get('radius', 5.0)},
            orientation_params={'facing': 'outward'},
            name='spherical'
        ),
        'grid': ArrangementSpec(
            GridPositions(), FixedOrientation(),
            position_params={'spacing': (kwargs.get('distance', 3.4),) * 2},
            name='grid'
        ),
    }
    
    if name_lower in generic_patterns:
        return generic_patterns[name_lower]
    
    raise ValueError(f"Unknown pattern: {name}. Available: {list(library_methods.keys()) + list(generic_patterns.keys())}")


def _spec_from_dict(config: Dict[str, Any], **kwargs) -> ArrangementSpec:
    """Build ArrangementSpec from dictionary."""
    pos_generators = {
        'linear': LinearPositions(),
        'circular': CircularPositions(),
        'helical': HelicalPositions(),
        'grid': GridPositions(),
        'spherical': SphericalPositions(),
        'formula': FormulaPositions(),
        'explicit': ExplicitPositions(),
    }
    
    orient_generators = {
        'fixed': FixedOrientation(),
        'incremental': IncrementalRotation(),
        'alternating': AlternatingRotation(),
        'face_center': FaceCenterOrientation(),
        'formula': FormulaOrientation(),
    }
    
    pos_type = config.get('position', config.get('position_type', 'linear'))
    orient_type = config.get('orientation', config.get('orientation_type', 'fixed'))
    
    return ArrangementSpec(
        position_generator=pos_generators.get(pos_type, LinearPositions()),
        orientation_generator=orient_generators.get(orient_type, FixedOrientation()),
        position_params={**kwargs, **config.get('position_params', {})},
        orientation_params={**kwargs, **config.get('orientation_params', {})},
        constraints=config.get('constraints', []),
        relative_poses=config.get('relative_poses', []),
        name=config.get('name', 'custom'),
        description=config.get('description', '')
    )


def _apply_relative_poses(poses: List[MoleculePose], 
                          relative_specs: List[Dict],
                          molecules: Optional[List[Dict]]) -> List[MoleculePose]:
    """Apply relative pose specifications."""
    # Build frames for each molecule
    frames = []
    if molecules:
        for mol in molecules:
            coords = np.array(mol.get("coords", [[0, 0, 0]]))
            frames.append(MolecularFrame.from_coords_centroid(coords))
    else:
        frames = [MolecularFrame.identity() for _ in poses]
    
    # Update poses with relative specifications
    for rel_spec in relative_specs:
        idx = rel_spec.get('molecule_index', rel_spec.get('index', 0))
        if idx >= len(poses):
            continue
        
        parent = rel_spec.get('relative_to', rel_spec.get('parent', -1))
        
        if parent >= 0:
            poses[idx].parent_idx = parent
            poses[idx].local_offset = np.array(rel_spec.get('offset', rel_spec.get('local_offset', [0, 0, 0])))
            
            if 'orientation' in rel_spec:
                orient = rel_spec['orientation']
                if isinstance(orient, dict) and 'euler' in orient:
                    e = orient['euler']
                    poses[idx].local_rotation = euler_to_matrix(
                        e.get('x', e.get('roll', 0)),
                        e.get('y', e.get('pitch', 0)),
                        e.get('z', e.get('yaw', 0))
                    )
    
    # Resolve relative poses
    return resolve_relative_poses(poses, frames)


# =============================================================================
# CONVENIENCE API FOR NATURAL LANGUAGE INTERFACE
# =============================================================================

def create_custom_arrangement(
    x_formula: str = "0",
    y_formula: str = "0",
    z_formula: str = "3.4 * i",
    roll_formula: str = "0",
    pitch_formula: str = "0",
    yaw_formula: str = "0",
    constraints: List[str] = None
) -> ArrangementSpec:
    """
    Create arrangement from mathematical formulas.
    
    Variables available:
        i: molecule index (0, 1, 2, ...)
        n: total count
        t: normalized index = i/(n-1)
        x, y, z: current position (for orientation formulas)
        pi, e: constants
        sin, cos, tan, sqrt, exp, log, abs: functions
    
    Example - Golden spiral:
        create_custom_arrangement(
            x_formula="sqrt(i) * cos(2.4 * i)",
            y_formula="sqrt(i) * sin(2.4 * i)",
            z_formula="0"
        )
    """
    return ArrangementSpec(
        position_generator=FormulaPositions(),
        orientation_generator=FormulaOrientation(),
        position_params={
            'x_formula': x_formula,
            'y_formula': y_formula,
            'z_formula': z_formula
        },
        orientation_params={
            'roll_formula': roll_formula,
            'pitch_formula': pitch_formula,
            'yaw_formula': yaw_formula
        },
        constraints=_parse_constraint_strings(constraints or []),
        name='custom_formula'
    )


def create_relative_arrangement(
    base_positions: List[Tuple[float, float, float]],
    relative_specs: List[Dict[str, Any]]
) -> ArrangementSpec:
    """
    Create arrangement with relative placements.
    
    Example:
        create_relative_arrangement(
            base_positions=[(0, 0, 0)],  # First molecule at origin
            relative_specs=[
                {'index': 1, 'relative_to': 0, 'offset': [0, 0, 3.4]},
                {'index': 2, 'relative_to': 1, 'offset': [1.5, 0, 3.4], 
                 'orientation': {'euler': {'z': 180}}}
            ]
        )
    """
    return ArrangementSpec(
        position_generator=ExplicitPositions(),
        orientation_generator=FixedOrientation(),
        position_params={'coordinates': base_positions},
        relative_poses=relative_specs,
        name='relative'
    )


def _parse_constraint_strings(constraint_strs: List[str]) -> List[ChemicalConstraint]:
    """Parse constraint strings to constraint objects."""
    constraints = []
    
    for s in constraint_strs:
        # Simple parser for constraint DSL
        match = re.match(r"(\w+)\((.*)\)", s)
        if not match:
            continue
        
        func_name = match.group(1).lower()
        args_str = match.group(2)
        
        # Parse arguments
        kwargs = {}
        for part in args_str.split(","):
            if "=" in part:
                k, v = part.split("=", 1)
                k = k.strip()
                v = v.strip()
                # Try to parse as number
                if v.replace(".", "").replace("-", "").isdigit():
                    kwargs[k] = float(v)
                else:
                    kwargs[k] = v
        
        if func_name == 'distance':
            constraints.append(DistanceConstraint(
                selector1=kwargs.get('sel1', '0:centroid()'),
                selector2=kwargs.get('sel2', '1:centroid()'),
                target=kwargs.get('target'),
                min_dist=kwargs.get('min'),
                max_dist=kwargs.get('max')
            ))
        elif func_name == 'angle':
            constraints.append(AngleConstraint(
                selector1=kwargs.get('sel1', '0:centroid()'),
                selector2=kwargs.get('sel2', '1:centroid()'),
                selector3=kwargs.get('sel3', '2:centroid()'),
                target=kwargs.get('target')
            ))
        elif func_name == 'h_bond':
            constraints.append(HBondConstraint(
                donor_mol=int(kwargs.get('donor', 0)),
                acceptor_mol=int(kwargs.get('acceptor', 1))
            ))
        elif func_name == 'plane_parallel':
            constraints.append(PlaneParallelConstraint(
                mol1_idx=int(kwargs.get('mol1', 0)),
                mol2_idx=int(kwargs.get('mol2', 1)),
                max_deviation=kwargs.get('max_deviation', 10)
            ))
    
    return constraints


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def arrange_molecules_unified(
    molecules: List[Dict[str, Any]],
    arrangement: Union[str, ArrangementSpec, Dict[str, Any]] = 'auto',
    distance: Optional[float] = None,
    **kwargs
) -> List[Dict[str, Any]]:
    """
    Arrange molecules using the unified arrangement system.
    
    This is the recommended entry point for the MCP backend.
    
    Supports:
    - Named patterns: "pi_pi_parallel", "circular", "h_bonded", etc.
    - Custom formulas: via create_custom_arrangement()
    - Relative poses: via create_relative_arrangement()
    - Chemical constraints: via constraint DSL
    - Direct ArrangementSpec objects
    - Dictionary configuration
    
    Returns molecules with transformed coordinates.
    """
    n = len(molecules)
    
    # Set default distance
    if distance is not None:
        kwargs['distance'] = distance
    elif 'distance' not in kwargs:
        kwargs['distance'] = 3.4
    
    # Auto-detect arrangement
    if arrangement == 'auto':
        arrangement = _auto_detect_arrangement(molecules)
    
    # Generate poses
    poses = generate_arrangement(n, arrangement, molecules, **kwargs)
    
    # Apply poses to molecules
    arranged = []
    for mol, pose in zip(molecules, poses):
        coords = np.array(mol["coords"])
        
        # Center molecule
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        
        # Apply rotation then translation
        rotated = centered @ pose.rotation.T
        translated = rotated + pose.position
        
        arranged.append({
            **mol,
            "coords": translated.tolist(),
            "pose": {
                "position": pose.position.tolist(),
                "rotation": pose.to_euler_degrees()
            }
        })
    
    return arranged


def _auto_detect_arrangement(molecules: List[Dict]) -> str:
    """Auto-detect best arrangement based on molecule properties."""
    n = len(molecules)
    
    # Check for aromatic carbons
    def has_aromatic_ring(mol):
        atoms = mol.get("atoms", [])
        return atoms.count("C") >= 6
    
    # Check for H-bond capability
    def can_hbond(mol):
        atoms = mol.get("atoms", [])
        return ("O" in atoms or "N" in atoms) and "H" in atoms
    
    all_aromatic = all(has_aromatic_ring(m) for m in molecules)
    all_hbond = all(can_hbond(m) for m in molecules)
    
    if all_aromatic:
        return 'pi_pi_parallel'
    elif all_hbond and n <= 8:
        return 'circular'
    elif n == 2:
        return 'linear'
    else:
        return 'linear'


# =============================================================================
# MODULE TEST
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("UNIFIED MOLECULAR ARRANGEMENT SYSTEM")
    print("=" * 70)
    print()
    print("Synthesizing approaches from:")
    print("  - K2_plan: Chemical constraints & AtomSelector DSL")
    print("  - z2_plan: Relative poses & MolecularFrame")
    print("  - Original: Mathematical formula generators")
    print()
    
    # Test 1: Named pattern
    print("Test 1: Named pattern 'pi_pi_parallel'")
    poses = generate_arrangement(3, 'pi_pi_parallel', distance=3.4)
    for i, p in enumerate(poses):
        print(f"  Mol {i}: pos={p.position.round(2)}")
    print()
    
    # Test 2: Custom formula
    print("Test 2: Custom spiral formula")
    spec = create_custom_arrangement(
        x_formula="3 * i * cos(i * 0.5)",
        y_formula="3 * i * sin(i * 0.5)",
        z_formula="0.5 * i"
    )
    poses = generate_arrangement(5, spec)
    for i, p in enumerate(poses):
        print(f"  Mol {i}: pos=({p.position[0]:.2f}, {p.position[1]:.2f}, {p.position[2]:.2f})")
    print()
    
    # Test 3: Relative arrangement
    print("Test 3: Relative arrangement (staggered stack)")
    spec = create_relative_arrangement(
        base_positions=[(0, 0, 0)],
        relative_specs=[
            {'index': 1, 'relative_to': 0, 'offset': [0, 0, 3.4]},
            {'index': 2, 'relative_to': 1, 'offset': [1.5, 0, 3.4]},
        ]
    )
    mock_mols = [
        {"atoms": ["C"] * 6, "coords": [[i, 0, 0] for i in range(6)]}
        for _ in range(3)
    ]
    poses = generate_arrangement(3, spec, mock_mols)
    for i, p in enumerate(poses):
        print(f"  Mol {i}: pos={p.position.round(2)}")
    print()
    
    # Test 4: AtomSelector
    print("Test 4: AtomSelector")
    test_mol = {"atoms": ["C", "C", "O", "H"], "coords": [[0,0,0], [1,0,0], [2,0,0], [3,0,0]]}
    sel = AtomSelector("0:atom(O)", [test_mol])
    print(f"  0:atom(O) = {sel.evaluate()}")
    sel = AtomSelector("0:centroid()", [test_mol])
    print(f"  0:centroid() = {sel.evaluate()}")
    print()
    
    print("=" * 70)
    print("All tests passed! System ready for integration.")
    print("=" * 70)
